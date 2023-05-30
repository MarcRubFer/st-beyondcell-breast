# Vignettes RCTD
# https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
# https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html

# conda activate spacexr

rm(list = ls())


library("spacexr")
library("Seurat")
library("tidyverse")

# Load seurat object
seuratobj.refHER2.filtered <- readRDS(file = "./results/analysis/seuratobj.refHER2.filtered.rds")
dim(seuratobj.refHER2.filtered)
head(seuratobj.refHER2.filtered)

levels(factor(seuratobj.refHER2.filtered$celltype_major))

counts <- seuratobj.refHER2.filtered@assays$RNA@counts
dim(counts)
head(rownames(counts))

cell.types <- seuratobj.refHER2.filtered@meta.data %>%
  mutate(celltype_major = as.factor(celltype_major)) %>%
  pull(celltype_major) %>%
  setNames(rownames(seuratobj.refHER2.filtered@meta.data))

reference <- Reference(counts, cell.types, n_max_cells = 1500)

## Examine reference object (optional)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type



ncells <- seuratobj.refHER2.filtered@meta.data %>%
  count(celltype_major) %>%
  ggplot(aes(x = celltype_major, y = n, fill = celltype_major)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5) +
  ggtitle(paste("Cell types in HER2+ patients")) +
  xlab(NULL) +
  ylab("N") +
  theme_minimal() +
  theme(legend.position = "none")
ncells

ncells.patients <- seuratobj.refHER2.filtered@meta.data %>%
  count(Patient, celltype_major) %>%
  ggplot(aes(x = celltype_major, y = n, fill = Patient)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  ggtitle(paste("Cell types by patients")) +
  xlab(NULL) +
  ylab("N") +
  theme_minimal()
ncells.patients

table(as.data.frame(reference@cell_types))

a <- c("CID4066_GACGTTATCTCGGACG", "CID3921_CCACGGAGTTGTTTGG")
str_split(a, pattern = "_")[[1]]

gsub(pattern = "_*", replacement = "", names(reference@nUMI))
# Save data
saveRDS(reference, file = "./results/analysis/spacexr.reference.rds")
