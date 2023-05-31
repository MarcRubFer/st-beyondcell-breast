# Vignettes RCTD
# https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
# https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html

# conda activate spacexr

rm(list = ls())


library("spacexr")
library("Seurat")
library("tidyverse")
library("patchwork")

# Load seurat object
seuratobj.refHER2.filtered <- readRDS(file = "./results/analysis/seuratobj.refHER2.filtered.rds")
dim(seuratobj.refHER2.filtered)
head(seuratobj.refHER2.filtered)

levels(factor(seuratobj.refHER2.filtered$celltype_major))

counts <- seuratobj.refHER2.filtered@assays$RNA@counts
dim(counts)
head(colnames(counts))
head(rownames(counts))

cell.types <- seuratobj.refHER2.filtered@meta.data %>%
  mutate(celltype_major = as.factor(celltype_major)) %>%
  pull(celltype_major) %>%
  setNames(rownames(seuratobj.refHER2.filtered@meta.data))

reference <- Reference(counts, cell.types, n_max_cells = 1500)

## Examine reference object (optional)
dim(reference@counts) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type

names(reference@cell_types)
reference


d <- as.data.frame(reference@cell_types)
colnames(d) <- c("cell_types")
names(d)
patientid <- gsub(pattern = "_.*", replacement = "", names(reference@cell_types))
d$patientid <- patientid

ncells_patients <- ggplot(data = d, aes(x=patientid)) +
  geom_bar(aes(fill = cell_types)) +
  ggtitle(paste("Number of cell types / Distribution by patients")) +
  theme_minimal() +
  theme(axis.title.x = element_blank())

dir.create(path = "./results/plots")
ggsave(filename = "reference.ncells_by_patients.pdf", plot = ncells_patients, path = "./results/plots/")


# Save data
saveRDS(reference, file = "./results/analysis/spacexr.reference.rds")
