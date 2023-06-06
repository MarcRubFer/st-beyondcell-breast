# Vignettes RCTD
# https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
# https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html

# conda activate spacexr

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

rm(list = ls())

# Packages
library("spacexr")
library("Seurat")
library("tidyverse")
library("patchwork")
library("ggsci")

# Load seurat object
seuratobj.refHER2.filtered <- readRDS(file = "./results/analysis/reference/seuratobj.refHER2.filtered.rds")
head(seuratobj.refHER2.filtered)

counts <- seuratobj.refHER2.filtered@assays$RNA@counts
dim(counts)
head(colnames(counts))
head(rownames(counts))

cell.types <- seuratobj.refHER2.filtered@meta.data %>%
  mutate(celltype_major = as.factor(celltype_major)) %>%
  pull(celltype_major) %>%
  setNames(rownames(seuratobj.refHER2.filtered@meta.data))

# Create pre-reference with n_max_cells default 10.000
pre.reference <- Reference(counts, cell.types)
table(pre.reference@cell_types) #number of occurences for each cell type
barplot(table(pre.reference@cell_types), cex.names = 0.8)
rm(pre.reference)

# N_max_cells stablish to 1500 after pre-analysis
reference <- Reference(counts, cell.types, n_max_cells = 1500)

## Examine reference object (optional)
dim(reference@counts) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type


# Plots cell types and patients
df.celltypes <- as.data.frame(reference@cell_types)
colnames(df.celltypes) <- c("cell_types")
names(df.celltypes)
patientid <- gsub(pattern = "_.*", replacement = "", names(reference@cell_types))
df.celltypes$patientid <- patientid

ncells_patients <- ggplot(data = df.celltypes, aes(x=patientid)) +
  geom_bar(aes(fill = cell_types), colour = "#333333") +
  ggtitle(paste("Number of cell types / Distribution by patients")) +
  scale_fill_igv() +
  theme(axis.title.x = element_blank())

ncells_patients_wrap <- ggplot(data = as.data.frame(table(df.celltypes)), 
                               aes(x = patientid, y = Freq)) +
  geom_col(aes(fill = cell_types), 
           position = position_dodge(), 
           color = "black") +
  facet_wrap(~cell_types) +
  scale_fill_igv()

ncells_patients_wrap.2 <- ggplot(data = as.data.frame(table(df.celltypes)), 
       aes(x = cell_types, y = Freq)) +
  geom_col(aes(fill = cell_types), 
           position = position_dodge(), 
           color = "black") +
  facet_wrap(~patientid) +
  scale_fill_igv() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85,0.25))

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/reference"), recursive = TRUE)
ggsave(filename = "reference.ncells_by_patients.png", 
       plot = ncells_patients, 
       path = paste0(out.dir,"/plots/reference"))
ggsave(filename = "reference.ncells_by_patients_wrap.png", 
       plot = ncells_patients_wrap, 
       path = paste0(out.dir,"/plots/reference"))
ggsave(filename = "reference.ncells_by_patients_wrap_2.png", 
       plot = ncells_patients_wrap.2, 
       path = paste0(out.dir,"/plots/reference"))

dir.create(path = paste0(out.dir,"/ggplots/reference"), recursive = TRUE)
all.plots <- list(ncells_patients, ncells_patients_wrap, ncells_patients_wrap.2)
save(all.plots, file = paste0("./results/ggplots/reference/reference.celltypes.RData"))

# Save data
dir.create(path = paste0(out.dir,"/analysis/reference"), recursive = TRUE)
saveRDS(reference, file = "./results/analysis/reference/spacexr.reference.rds")
