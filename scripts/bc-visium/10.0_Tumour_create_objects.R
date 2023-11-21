rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")

# stablish seed
set.seed(1)

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")
bc.recomputed <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")

# Subset for TME-TCs 
# Extract TC-1 and TC-2 spots
TCs.Tumour <- levels(bc.recomputed@meta.data$TCs_res.0.3)[3:6]
cells.TCs.Tumour <- bc.recomputed@meta.data %>%
  filter(TCs_res.0.3 %in% TCs.Tumour) %>%
  rownames_to_column(var = "spots") %>%
  pull(spots)

# Subset Seurat object (seuratobj.tcs)
seurat.tumour <- subset(seuratobj.tcs, cells = cells.TCs.Tumour)

# After subset, levels of TCs were not update. Re-categorize
seurat.tumour@meta.data$TCs_res.0.3 <- as.character(seurat.tumour@meta.data$TCs_res.0.3)
seurat.tumour@meta.data$TCs_res.0.3 <- as.factor(seurat.tumour@meta.data$TCs_res.0.3)

# Subset Beyondcell object (bc.recomputed)
bc.TCs.Tumour <- bcSubset(bc.recomputed, cells = cells.TCs.Tumour)

# After subset, levels of TCs were not update. Re-categorize
bc.TCs.Tumour@meta.data$TCs_res.0.3 <- as.character(bc.TCs.Tumour@meta.data$TCs_res.0.3)
bc.TCs.Tumour@meta.data$TCs_res.0.3 <- as.factor(bc.TCs.Tumour@meta.data$TCs_res.0.3)


# Save new objects
saveRDS(object = seurat.tumour, 
        file = "./results/analysis/seuratobj.Tumour-TCs.rds")
saveRDS(object = bc.TCs.Tumour,
        file = "./results/analysis/beyondcell.Tumour-TCs.rds")
