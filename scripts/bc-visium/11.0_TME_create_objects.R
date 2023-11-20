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
TCs.TME <- levels(bc.recomputed@meta.data$TCs_res.0.3)[1:2]
cells.TCs.TME <- bc.recomputed@meta.data %>%
  filter(TCs_res.0.3 %in% TCs.TME) %>%
  rownames_to_column(var = "spots") %>%
  pull(spots)

# Subset Seurat object (seuratobj.tcs)
seurat.TME <- subset(seuratobj.tcs, cells = cells.TCs.TME)

# Subset Beyondcell object (bc.recomputed)
bc.TCs.TME <- bcSubset(bc.recomputed, cells = cells.TCs.TME)

# After subset, levels of TCs were not update. Re-categorize
bc.TCs.TME@meta.data$TCs_res.0.3 <- as.character(bc.TCs.TME@meta.data$TCs_res.0.3)
bc.TCs.TME@meta.data$TCs_res.0.3 <- as.factor(bc.TCs.TME@meta.data$TCs_res.0.3)


# Save new objects
saveRDS(object = seurat.TME, 
        file = "./results/analysis/seuratobj.TME-TCs.rds")
saveRDS(object = bc.TCs.TME,
        file = "./results/analysis/beyondcell.TME-TCs.rds")
