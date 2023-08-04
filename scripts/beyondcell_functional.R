rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

set.seed(1)
# Read single-cell experiment.
seuratobj.distances <- readRDS("./results/analysis/seuratobj.distances.rds")

# Set Assay.
DefaultAssay(seuratobj.distances) <- "SCT"

# Generate geneset with functional gmt.
gs.functional <- GenerateGenesets(x = "./data/gmts/functional.gmt")

# Compute score for the SSc. This might take a few minutes depending on the size of your dataset.
bc.functional <- bcScore(seuratobj.distances, gs.functional, expr.thres = 0.1)

# Number of NAs
n.NA <- data.frame(nNAs = colSums(is.na(bc.functional@normalized)),
                   row.names = colnames(bc.functional@normalized))
bc.functional <- bcAddMetadata(bc.functional, n.NA)
bc.functional@meta.data
# Filter out spots with a high percentage of NAs
bc.func.filtered <- bcSubset(bc.functional, nan.cells = 0.95)

# Replace NAs by 0s
bc.func.filtered@normalized[is.na(bc.func.filtered@normalized)] <- 0
bc.func.recomputed <- bcRecompute(bc.func.filtered, slot = "normalized")

View(bc.func.recomputed@scaled)
func.matrix <- bc.func.recomputed@normalized

pheatmap(mat = func.matrix,
         show_colnames = FALSE)
bcSignatures(bc.func.recomputed, UMAP = "Seurat", signatures = list(values="SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE"), pt.size = 0.5)
bcSignatures(bc.func.recomputed, UMAP = "Seurat", signatures = list(values="BRCANESS"), pt.size = 0.5)

bc4Squares(bc.func.recomputed, idents = "SCT_snn_res.0.3")

# Save Data
saveRDS(bc.func.recomputed, file = "./results/analysis/beyondcell_functional.rds")
