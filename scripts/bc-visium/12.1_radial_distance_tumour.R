seurat.tumor <- readRDS("./results/analysis/seuratobj.Tumour-TCs.rds")

seurat.tumor@meta.data

SpatialDimPlot(seurat.tumor, group.by = "spot.dual") / SpatialDimPlot(seurat.tumor, group.by = "TCs_res.0.3")
