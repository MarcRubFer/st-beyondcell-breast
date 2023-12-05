rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE, showWarnings = FALSE)

cell.colors <- c("B-cells" = "cornflowerblue",
                 "CAFs" = "goldenrod2",
                 "Cancer Epithelial" = "red3",
                 "Endothelial" = "seagreen4",
                 "Myeloid" = "darkorchid",
                 "Normal Epithelial" = "hotpink",
                 "Plasmablasts" = "greenyellow",
                 "PVL" = "darkorange1",
                 "T-cells" = "steelblue4")

cell.cycle.colors <- c("G1" = "tomato",
                       "G2M" = "royalblue",
                       "S" = "green3")

## CODE ##
set.seed(1)

# Load seurat reference HER2 patients
seuratobj.refHER2.filtered <- readRDS("./results/analysis/reference/seuratobj.refHER2.filtered.rds")
head(seuratobj.refHER2.filtered)

seurat.refHER2.SCT <- SCTransform(seuratobj.refHER2.filtered, 
                        assay = "RNA", 
                        return.only.var.genes = TRUE, 
                        verbose = TRUE,
                        min_cells = 100,
                        vst.flavor = "v1")
seurat.refHER2.SCT <- RunPCA(seurat.refHER2.SCT, 
                   assay = "SCT", 
                   npcs = 50, 
                   features = VariableFeatures(seurat.refHER2.SCT))

dimplot <- DimPlot(seurat.refHER2.SCT, group.by = "celltype_major", cols = cell.colors)

# Cell cycle effect
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes

seurat.refHER2.phase <- CellCycleScoring(seurat.refHER2.SCT,
                              g2m.features = g2m.genes, 
                              s.features = s.genes, 
                              assay = "SCT")
seurat.refHER2.phase <- RunPCA(seurat.refHER2.phase, 
                    assay = "SCT", 
                    npcs = 50, 
                    features = VariableFeatures(seurat.refHER2.phase))

Idents(seurat.refHER2.phase) <- "Phase"
cell.cycle.without <- DimPlot(seurat.refHER2.phase, cols = cell.cycle.colors) + ggtitle("Non-regressed cell cycle")
cell.cycle.phase.without <- DimPlot(seurat.refHER2.phase, split.by = "Phase", cols = cell.cycle.colors)
patch.cell.cycle.without <- cell.cycle.without + cell.cycle.phase.without

elbow <- ElbowPlot(seurat.refHER2.phase, 
                   ndims = 50, 
                   reduction = "pca")


seurat.refHER2.clusters <- FindNeighbors(seurat.refHER2.phase, 
                                    reduction = "pca", 
                                    dims = 1:20, 
                                    k.param = 20)
res <- c(0.1,0.2,0.3,0.4,0.5)
seurat.refHER2.clusters <- FindClusters(seurat.refHER2.clusters, 
                                   resolution = res, 
                                   verbose = FALSE)

# Run UMAP
seurat.refHER2.clusters <- RunUMAP(seurat.refHER2.clusters, 
                              reduction = "pca", 
                              dims = 1:20, 
                              n.components = 2)

# Clustree plots
clustree.plot <- clustree(seurat.refHER2.clusters, 
                          prefix = "SCT_snn_res.",
                          node_colour = "sc3_stability") 

UMAP <- DimPlot(seurat.refHER2.clusters, reduction = "umap", group.by = "celltype_major", cols = cell.colors) + NoLegend()
LabelClusters(plot = UMAP, id = "celltype_major", box = TRUE)
