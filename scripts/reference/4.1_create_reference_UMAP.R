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
                        vst.flavor = "v1")
seurat.refHER2.SCT <- RunPCA(seurat.refHER2.SCT, 
                   assay = "SCT", 
                   npcs = 50, 
                   features = VariableFeatures(seurat.refHER2.SCT))

head(seurat.refHER2.SCT)
dimplot <- DimPlot(seurat.refHER2.SCT, group.by = "celltype_major", cols = cell.colors)
DimPlot(seurat.refHER2.SCT, group.by = "orig.ident") | DimPlot(seurat.refHER2.patient, group.by = "orig.ident") 
# Regression by patient (orig.ident)
seurat.refHER2.patient <- SCTransform(seurat.refHER2.SCT, 
                                      assay = "RNA", 
                                      vars.to.regress = "orig.ident",
                                      return.only.var.genes = TRUE, 
                                      verbose = TRUE)
seurat.refHER2.patient <- RunPCA(seurat.refHER2.patient, 
                             assay = "SCT", 
                             npcs = 50, 
                             features = VariableFeatures(seurat.refHER2.patient))
dimplot.patient <- DimPlot(seurat.refHER2.patient, group.by = "celltype_major", cols = cell.colors) +
  ggtitle("patient regressed")

dimplot | dimplot.patient

# Cell cycle effect
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes

seurat.refHER2.pre.phase <- CellCycleScoring(seurat.refHER2.patient,
                              g2m.features = g2m.genes, 
                              s.features = s.genes, 
                              assay = "SCT")
seurat.refHER2.pre.phase <- RunPCA(seurat.refHER2.pre.phase, 
                    assay = "SCT", 
                    npcs = 50, 
                    features = VariableFeatures(seurat.refHER2.pre.phase))

Idents(seurat.refHER2.pre.phase) <- "Phase"
cell.cycle.without <- DimPlot(seurat.refHER2.pre.phase, cols = cell.cycle.colors) + ggtitle("Non-regressed cell cycle")
cell.cycle.phase.without <- DimPlot(seurat.refHER2.pre.phase, split.by = "Phase", cols = cell.cycle.colors)
patch.cell.cycle.without <- cell.cycle.without + cell.cycle.phase.without

# Cell cycle regression
seurat.refHER2.phase <- SCTransform(seurat.refHER2.pre.phase, 
                                      assay = "RNA", 
                                      vars.to.regress = c("S.Score", "G2M.Score"),
                                      return.only.var.genes = TRUE, 
                                      verbose = TRUE)
seurat.refHER2.phase <- RunPCA(seurat.refHER2.phase, 
                                   assay = "SCT", 
                                   npcs = 50, 
                                   features = VariableFeatures(seurat.refHER2.phase))

Idents(seurat.refHER2.phase) <- "Phase"
cell.cycle.with <- DimPlot(seurat.refHER2.phase, cols = cell.cycle.colors) + ggtitle("Regressed cell cycle")
cell.cycle.phase.with <- DimPlot(seurat.refHER2.phase, split.by = "Phase", cols = cell.cycle.colors)
patch.cell.cycle.with <- cell.cycle.with + cell.cycle.phase.with


patch.cell.cycle.without | patch.cell.cycle.with


seurat.refHER2.UMAP <- RunUMAP(seurat.refHER2.phase, 
                                   reduction = "pca", 
                                   dims = 1:20, 
                                   n.components = 2)

UMAP.celltype <- DimPlot(seurat.refHER2.UMAP, reduction = "umap", group.by = "celltype_major", cols = cell.colors, label = T, label.box = T) + NoLegend() + ggtitle("Cell type")
UMAP.patient <- DimPlot(seurat.refHER2.UMAP, reduction = "umap", group.by = "orig.ident", label = T, label.box = T) + NoLegend() + ggtitle("Patient")

UMAP.patch <- UMAP.celltype | UMAP.patient
UMAP.patch

ggsave(filename = "UMAP_reference_celltypes.png",
       plot = UMAP.celltype,
       path = "./results/plots/reference/",
       scale = 1,
       width = 12,
       height = 12,
       units = "in")
