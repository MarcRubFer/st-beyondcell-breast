# Create an integrated object (seurat v.4.3)
# vignette: https://satijalab.org/seurat/archive/v4.3/integration_rpca

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

# Split the dataset into a list of length = n.patients
patient.list <- SplitObject(seuratobj.refHER2.filtered, split.by = "orig.ident")

# Normalize each object with SCTransform
patient.list.norm <- lapply(X = patient.list, FUN = SCTransform, method = "glmGamPoi")

# Select features for integration
features <- SelectIntegrationFeatures(object.list = patient.list.norm, nfeatures = 3000)

# Prepare objects for integration and run PCAs
patient.list.prep <- PrepSCTIntegration(object.list = patient.list.norm, anchor.features = features)
patient.list.prep <- lapply(X = patient.list.prep, FUN = RunPCA, features = features)

# Search anchors for integration process
patient.anchors <- FindIntegrationAnchors(object.list = patient.list.prep, 
                                          normalization.method = "SCT",
                                          anchor.features = features, 
                                          dims = 1:30, 
                                          reduction = "rpca", 
                                          k.anchor = 20)
# Integrate data
seurat.refHER2.integrated <- IntegrateData(anchorset = patient.anchors, 
                                           normalization.method = "SCT", 
                                           dims = 1:30)

# Run PCA and UMAP
seurat.refHER2.integrated <- RunPCA(seurat.refHER2.integrated, verbose = FALSE)
seurat.refHER2.integrated <- RunUMAP(seurat.refHER2.integrated, reduction = "pca", dims = 1:30)

cell.colors <- c("B-cells" = "cornflowerblue",
                 "CAFs" = "goldenrod2",
                 "Cancer Epithelial" = "red3",
                 "Endothelial" = "seagreen4",
                 "Myeloid" = "darkorchid",
                 "Normal Epithelial" = "hotpink",
                 "Plasmablasts" = "greenyellow",
                 "PVL" = "darkorange1",
                 "T-cells" = "steelblue4")

dimplot.patients.integrated <- DimPlot(seurat.refHER2.integrated, reduction = "umap", group.by = "orig.ident" )
dimplot.celltype.integrated <- DimPlot(seurat.refHER2.integrated, reduction = "umap", group.by = "celltype_major", cols = cell.colors, label = T, label.box = T)

all.plots.integrated <- list(dimplot.celltype.integrated, dimplot.patients.integrated)
save(all.plots.integrated, file = paste0("./results/ggplots/reference/integrated_plots.RData"))
saveRDS(seurat.refHER2.integrated, file = "./results/analysis/reference/seurat.refHER2.integrated.rds")
