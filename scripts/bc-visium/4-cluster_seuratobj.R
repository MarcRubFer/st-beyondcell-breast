rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")

# Load RDS

seuratobj.phase <- readRDS(file = "./results/analysis/seuratobj.phase.rds")

# RunPCA made in previous script
# Elbow plot
elbow <- ElbowPlot(seuratobj.phase, 
                   ndims = 50, 
                   reduction = "pca")
elbow
#ggsave(elbow, filename = paste0(out.dir, "/plots/elbow.pdf")) 

# Clustering
seuratobj.clusters <- FindNeighbors(seuratobj.phase, 
                                    reduction = "pca", 
                                    dims = 1:20, 
                                    k.param = 20)

seuratobj.clusters <- FindClusters(seuratobj.clusters, 
                                   resolution = 0.1, 
                                   verbose = FALSE)

# Run UMAP
seuratobj.clusters <- RunUMAP(seuratobj.clusters, 
                              reduction = "pca", 
                              dims = 1:20, 
                              n.components = 2)

umap.plot <- DimPlot(seuratobj.clusters)
spatial.plot <- SpatialDimPlot(seuratobj.clusters, group = "ident")

(umap.plot | spatial.plot)

SpatialFeaturePlot(seuratobj.clusters, features = c("ESR1", "PGR", "ERBB2"))

SpatialFeaturePlot(seuratobj.clusters, features = c("ESR1", "PGR", "ERBB2")) / 
  (umap.plot | spatial.plot)

# Save data
saveRDS(seuratobj.clusters, file = "./results/analysis/seuratobj.clusters.rds")
