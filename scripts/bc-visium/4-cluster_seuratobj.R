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
res <- c(0.1,0.2,0.3,0.4,0.5)
seuratobj.clusters <- FindClusters(seuratobj.clusters, 
                                   resolution = res, 
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
