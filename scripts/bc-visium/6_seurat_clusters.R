rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("clustree")
library("tidygraph")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Load RDS

seuratobj.phase <- readRDS(file = "./results/analysis/seuratobj.phase.rds")

# RunPCA made in previous script
# Elbow plot
elbow <- ElbowPlot(seuratobj.phase, 
                   ndims = 50, 
                   reduction = "pca")
elbow
#ggsave(elbow, filename = paste0(out.dir, "/plots/elbow.pdf")) 

# Selection of kNN at different resolution levels

k.param <- c(10,20,30,40,50)
l <- lapply(X = k.param, FUN = function(x){
  res <- c(0.1,0.2,0.3,0.4,0.5)
  a <- FindNeighbors(seuratobj.phase, 
                     reduction = "pca", 
                     dims = 1:20, 
                     k.param = x)
  a <- FindClusters(a, 
                    resolution = res, 
                    verbose = FALSE)
  # Run UMAP
  a <- RunUMAP(a, 
               reduction = "pca", 
               dims = 1:20, 
               n.components = 2)
  clustree.graph <- clustree(a, 
                             prefix = "SCT_snn_res.",
                             prop_filter = 0.1, 
                             node_colour = "sc3_stability",
                             return = "graph")
  max.stability <- clustree.graph %>%
    activate(nodes) %>%
    as.data.frame() %>%
    group_by(SCT_snn_res.) %>%
    summarise(median.stability = median(sc3_stability),
              n.clusters = length(unique(cluster))) %>%
    mutate(k.param = x)
  return(max.stability)
})

l <- l %>%
  bind_rows() %>%
  mutate(k.param = as.factor(k.param))

clustersc3.kparam <- ggplot(data=l, aes(x=SCT_snn_res., y=median.stability, group = k.param, color = k.param)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability")


# Recompute Clustering at selected K.param
# (in this case k.param = 50)
seuratobj.clusters <- FindNeighbors(seuratobj.phase, 
                                    reduction = "pca", 
                                    dims = 1:20, 
                                    k.param = 20)
res <- c(0.1,0.2,0.3,0.4,0.5)
seuratobj.clusters <- FindClusters(seuratobj.clusters, 
                                   resolution = res, 
                                   verbose = FALSE)

# ReRun UMAP
seuratobj.clusters <- RunUMAP(seuratobj.clusters, 
                              reduction = "pca", 
                              dims = 1:20, 
                              n.components = 2)
  
# Clustree plots
clustree.plot <- clustree(seuratobj.clusters, 
                          prefix = "SCT_snn_res.",
                          node_colour = "sc3_stability") 

clustree.graph <- clustree(seuratobj.clusters, 
                           prefix = "SCT_snn_res.",
                           prop_filter = 0.1, 
                           node_colour = "sc3_stability",
                           return = "graph")

max.stability <- clustree.graph %>%
  activate(nodes) %>%
  as.data.frame() %>%
  group_by(SCT_snn_res.) %>%
  summarise(median.stability = median(sc3_stability),
            n.clusters = length(unique(cluster)))

max.stability.plot <- ggplot(data=max.stability, aes(x=SCT_snn_res., y=median.stability, group = 1)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability")

clustree.analysis <- (clustree.plot | max.stability.plot)

# Boxplot by cluster, celltypes proportion

cell.types <- names(seuratobj.clusters@meta.data)[6:14]
df.celltypes.clusters <- seuratobj.clusters@meta.data[6:14]
df.celltypes.clusters$cluster <- seuratobj.clusters@meta.data$SCT_snn_res.0.2
df.celltypes.clusters$cluster <- seuratobj.clusters@meta.data$SCT_snn_res.0.1
df.celltypes.clusters$cluster <- as.factor(df.celltypes.clusters$cluster)

df.celltypes.clusters <- as.data.frame(df.celltypes.clusters)
df.celltypes.clusters.pivot <- pivot_longer(df.celltypes.clusters, 
                                             cols = -cluster, 
                                             names_to = "cell.type", 
                                             values_to = "cell.prop")

boxplot.celltypes.clusters <- ggplot(df.celltypes.clusters.pivot, aes(x = cell.type, y = cell.prop)) +
  geom_boxplot(aes(fill = cell.type)) + 
  facet_wrap(~cluster) +
  ggtitle(label = "Cell types proportions within each cluster") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) 

spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2")
dim.clusters <- DimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2")

spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1")
dim.clusters <- DimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1")
boxplot.celltypes.clusters / (spatial.clusters | dim.clusters)

breastcancermarkers.spatial <- SpatialFeaturePlot(seuratobj.clusters, features = c("ESR1", "PGR", "ERBB2"), ncol = 4)
breastcancermarkers.spatial / (spatial.clusters | dim.clusters)

DimPlot(seuratobj.clusters, group.by = "orig.ident")


# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/clustering/"), recursive = TRUE)
ggsave(filename = "SC3stability_by_kparams.png",
       plot = clustersc3.kparam,
       path = paste0(out.dir,"/plots/clustering/"))
ggsave(filename = "clustree_analysis.png",
       plot = clustree.analysis,
       path = paste0(out.dir,"/plots/clustering/"))
ggsave(filename = "boxplot_celltypes_by_clusters.png",
       plot = boxplot.celltypes.clusters,
       path = paste0(out.dir,"/plots/clustering/"))
ggsave(filename = "clusters_on_slides.png",
       plot = spatial.clusters,
       path = paste0(out.dir,"/plots/clustering/"))
ggsave(filename = "dimplot_clusters.png",
       plot = dim.clusters,
       path = paste0(out.dir,"/plots/clustering/"))
ggsave(filename = "spatial_clusters.png",
       plot = spatial.clusters,
       path = paste0(out.dir,"/plots/clustering/"))
ggsave(filename = "breastcancermarkers_on_slides.png",
       plot = breastcancermarkers.spatial,
       path = paste0(out.dir,"/plots/clustering/"))

all.plots <- list(clustersc3.kparam, clustree.analysis, boxplot.celltypes.clusters,spatial.clusters, dim.clusters, breastcancermarkers.spatial)
save(all.plots, file = paste0("./results/ggplots/clustering.RData"))

# Save data
saveRDS(seuratobj.clusters, file = "./results/analysis/seuratobj.clusters.rds")







