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
seuratobj.phase.alt <- readRDS(file = "./results/analysis/seuratobj.phase_alternative.rds")

# RunPCA made in previous script
# Elbow plot
elbow <- ElbowPlot(seuratobj.phase, 
                   ndims = 50, 
                   reduction = "pca")
elbow.alt <- ElbowPlot(seuratobj.phase.alt, 
                       ndims = 50, 
                       reduction = "pca")
elbow | elbow.alt
#ggsave(elbow, filename = paste0(out.dir, "/plots/elbow.pdf")) 

# Selection of kNN at different resolution levels

k.param <- c(10,20,30,40,50)
res <- c(0.1,0.2,0.3,0.4,0.5)
# Function to create a data.frame with SC3 stability data
# for each k.parameters and resolution.
selection.kparam.res <- function(object,k.param,res) {
  list.data <- lapply(X = k.param, FUN = function(x){
    res <- res
    obj.temp <- FindNeighbors(object, 
                              reduction = "pca", 
                              dims = 1:20, 
                              k.param = x)
    obj.temp <- FindClusters(obj.temp, 
                             resolution = res, 
                             verbose = FALSE)
    # Run UMAP
    obj.temp <- RunUMAP(obj.temp, 
                        reduction = "pca", 
                        dims = 1:20, 
                        n.components = 2)
    clustree.graph <- clustree(obj.temp, 
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
  
  df.data <- list.data %>%
    bind_rows() %>%
    mutate(k.param = as.factor(k.param))
  return(df.data)
}

sc3.phase.data <- selection.kparam.res(object = seuratobj.phase, k.param = k.param, res = res)
sc3.phase.graph <- ggplot(data=sc3.phase.data, aes(x=SCT_snn_res., y=median.stability, group = k.param, color = k.param)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability") +
  ylim(0.0,0.8) +
  ggtitle(label = "SC3 stability across resolution and k.params", subtitle = "Cell cycle regression: standard") +
  theme_light()

sc3.phase.alt.data <- selection.kparam.res(object = seuratobj.phase.alt, k.param = k.param, res = res)
sc3.phase.alt.graph <- ggplot(data=sc3.phase.alt.data, aes(x=SCT_snn_res., y=median.stability, group = k.param, color = k.param)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability") +
  ylim(0.0,0.8) +
  ggtitle(label = "SC3 stability across resolution and k.params", subtitle = "Cell cycle regression: alternative") +
  theme_light() +
  theme(plot.subtitle = element_text(colour = "red"))

sc3.phase.graph | sc3.phase.alt.graph


# Recompute Clustering at selected K.param
# (in this case k.param = 50)
seuratobj.clusters <- FindNeighbors(seuratobj.phase, 
                                    reduction = "pca", 
                                    dims = 1:20, 
                                    k.param = 50)
seuratobj.clusters.alt <- FindNeighbors(seuratobj.phase.alt, 
                                        reduction = "pca", 
                                        dims = 1:20, 
                                        k.param = 50)
res <- c(0.1,0.2,0.3,0.4,0.5)
seuratobj.clusters <- FindClusters(seuratobj.clusters, 
                                   resolution = res, 
                                   verbose = FALSE)
seuratobj.clusters.alt <- FindClusters(seuratobj.clusters.alt, 
                                   resolution = res, 
                                   verbose = FALSE)

# ReRun UMAP
seuratobj.clusters <- RunUMAP(seuratobj.clusters, 
                              reduction = "pca", 
                              dims = 1:20, 
                              n.components = 2)
seuratobj.clusters.alt <- RunUMAP(seuratobj.clusters.alt, 
                              reduction = "pca", 
                              dims = 1:20, 
                              n.components = 2)

# Clustree plots
## Phase
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

## Alternative phase
clustree.plot.alt <- clustree(seuratobj.clusters.alt, 
                          prefix = "SCT_snn_res.",
                          node_colour = "sc3_stability") 

clustree.analysis <- (clustree.plot | clustree.plot.alt)

# Boxplot by cluster, celltypes proportion
df.celltypes.clusters <- seuratobj.clusters@meta.data %>%
  select(B.cells:T.cells,SCT_snn_res.0.2) %>%
  rename(clusters = SCT_snn_res.0.2) %>%
  mutate(clusters = as.factor(clusters)) %>%
  pivot_longer(cols = -clusters,
               names_to = "cell.type",
               values_to = "cell.prop")

boxplot.celltypes.clusters <- ggplot(df.celltypes.clusters, aes(x = cell.type, y = cell.prop)) +
  geom_boxplot(aes(fill = cell.type)) + 
  facet_wrap(~clusters) +
  ggtitle(label = "Cell types proportions within each cluster", subtitle = "Cell cycle regression: standard") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) 

spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2")
dim.clusters <- DimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2")

df.celltypes.clusters.alt <- seuratobj.clusters.alt@meta.data %>%
  select(B.cells:T.cells,SCT_snn_res.0.2) %>%
  rename(clusters = SCT_snn_res.0.2) %>%
  mutate(clusters = as.factor(clusters)) %>%
  pivot_longer(cols = -clusters,
               names_to = "cell.type",
               values_to = "cell.prop")

boxplot.celltypes.clusters.alt <- ggplot(df.celltypes.clusters.alt, aes(x = cell.type, y = cell.prop)) +
  geom_boxplot(aes(fill = cell.type)) + 
  facet_wrap(~clusters) +
  ggtitle(label = "Cell types proportions within each cluster", subtitle = "Cell cycle regression: alternative") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        plot.subtitle = element_text(colour = "red")) 

spatial.clusters.alt <- SpatialDimPlot(seuratobj.clusters.alt, group.by = "SCT_snn_res.0.2")
dim.clusters.alt <- DimPlot(seuratobj.clusters.alt, group.by = "SCT_snn_res.0.2")



boxplot.celltypes.clusters | boxplot.celltypes.clusters.alt
spatial.clusters / spatial.clusters.alt
dim.clusters | dim.clusters.alt

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
saveRDS(seuratobj.clusters.alt, file = "./results/analysis/seuratobj.clusters.alt.rds")






