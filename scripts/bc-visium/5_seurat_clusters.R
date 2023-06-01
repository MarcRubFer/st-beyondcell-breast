rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("clustree")
library("tidygraph")

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
head(seuratobj.clusters@meta.data$seurat_clusters)
# Plot clusters
#  for (r in res) {
#    Idents(seuratobj.clusters) <- paste0("SCT_snn_res.", r)
#    umap.plot <- DimPlot(seuratobj.clusters, group.by = c("ident", "Phase"))
#    spatial.plot <- SpatialDimPlot(seuratobj.clusters, group = "ident")[[1]]
#    spatial.plot.splitted <- SpatialDimPlotSplitted(seuratobj.clusters)
#    ggsave(umap.plot, width = 14, 
#           filename = paste0(out.dir, "/plots/", r, 
#                             "/expression_clusters_umap_phase.png"))
#    ggsave(spatial.plot,
#           filename = paste0(out.dir, "/plots/", r, 
#                             "/expression_clusters_spatial.png"))
#    ggsave(spatial.plot.splitted, width = 17, height = 10, 
#           filename = paste0(out.dir, "/plots/", r,
#                             "/expression_clusters_spatial_splitted.png"))
#    # Save plots
#    all.plots <- list(elbow, umap.plot, spatial.plot, spatial.plot.splitted)
#    save(all.plots, file = paste0(out.dir, "/ggplots/", r, 
#                                  "/seurat_clusters.RData"))
#  }

SpatialFeaturePlot(seuratobj.clusters, features = c("ESR1", "PGR", "ERBB2"))

SpatialFeaturePlot(seuratobj.clusters, features = c("ESR1", "PGR", "ERBB2")) / 
  (umap.plot | spatial.plot)

# Clustree: selection of resolution

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

(clustree.plot | max.stability.plot)
max.stability + max.stability.plot

# Boxplot by cluster, celltypes proportion

cell.types <- names(seuratobj.clusters@meta.data)[6:14]
d <- seuratobj.clusters@meta.data[6:14]
d$cluster <- seuratobj.clusters@meta.data$SCT_snn_res.0.2
d$cluster <- as.factor(d$cluster)

d <- as.data.frame(d)
d.pivot <- pivot_longer(d, cols = -cluster, names_to = "cell.type", values_to = "cell.prop")

ggplot(d.pivot, aes(x = cell.type, y = cell.prop)) +
  geom_boxplot(aes(fill = cell.type)) + 
  facet_wrap(~cluster) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())
# Save data
saveRDS(seuratobj.clusters, file = "./results/analysis/seuratobj.clusters.rds")
