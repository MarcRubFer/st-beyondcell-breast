rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

set.seed(1)
# Read single-cell experiment.
seuratobj.distances <- readRDS("./results/analysis/seuratobj.distances.rds")
bc.object <- readRDS("./results/analysis/beyondcellobject.rds")
head(bc.object@meta.data)

# Boxplot by terapheutc cluster, celltypes proportion

cell.types <- names(bc.object@meta.data)[6:14]
df.celltypes.clusters <- bc.object@meta.data[6:14]
df.celltypes.clusters$cluster <- bc.object@meta.data$bc_clusters_res.0.3
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


