rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")

set.seed(1)
# Read single-cell experiment.
seuratobj.clusters <- readRDS("./results/analysis/seuratobj.clusters.rds")

# Set Assay.
DefaultAssay(seuratobj.clusters) <- "SCT"

# Generate geneset object with one of the ready to use signature collections.
gs.ssc <- GetCollection(SSc)

# Compute score for the PSc. This might take a few minutes depending on the size of your dataset.
bc <- bcScore(seuratobj.clusters, gs.ssc, expr.thres = 0.1) 

# Number of NAs
n.NA <- data.frame(nNAs = colSums(is.na(bc@normalized)),
                   row.names = colnames(bc@normalized))
bc <- bcAddMetadata(bc, n.NA)

# Filter out spots with a high percentage of NAs
bc.filtered <- bcSubset(bc, nan.cells = 0.95)

# Replace NAs by 0s
bc.filtered@normalized[is.na(bc.filtered@normalized)] <- 0
bc.recomputed <- bcRecompute(bc.filtered, slot = "normalized")

# Run the UMAP reduction. 
res <- c(0.1,0.2,0.3,0.4,0.5)
bc.recomputed <- bcUMAP(bc.recomputed, k.neighbors = 4, res = res)

#Selection of k parameters
k.param <- c(4,8,12,16,20,40,60)
l <- lapply(X = k.param, FUN = function(x){
  res <- c(0.1,0.2,0.3,0.4,0.5)
  a <- bcUMAP(bc.recomputed, pc =20, k.neighbors = x, res = res)
  clustree.graph <- clustree(a@meta.data, 
                             prefix = "bc_clusters_res.",
                             prop_filter = 0.1, 
                             node_colour = "sc3_stability",
                             return = "graph")
  
  max.stability <- clustree.graph %>%
    activate(nodes) %>%
    as.data.frame() %>%
    group_by(bc_clusters_res.) %>%
    summarise(median.stability = median(sc3_stability),
              n.clusters = length(unique(cluster))) %>%
    mutate(k.param = x)
  return(max.stability)
})

l <- l %>%
  bind_rows() %>%
  mutate(k.param = as.factor(k.param))

clustersc3.kparam <- ggplot(data=l, aes(x=bc_clusters_res., y=median.stability, group = k.param, color = k.param)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability")

# Run the bcUMAP function again, specifying the number of principal components you want to use.
bc.recomputed <- bcUMAP(bc.recomputed, pc = 20, k.neighbors = 4, res = res)


clustree.plot <- clustree(bc.recomputed@meta.data, 
                          prefix = "bc_clusters_res.",
                          node_colour = "sc3_stability") 

clustree.graph <- clustree(bc.recomputed@meta.data, 
                           prefix = "bc_clusters_res.",
                           prop_filter = 0.1, 
                           node_colour = "sc3_stability",
                           return = "graph")

max.stability <- clustree.graph %>%
  activate(nodes) %>%
  as.data.frame() %>%
  group_by(bc_clusters_res.) %>%
  summarise(median.stability = median(sc3_stability),
            n.clusters = length(unique(cluster)))

max.stability.plot <- ggplot(data=max.stability, aes(x=bc_clusters_res., y=median.stability, group = 1)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability")

(clustree.plot | max.stability.plot)

# Visualize whether the cells are clustered based on the number of genes detected per each cell.
bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "nFeature_Spatial", pt.size = 1.5, factor.col = FALSE)
bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "nCount_Spatial", pt.size = 1.5, factor.col = FALSE)
bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "Phase", pt.size = 1.5)

a <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.4", pt.size = 1.5) +
  ggtitle("BeyondCell Clusters (res 0.3)")
b <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_res.0.4", pt.size = 1.5) +
  ggtitle("Seurat Clusters (res 0.2)")

a | b

bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.3", pt.size = 1.5, spatial = TRUE)

saveRDS(bc, file = "./results/analysis/beyondcellobject.rds")
