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
       y = "Median of SC3 stability") +
  ggtitle(label = "SC3 stability at different kparam and resolutions")

# Run the bcUMAP function again, specifying the k.params
res <- c(0.1,0.2,0.3,0.4,0.5)
bc.recomputed <- bcUMAP(bc.recomputed, pc = 20, k.neighbors = 16, res = res)

# Reanalysis of sc3stability for plotting
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

clustree.plot <- clustree.plot + ggtitle(label = "Clustree with 16 k.param")

patch.clustree <- clustersc3.kparam | clustree.plot

# Visualize whether the cells are clustered based on the number of genes detected per each cell.
bc.nFeatureRNA <- bcClusters(bc.recomputed, 
                             UMAP = "beyondcell", 
                             idents = "nFeature_Spatial", 
                             pt.size = 1.5, 
                             factor.col = FALSE) +
  ggtitle("nFeatureRNA")
bc.nCountRNA <- bcClusters(bc.recomputed, 
                           UMAP = "beyondcell", 
                           idents = "nCount_Spatial", 
                           pt.size = 1.5, 
                           factor.col = FALSE) +
  ggtitle("nCountsRNA")
bc.phases <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "Phase", pt.size = 1.5)

bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.2)")
bc.clusters.seurat <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (res 0.2)")

bc.nFeatureRNA | bc.nCountRNA

patch.bc.clusters <- bc.clusters | bc.clusters.seurat

spatial.bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters <- spatial.bc.clusters +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters")

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/beyondcell_create"), recursive = TRUE)

ggsave(filename = "analysis_sc3_kparam.png",
       plot = clustersc3.kparam,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "clustree_beyondcell_k16_res02.png",
       plot = clustree.plot,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "patch_clustree_analysis.png",
       plot = patch.clustree,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "bc_nFeature_dimplot.png",
       plot = bc.nFeatureRNA,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "bc_nCounts_dimplot.png",
       plot = bc.nCountRNA,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "bc_phases_dimplot.png",
       plot = bc.phases,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "beyondcell_clusters_dimplot.png",
       plot = bc.clusters,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "bc_clusters_dimplot_UMAPseurat.png",
       plot = bc.clusters.seurat,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "patch_beyondcell_clusters.png",
       plot = patch.bc.clusters,
       path = paste0(out.dir,"/plots/beyondcell_create/"))
ggsave(filename = "spatial_distr_bc_clusters.png",
       plot = spatial.bc.clusters,
       path = paste0(out.dir,"/plots/beyondcell_create/"))

all.plots <- list(clustersc3.kparam, clustree.plot, patch.clustree, bc.nFeatureRNA ,bc.nCountRNA, bc.phases,
                  bc.clusters, bc.clusters.seurat, patch.bc.clusters, spatial.bc.clusters)
save(all.plots, file = paste0("./results/ggplots/beyondcell_create.RData"))

# Save Data
saveRDS(bc.recomputed, file = "./results/analysis/beyondcellobject.rds")
