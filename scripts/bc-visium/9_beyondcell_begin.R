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
seuratobj.distances.alt <- readRDS("./results/analysis/seuratobj.distances.alt.rds")

# Set Assay.
DefaultAssay(seuratobj.distances) <- "SCT"
DefaultAssay(seuratobj.distances.alt) <- "SCT"
# Generate geneset object with one of the ready to use signature collections.
gs.ssc <- GetCollection(SSc)

# Compute score for the SSc. This might take a few minutes depending on the size of your dataset.
bc <- bcScore(seuratobj.distances, gs.ssc, expr.thres = 0.1) 
bc.alt <- bcScore(seuratobj.distances.alt, gs.ssc, expr.thres = 0.1) 
# Number of NAs
n.NA <- data.frame(nNAs = colSums(is.na(bc@normalized)),
                   row.names = colnames(bc@normalized))
bc <- bcAddMetadata(bc, n.NA)
bc@meta.data
n.NA.alt <- data.frame(nNAs = colSums(is.na(bc.alt@normalized)),
                       row.names = colnames(bc.alt@normalized))
bc.alt <- bcAddMetadata(bc.alt, n.NA.alt)
# Filter out spots with a high percentage of NAs
bc.filtered <- bcSubset(bc, nan.cells = 0.95)
bc.filtered.alt <- bcSubset(bc.alt, nan.cells = 0.95)

# Replace NAs by 0s
bc.filtered@normalized[is.na(bc.filtered@normalized)] <- 0
bc.recomputed <- bcRecompute(bc.filtered, slot = "normalized")

bc.filtered.alt@normalized[is.na(bc.filtered.alt@normalized)] <- 0
bc.recomputed.alt <- bcRecompute(bc.filtered.alt, slot = "normalized")

#Selection of k parameters
k.param <- c(10, 20, 30, 40, 50)
res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
bcKparRes <- function(bc.object, k.param, res){
  list.data <- lapply(X = k.param, FUN = function(x){
    umap <- bcUMAP(bc.object, pc =20, k.neighbors = x, res = res)
    clustree.graph <- clustree(umap@meta.data, 
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
  
  list.data <- list.data %>%
    bind_rows() %>%
    mutate(k.param = as.factor(k.param))
  return(list.data)
  
}
#l <- lapply(X = k.param, FUN = function(x){
#  res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
#  a <- bcUMAP(bc.recomputed, pc =20, k.neighbors = x, res = res)
#  clustree.graph <- clustree(a@meta.data, 
#                             prefix = "bc_clusters_res.",
#                             prop_filter = 0.1, 
#                             node_colour = "sc3_stability",
#                             return = "graph")
#  
#  max.stability <- clustree.graph %>%
#    activate(nodes) %>%
#    as.data.frame() %>%
#    group_by(bc_clusters_res.) %>%
#    summarise(median.stability = median(sc3_stability),
#              n.clusters = length(unique(cluster))) %>%
#    mutate(k.param = x)
#  return(max.stability)
#})
#
#l <- l %>%
#  bind_rows() %>%
#  mutate(k.param = as.factor(k.param))

sc3.kparam <- bcKparRes(bc.object = bc.recomputed, k.param = k.param, res = res)
sc3.kparam.plot <- ggplot(data=sc3.kparam, aes(x=bc_clusters_res., y=median.stability, group = k.param, color = k.param)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability") +
  ggtitle(label = "SC3 stability at different kparam and resolutions",
          subtitle = "Cell cycle regression: standard") +
  theme_light()

sc3.kparam.alt <- bcKparRes(bc.object = bc.recomputed.alt, k.param = k.param, res = res)
sc3.kparam.alt.plot <- ggplot(data=sc3.kparam.alt, aes(x=bc_clusters_res., y=median.stability, group = k.param, color = k.param)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability") +
  ggtitle(label = "SC3 stability at different kparam and resolutions",
          subtitle = "Cell cycle regression: alternative") +
  theme_light() +
  theme(plot.subtitle = element_text(colour = "red"))

sc3.kparam.plot | sc3.kparam.alt.plot

# Run the bcUMAP function again, specifying the k.params
res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
bc.recomputed <- bcUMAP(bc.recomputed, pc = 20, k.neighbors = 40, res = res)
head(bc.recomputed@meta.data)
bc.recomputed.alt <- bcUMAP(bc.recomputed.alt, pc = 20, k.neighbors = 40, res = res)
head(bc.recomputed.alt@meta.data)
# Reanalysis of sc3stability for plotting
clustree.plot <- clustree(bc.recomputed@meta.data, 
                          prefix = "bc_clusters_res.",
                          node_colour = "sc3_stability") 

clustree.plot.alt <- clustree(bc.recomputed.alt@meta.data, 
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

clustree.plot <- clustree.plot + ggtitle(label = "Clustree with 40 k.param")

patch.clustree <- sc3.kparam.plot | clustree.plot

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

bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.1)")
bc.clusters.seurat <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_res.0.3", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (res 0.3)")

bc.nFeatureRNA | bc.nCountRNA

patch.bc.clusters <- bc.clusters | bc.clusters.seurat

spatial.bc.clusters.01 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.02 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.03 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.3", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))

clustree.plot | (spatial.bc.clusters.01 / spatial.bc.clusters.02 / spatial.bc.clusters.03)

spatial.bc.clusters.01.alt <- bcClusters(bc.recomputed.alt, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.02.alt <- bcClusters(bc.recomputed.alt, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.03.alt <- bcClusters(bc.recomputed.alt, UMAP = "beyondcell", idents = "bc_clusters_res.0.3", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))

clustree.plot.alt | (spatial.bc.clusters.01.alt / spatial.bc.clusters.02.alt / spatial.bc.clusters.03.alt)

spatial.bc.clusters.03 / spatial.bc.clusters.03.alt


spatial.bc.clusters <- spatial.bc.clusters +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters")

# Sankey diagrams?
library(ggsankey)
head(bc.recomputed@meta.data)
df <- bc.recomputed@meta.data %>%
  select(spot.composition.collapse, bc_clusters_res.0.3) %>%
  make_long(spot.composition.collapse,bc_clusters_res.0.3)
pl <- ggplot(df, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl <- pl +geom_sankey(flow.alpha = 0.5
                      , node.color = "black"
                      ,show.legend = T)
pl <- pl +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = 1)
pl + ggtitle("Cell cycle standard")+
  theme_void() + 
  theme(legend.position = "none")


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
saveRDS(bc.recomputed.alt, file = "./results/analysis/beyondcellobject.alt.rds")
