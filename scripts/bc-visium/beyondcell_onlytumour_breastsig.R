rm(list = ls())

#devtools::install_github("igrabski/sc-SHC")
library("scSHC")

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("viridis")
library("RColorBrewer")


out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# stablish seed
set.seed(1)

# Read seurat object.
seuratobj.distances <- readRDS("./results/analysis/seuratobj.distances.rds")

# Reformat spot composition in two levels: Pure Tumour and Rest

## REVIEW: is it needed to reformat? Direct subset of Pure_Tumour has same dimensions ##
spot.dual <- seuratobj.distances@meta.data %>%
  select(spot.composition.filter) %>%
  mutate(spot.dual = case_when(spot.composition.filter == "Pure_Tumour" ~ "Pure_Tumour",
                               TRUE ~ "Rest"),
         spot.composition.filter = NULL)

seuratobj.distances <- AddMetaData(seuratobj.distances, metadata = spot.dual)
seuratobj.distances@meta.data

seuratobj.pure <- subset(seuratobj.distances, subset = spot.dual == "Pure_Tumour")
seuratobj.pure@meta.data

pure <- subset(seuratobj.distances, subset = (spot.composition.filter == "Pure_Tumour"))
pure@meta.data

# Generate new geneset with only breast signatures
gs.breast <- GenerateGenesets(x = "./data/gmts/drug_signatures_classic_nodup.gmt")

# Compute beyondcell score with breast signatures
bc <- bcScore(seuratobj.pure, gs.breast, expr.thres = 0.1) 

# Number of NAs
n.NA <- data.frame(nNAs = colSums(is.na(bc@normalized)),
                   row.names = colnames(bc@normalized))
bc <- bcAddMetadata(bc, n.NA)
bc@meta.data

# Filter out spots with a high percentage of NAs
bc.filtered <- bcSubset(bc, nan.cells = 0.95)

# Replace NAs by 0s
bc.filtered@normalized[is.na(bc.filtered@normalized)] <- 0
bc.recomputed <- bcRecompute(bc.filtered, slot = "normalized")

#Selection of k parameters
k.param <- c(10, 20, 30, 40, 50)
l <- lapply(X = k.param, FUN = function(x){
  res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
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
res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
bc.recomputed <- bcUMAP(bc.recomputed, pc = 20, k.neighbors = 30, res = res)
head(bc.recomputed@meta.data)

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

clustree.plot <- clustree.plot + ggtitle(label = "Clustree with 30 k.param")



# Test significance of clusters (sc-SHC package)

# Round scaled data
scaled <- round(bc.recomputed@scaled * 100, 0)

## RESOLUTION 0.1
# Compute Therapeutic Clusters (TCs)
original.clusters.res01 <- bc.recomputed@meta.data[colnames(scaled), ] %>%
  pull(all_of("bc_clusters_res.0.1")) %>%
  as.character()
new.clusters.01 <- testClusters(scaled, 
                             cluster_ids = original.clusters.res01,
                             batch = NULL, 
                             alpha = 0.05, 
                             num_features = 100, 
                             num_PCs = 20, 
                             parallel = FALSE)
# Add metadata
TCs <- as.data.frame(new.clusters.01[[1]])
colnames(TCs) <- "bc_clusters_new_res_0.1"
rownames(TCs) <- colnames(scaled)
TCs <- TCs %>%
  mutate(bc_clusters_new_res_0.1 = str_remove(bc_clusters_new_res_0.1, pattern = "^new"),
         bc_clusters_new_renamed_res_0.1 = case_when(bc_clusters_new_res_0.1 == "2" ~ "1",
                                             bc_clusters_new_res_0.1 == "3" ~ "2",
                                             bc_clusters_new_res_0.1 == "1" ~ "3",
                                             TRUE ~ bc_clusters_new_res_0.1),
         bc_clusters_new_res_0.1 = as.factor(bc_clusters_new_res_0.1),
         bc_clusters_new_renamed_res_0.1 = as.factor(bc_clusters_new_renamed_res_0.1))
print(head(TCs))


bc.recomputed <- bcAddMetadata(bc.recomputed, TCs)
head(bc.recomputed@meta.data)

# Plot UMAPs Beyondcell and Seurat clustering
bc.clusters.res.0.1 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.1)")
bc.clusters.seurat.res.0.1 <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_res.0.1", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (res 0.1)")
patch.bc.clusters.res.0.1 <- bc.clusters.res.0.1 | bc.clusters.seurat.res.0.1

bc.clusters.sc.res.0.1 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.1", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (sc-SHC)")
bc.clusters.seurat.sc.res.0.1 <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_new_renamed_res_0.1", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (sc-SHC)")

patch.bc.clusters.sc <- bc.clusters.sc.res.0.1 | bc.clusters.seurat.sc.res.0.1

patch.clustering <- bc.clusters.res.0.1 | bc.clusters.sc.res.0.1



# Plot spatial distribution
spatial.bc.clusters.res.0.1 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.res.0.1 <- spatial.bc.clusters.res.0.1 +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters")

spatial.bc.clusters.sc.res.0.1 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.1", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.sc.res.0.1 <- spatial.bc.clusters.sc.res.0.1 +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters (sc-SHC)")

patch.spatial.clustering <- spatial.bc.clusters.res.0.1 / spatial.bc.clusters.sc.res.0.1



## RESOLUTION 0.2
# Compute Therapeutic Clusters (TCs)
original.clusters.res02 <- bc.recomputed@meta.data[colnames(scaled), ] %>%
  pull(all_of("bc_clusters_res.0.2")) %>%
  as.character()
new.clusters.02 <- testClusters(scaled, 
                                cluster_ids = original.clusters.res02,
                                batch = NULL, 
                                alpha = 0.05, 
                                num_features = 100, 
                                num_PCs = 20, 
                                parallel = FALSE)
# Add metadata
TCs <- as.data.frame(new.clusters.02[[1]])
colnames(TCs) <- "bc_clusters_new_res_0.2"
rownames(TCs) <- colnames(scaled)
TCs <- TCs %>%
  mutate(bc_clusters_new_res_0.2 = str_remove(bc_clusters_new_res_0.2, pattern = "^new"),
         bc_clusters_new_renamed_res_0.2 = case_when(bc_clusters_new_res_0.2 == "3" ~ "1",
                                                     bc_clusters_new_res_0.2 == "1" ~ "2",
                                                     bc_clusters_new_res_0.2 == "2" ~ "3",
                                                     TRUE ~ bc_clusters_new_res_0.2),
         bc_clusters_new_res_0.2 = as.factor(bc_clusters_new_res_0.2),
         bc_clusters_new_renamed_res_0.2 = as.factor(bc_clusters_new_renamed_res_0.2))
print(head(TCs))


bc.recomputed <- bcAddMetadata(bc.recomputed, TCs)
head(bc.recomputed@meta.data)

# Plot UMAPs Beyondcell and Seurat clustering
bc.clusters.res.0.2 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.2)")
bc.clusters.seurat.res.0.2 <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (res 0.2)")

patch.bc.clusters.res.0.2 <- bc.clusters.res.0.2 | bc.clusters.seurat.res.0.2

bc.clusters.sc.res.0.2 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (sc-SHC)")
bc.clusters.seurat.sc.res.0.2 <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_new_renamed_res_0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (sc-SHC)")

patch.bc.clusters.sc.res.0.2 <- bc.clusters.sc.res.0.2 | bc.clusters.seurat.sc.res.0.2

patch.clustering.res.0.2 <- bc.clusters.res.0.2 | bc.clusters.sc.res.0.2


# Plot spatial distribution
spatial.bc.clusters.res.0.2 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.res.0.2 <- spatial.bc.clusters.res.0.2 +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters")

spatial.bc.clusters.sc.res.0.2 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.sc.res.0.2 <- spatial.bc.clusters.sc.res.0.2 +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters (sc-SHC)")

patch.spatial.clustering.res.0.2 <- spatial.bc.clusters.res.0.2 / spatial.bc.clusters.sc.res.0.2



#SAVE PLOTS
ggsave(filename = "bc_clustersc3_kparam.png",
       plot = clustersc3.kparam,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "clustreeplot_beyondcell.png",
       plot = clustree.plot,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "bc_patch_clustering_res_0.1.png",
       plot = patch.clustering,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "bc_patch_spatial_clustering_res_0.1.png",
       plot = patch.spatial.clustering,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "bc_patch_clustering_res_0.2.png",
       plot = patch.clustering.res.0.2,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "bc_patch_spatial_clustering_res_0.2.png",
       plot = patch.spatial.clustering.res.0.2,
       path = "./results/plots/beyondcell_pure_breast/")

# SAVE DATA
saveRDS(bc.recomputed, file = paste0("./results/analysis/beyondcell_recomputed.rds"))

