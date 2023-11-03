rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("RColorBrewer")
#devtools::install_github("igrabski/sc-SHC")
library("scSHC")
#install.packages("remotes")
#library("remotes")
#remotes::install_github("davidsjoberg/ggsankey")
library("ggsankey")
library("cowplot")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Set seed
set.seed(1)

# Read single-cell experiment.
seuratobj.aligned <- readRDS(file = "./results/analysis/seuratobj.aligned.rds")

# Set Assay.
DefaultAssay(seuratobj.aligned) <- "SCT"

# Generate geneset object with breast specific signature
gs.breast <- GenerateGenesets(x = "./data/gmts/drug_signatures_classic_nodup.gmt")

# Compute score for the SSc. This might take a few minutes depending on the size of your dataset.
bc <- bcScore(seuratobj.aligned, gs.breast, expr.thres = 0.1) 

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

# Update 0ct23: we fixed dims and k.param
# Step 1:
# Compute bcUMAP with no 'pc' parameter to check elbow diagram and select num.
  # of PCs. Rest of parameters default.
bc.recomputed <- bcUMAP(bc.recomputed, 
                        k.neighbors = 20, 
                        res = 0.2)

# Once analyzed elbow we establish pc = 10, k.parama=default(20) and test 
# different resolutions
res <- c(0.1, 0.2, 0.3, 0.4, 0.5)
bc.recomputed <- bcUMAP(bc.recomputed, 
                        pc = 10, 
                        k.neighbors = 20, 
                        res = res)
head(bc.recomputed@meta.data)

# Reanalysis of sc3stability for plotting
clustree.plot <- clustree(bc.recomputed@meta.data, 
                          prefix = "bc_clusters_res.",
                          node_colour = "sc3_stability") 

# Save to svg to add resolutions labels
ggsave(filename = "clustree.svg",
       plot = clustree.plot,
       path = "./results/plots/Beyondcell_oct23_breastsig//")

# Test significance of clusters (sc-SHC package)
# Round scaled data
scaled <- round(bc.recomputed@scaled * 100, 0)

# Compute Therapeutic Clusters (TCs)
# We test significance at every resolution level previous selected. This parameter
# could be change to test some particular resolutions.
res <- c(0.1, 0.2, 0.3, 0.4, 0.5)
list.scSHC <- lapply(X = res, FUN = function(x){
  original.clusters <- bc.recomputed@meta.data[colnames(scaled), ] %>%
    pull(all_of(paste0("bc_clusters_res.",x))) %>%
    as.character()
  
  new.clusters <- testClusters(scaled, 
                               cluster_ids = original.clusters,
                               batch = NULL, 
                               alpha = 0.05, 
                               num_features = 100, 
                               num_PCs = 20, 
                               parallel = FALSE)
  return(as.data.frame(new.clusters[[1]]))
  
})
scSHCs <- do.call(cbind, list.scSHC)  
head(scSHCs)
colnames(scSHCs) <- paste0("bc_scSHC_res.",res)
rownames(scSHCs) <- colnames(scaled)
scSHCs <- scSHCs %>%
  mutate_all(.funs = ~str_replace_all(., "new", "bc_scSHC-")) %>%
  mutate_all(.funs = ~as.factor(.))
head(scSHCs)

# Add scSHC results to metadata
bc.recomputed <- bcAddMetadata(bc.recomputed, metadata = scSHCs)
head(bc.recomputed@meta.data)
# Plot sankey diagram to study transitions between bc-clusters and scSHCs
list.plots.scSHC <- lapply(res, function(x){
  cluster <- paste0("bc_clusters_res.",x)
  scSHC <- paste0("bc_scSHC_res.",x)
  df <- bc.recomputed@meta.data %>%
    select(cluster, scSHC)
  df <- df %>%
    make_long(cluster, scSHC)
  pl <- ggplot(df, aes(x = x
                       , next_x = next_x
                       , node = node
                       , next_node = next_node
                       , fill = factor(node)
                       , label = node)
  )
  pl <- pl +geom_sankey(flow.alpha = 0.5
                        ,flow.fill = "grey"
                        ,node.color = "black"
                        ,node.fill = "#aab756"
                        ,show.legend = T)
  pl <- pl +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = 1)
  pl <- pl + theme_sankey() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  return(pl)
})
wrap_plots(list.plots.scSHC)

#Plot of statistic significance ggsankey SCT vs scSHC at res selected 0.3
sankey.scSHC <- bc.recomputed@meta.data %>%
  select(bc_clusters_res.0.3, bc_scSHC_res.0.3) %>%
  make_long(bc_clusters_res.0.3, bc_scSHC_res.0.3) %>%
  ggplot(., aes(x = x
                , next_x = next_x
                , node = node
                , next_node = next_node
                , fill = factor(node)
                , label = node)) + 
  geom_sankey(flow.alpha = 0.5
              #,flow.fill = "grey"
              ,node.color = "black"
              ,node.fill = "#aab756"
              ,show.legend = T) + 
  geom_sankey_label(size = 3, color = "black", fill= "white", hjust = 0.5) + 
  scale_x_discrete(position = "top") + 
  theme_sankey() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold"))
sankey.scSHC

clustree.sankey <- (clustree.plot | sankey.scSHC) + plot_annotation(tag_levels = "A")
# Save to svg to add resolutions labels
ggsave(filename = "clustree_sankey.svg",
       plot = clustree.sankey,
       path = "./results/plots/Beyondcell_oct23_breastsig//")

# Once a resolution has been chosen. We create a new variable with 
# final denomination of therapeutic clusters (TCs)
head(bc.recomputed@meta.data)

TCs <- bc.recomputed@meta.data %>%
  select(bc_scSHC_res.0.3) %>%
  mutate(bc_scSHC_res.0.3 = str_replace(bc_scSHC_res.0.3,pattern = "bc_scSHC-",replacement = ""),
         TCs_res.0.3 = case_when(bc_scSHC_res.0.3 == 3 ~ "TC-1",
                                 bc_scSHC_res.0.3 == 4 ~ "TC-2",
                                 bc_scSHC_res.0.3 == 2 ~ "TC-3",
                                 bc_scSHC_res.0.3 == 1 ~ "TC-4",
                                 bc_scSHC_res.0.3 == 5 ~ "TC-5",
                                 bc_scSHC_res.0.3 == 6 ~ "TC-6",
                                 TRUE ~ bc_scSHC_res.0.3),
         TCs_res.0.3 = as.factor(TCs_res.0.3),
         bc_scSHC_res.0.3 = NULL)
head(TCs)

bc.recomputed <- bcAddMetadata(bc.recomputed, metadata = TCs)
head(bc.recomputed@meta.data)
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

bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "TCs_res.0.3", pt.size = 1, cols = TC.colors) 
bc.clusters.seurat <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "TCs_res.0.3", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (res 0.1)")

TC.colors <- c("TC-1" = "#00b2d7",
               "TC-2" = "#e3c26a",
               "TC-3" = "#903ca2",
               "TC-4" = "#3f8741",
               "TC-5" = "#ff7b00",
               "TC-6" = "#cb5c42")
spatial.bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "TCs_res.0.3", pt.size = 1.5, spatial = TRUE) 
spatial.bc.clusters <- lapply(spatial.bc.clusters, function(x){
  x + scale_fill_manual(values = TC.colors) 
})
spatial.bc.clusters <- wrap_plots(spatial.bc.clusters) + plot_layout(guides = "collect")
spatial.bc.clusters

# Dual tumor-tme
dual.colors <- c(
  "TUMOUR" = "#c4534e",
  "TME" = "#098cf0")
bc.dual <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "spot.dual", pt.size = 1.5, cols = dual.colors)

spatial.dual <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "spot.dual", pt.size = 1.5, spatial = TRUE) 
spatial.dual <- lapply(spatial.dual, function(x){
  x + scale_fill_manual(values = dual.colors) 
})
spatial.dual <- wrap_plots(spatial.dual) + plot_layout(guides = "collect")
spatial.dual


# Barplot distribution of celltype in Therapeutic Clusters (TCs)
df.barplot <- bc.recomputed@meta.data %>%
  select(spot.collapse,TCs_res.0.3) %>%
  group_by(TCs_res.0.3, spot.collapse) %>%
  summarise(n = n(),
  )

total_counts <- df.barplot %>%
  group_by(TCs_res.0.3) %>%
  summarise(total = sum(n))

df.barplot <- df.barplot %>%
  left_join(total_counts, by = "TCs_res.0.3")

df.barplot <- df.barplot %>%
  mutate(relat.prop = n / total)


head(df.barplot)
colors.categories <- toupper(c(#"#4fafe3",
  "MYELOID" ="#be7adc", #violet
  "CAFS" = "#dec36f", #ocre
  "ENDOTHELIAL" = "#549f42", #green
  "LYMPHOID" = "#f1703a", #orange
  "OTHERS" ="#79696B", #grey
  "TUMOUR" = "#c4534e")) #dark.red
names(colors.categories)
colors.categories <- colors[order(names(colors.categories))]
colors.categories
barplot.celltypes.TCs <- ggplot(df.barplot, aes(x=TCs_res.0.3, y=relat.prop, fill = spot.collapse)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = colors.categories)
barplot.celltypes.TCs

ggsave(filename = "barplot_celltype_vs_TCs.png",
       plot = barplot.celltypes.TCs,
       path = "./results/plots/Beyondcell_oct23_breastsig//")

((((bc.nFeatureRNA / bc.nCountRNA) | bc.clusters | bc.clusters.seurat) / barplot.celltypes.TCs) | (spatial.bc.clusters[[1]] / spatial.bc.clusters[[2]])) + plot_layout(widths = c(2,1))




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
saveRDS(bc.recomputed, file = "./results/analysis/beyondcell_allspots_breastsignature.rds")
saveRDS(bc.recomputed.alt, file = "./results/analysis/beyondcell_allspots_breastsignature.alt.rds")
