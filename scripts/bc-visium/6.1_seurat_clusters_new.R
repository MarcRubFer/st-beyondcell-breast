rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("clustree")
library("tidygraph")
#install.packages("remotes")
#library("remotes")
#remotes::install_github("davidsjoberg/ggsankey")
library("ggsankey")
#devtools::install_github("igrabski/sc-SHC")
library("scSHC")

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

#Find Neighboors and FindClusters at different resolutions
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

max.stability.plot | clustree.plot
# Test significance of clusters (sc-SHC package)
# Round scaled data
scaled <- round(seuratobj.clusters@assays$SCT@data * 100, 0)

# Compute Therapeutic Clusters (TCs)
# We test significance at every resolution level previous selected. This parameter
# could be change to test some particular resolutions.
list.scSHC <- lapply(X = res, FUN = function(x){
  original.clusters <- seuratobj.clusters@meta.data %>%
    pull(all_of(paste0("SCT_snn_res.",x))) %>%
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
ECs <- do.call(cbind, list.scSHC)  
head(ECs)
colnames(ECs) <- paste0("EC_scSHC_res.",res)
rownames(ECs) <- colnames(scaled)
ECs <- ECs %>%
  mutate_all(.funs = ~str_replace_all(., "new", "EC-scSHC-")) %>%
  mutate_all(.funs = ~as.factor(.))

seuratobj.clusters <- AddMetaData(seuratobj.clusters, metadata = ECs)
head(seuratobj.clusters@meta.data)

# Collapse all non-tumour cell types in TME category. Add metadata
dual.spot <- seuratobj.clusters@meta.data %>%
  select(spot.collapse) %>%
  mutate(spot.dual = case_when(spot.collapse == "TUMOUR" ~ spot.collapse,
                               TRUE ~ "TME"))
head(dual.spot)
seuratobj.clusters <- AddMetaData(seuratobj.clusters, metadata = dual.spot)
head(seuratobj.clusters@meta.data)
#Plot of statistic significance ggsankey SCT vs scSHC
list.plots.scSHC <- lapply(res, function(x){
  cluster <- paste0("SCT_snn_res.",x)
  scSHC <- paste0("EC_scSHC_res.",x)
  df <- seuratobj.clusters@meta.data %>%
    select(cluster, scSHC,spot.dual)
  df <- df %>%
    make_long(cluster, scSHC,spot.dual)
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

# Boxplot by cluster, celltypes proportion
df.celltypes.clusters <- seuratobj.clusters@meta.data %>%
  select(B.cells:T.cells,SCT_snn_res.0.3) %>%
  rename(clusters = SCT_snn_res.0.3) %>%
  mutate(clusters = as.factor(clusters)) %>%
  pivot_longer(cols = -clusters,
               names_to = "cell.type",
               values_to = "cell.prop")

boxplot.celltypes.clusters <- ggplot(df.celltypes.clusters, aes(x = cell.type, y = cell.prop)) +
  geom_boxplot(aes(fill = cell.type)) + 
  facet_wrap(~clusters) +
  #ggtitle(label = "Cell types proportions within each cluster", subtitle = "Cell cycle regression: standard") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) 
df <- seuratobj.clusters@meta.data %>%
  select(spot.collapse,SCT_snn_res.0.2) %>%
  group_by(SCT_snn_res.0.2, spot.collapse) %>%
  summarise(n = n(),
  )

total_counts <- df %>%
  group_by(SCT_snn_res.0.2) %>%
  summarise(total = sum(n))

df <- df %>%
  left_join(total_counts, by = "SCT_snn_res.0.2")

df <- df %>%
  mutate(relat.prop = n / total)


head(df)

colors.categories <- toupper(c(#"#4fafe3",
  "MYELOID" ="#be7adc", #violet
  "CAFS" = "#dec36f", #ocre
  "ENDOTHELIAL" = "#549f42", #green
  "LYMPHOID" = "#f1703a", #orange
  "OTHERS" ="#79696B", #grey
  "TUMOUR" = "#c4534e")) #dark.red

names(colors.categories)
colors.categories <- colors.categories[order(names(colors.categories))]
colors.categories
barplot.celltypes.clusters <- ggplot(df, aes(x=SCT_snn_res.0.2, y=relat.prop, fill = spot.collapse)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = colors.categories)

barplot.celltypes.clusters
spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1")
dim.clusters <- DimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1")

#df.celltypes.clusters.alt <- seuratobj.clusters.alt@meta.data %>%
#  select(B.cells:T.cells,SCT_snn_res.0.2) %>%
#  rename(clusters = SCT_snn_res.0.2) %>%
#  mutate(clusters = as.factor(clusters)) %>%
#  pivot_longer(cols = -clusters,
#               names_to = "cell.type",
#               values_to = "cell.prop")
#
#boxplot.celltypes.clusters.alt <- ggplot(df.celltypes.clusters.alt, aes(x = cell.type, y = cell.prop)) +
#  geom_boxplot(aes(fill = cell.type)) + 
#  facet_wrap(~clusters) +
#  ggtitle(label = "Cell types proportions within each cluster", subtitle = "Cell cycle regression: alternative") +
#  theme(axis.text.x = element_blank(),
#        axis.ticks = element_blank(), 
#        plot.subtitle = element_text(colour = "red")) 
#
#spatial.clusters.alt <- SpatialDimPlot(seuratobj.clusters.alt, group.by = "SCT_snn_res.0.3")
#dim.clusters.alt <- DimPlot(seuratobj.clusters.alt, group.by = "SCT_snn_res.0.2")


dim.clusters | (boxplot.celltypes.clusters / spatial.clusters)
DimPlot(seuratobj.clusters, group.by = "spot.collapse")

# Collapse all non-tumour cell types in TME category. Add metadata
dual.spot <- seuratobj.clusters@meta.data %>%
  select(spot.collapse) %>%
  mutate(spot.dual = case_when(spot.collapse == "TUMOUR" ~ spot.collapse,
                               TRUE ~ "TME"))
head(dual.spot)
seuratobj.clusters <- AddMetaData(seuratobj.clusters, metadata = dual.spot)

#DimPlot TME-Tumour
dual.colors <- c(
  "TUMOUR" = "#c4534e",
  "TME" = "#098cf0")
dim.dual <- DimPlot(seuratobj.clusters, group.by = "spot.dual", cols = dual.colors)
spatial.dual <- SpatialDimPlot(seuratobj.clusters, group.by = "spot.dual", cols = dual.colors)

dim.clusters <- DimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1")
dim.clusters.2 <- DimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2")

spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1",combine = F)
spatial.clusters <- lapply(seq_along(spatial.clusters), FUN = function(i){
  spatial.clusters[[i]] +
    scale_fill_manual(values = ECs.colors)
})
wrap_plots(spatial.clusters)
spatial.clusters.2 <- SpatialDimPlot(seuratobj.clusters, group.by = "EC_scSHC_res.0.2", combine = F) 
spatial.clusters.2 <- lapply(seq_along(spatial.clusters.2), FUN = function(x){
  spatial.clusters.2[[x]] +
    scale_fill_manual(values = ECs.colors)
})
wrap_plots(spatial.clusters)/wrap_plots(spatial.clusters.2)
ECs.colors <- c("#ce57ad",
                "#76b74b",
                "#8762c9",
                "#c59442",
                "#6d8cce",
                "#cb5842",
                "#4db497",
                "#c26083",
                "#6c823c")

panel.dimplots <- (dim.dual | dim.clusters.2) & ggtitle(label = NULL)
left.panels <- ((panel.dimplots / plot_spacer() / barplot.celltypes.clusters) + plot_layout(heights = c(1,0.1,1.5)))
right.panels <- ((spatial.dual[[1]] / spatial.dual[[2]]) | (spatial.clusters.2[[1]] / spatial.clusters.2[[2]])) 
(left.panels | right.panels) + plot_annotation(tag_levels = "A")

(clustree.plot / list.plots.scSHC[[2]]) | left.panels | right.panels

((clustree.plot / list.plots.scSHC[[2]]) + plot_layout(heights = c(2,1))) |  ((spatial.dual)/wrap_plots(spatial.clusters.2)/barplot.celltypes.clusters)


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






