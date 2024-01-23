# conda activate beyondcell2_EC

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
library("cowplot")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Set seed
set.seed(1)

# Load RDS
seuratobj.phase <- readRDS(file = "./results/analysis/seuratobj.phase.rds")

# RunPCA made in previous script
# Elbow plot
elbow <- ElbowPlot(seuratobj.phase, 
                   ndims = 50, 
                   reduction = "pca")
elbow
#ggsave(elbow, filename = paste0(out.dir, "/plots/elbow.pdf")) 

# Create SNN clustering graph using Seurat: (https://www.biostars.org/p/9572463/)
  # Find Neighboors: is a function that is used to find the nearest neighbors of your 
  #                   single cell data point within a dataset. It works by calculating the neighborhood 
  #                   overlap (Jaccard index) between every cell and its k. param nearest neighbors. 
  #                   It's often employed in various applications such as anomaly detection, and 
  #                   dimensionality reduction. 
  #                   The concept is that given a data point, you want to identify the closest data points 
  #                   to it based on some similarity metric, such as Euclidean distance or cosine similarity. 
  #                   This helps to identify similar points in the dataset, which can be useful for making 
  #                   predictions or understanding the distribution of the data.
  
  # FindClusters: is a function used for clustering data points into groups or 
  #               clusters based on their similarity. It uses a graph-based clustering approach 
  #               and a Louvain algorithm. Clustering is an unsupervised learning technique where 
  #               the algorithm groups similar cells together without any predefined labels. 
  #               The goal is to find patterns and structure in your data. 
  #               The number of clusters and the algorithm used can vary based on the problem and data characteristics

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




# Test significance of clusters (sc-SHC package)
# Round scaled data
scaled <- round(seuratobj.clusters@assays$SCT@data * 100, 0)

# Compute Therapeutic Clusters (TCs)
# We test significance at every resolution level previous selected. This parameter
# could be change to test some particular resolutions. (take some minutes)
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
scSHCs <- do.call(cbind, list.scSHC)  
head(scSHCs)
colnames(scSHCs) <- paste0("scSHC_res.",res)
rownames(scSHCs) <- colnames(scaled)
scSHCs <- scSHCs %>%
  mutate_all(.funs = ~str_replace_all(., "new", "scSHC-")) %>%
  mutate_all(.funs = ~as.factor(.))

# Add scSHCs data to metadata seurat object
seuratobj.clusters <- AddMetaData(seuratobj.clusters, metadata = scSHCs)
head(seuratobj.clusters@meta.data)

#Plot of statistic significance ggsankey SCT vs scSHC
sankey.scSHC <- seuratobj.clusters@meta.data %>%
  select(SCT_snn_res.0.2, scSHC_res.0.2) %>%
  make_long(SCT_snn_res.0.2, scSHC_res.0.2) %>%
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

# Create patch of stability and significance clustering
clustree.sankey <- (clustree.plot | sankey.scSHC) + plot_annotation(tag_levels = "A")

# Once a resolution has been chosen. We create a new variable with 
# final denomination of expression clusters (ECs)
head(seuratobj.clusters@meta.data)

ECs <- seuratobj.clusters@meta.data %>%
  select(scSHC_res.0.2) %>%
  mutate(scSHC_res.0.2 = str_replace(scSHC_res.0.2,pattern = "scSHC-",replacement = ""),
         ECs_res.0.2 = case_when(scSHC_res.0.2 == 5 ~ "EC-1",
                                 scSHC_res.0.2 == 4 ~ "EC-2",
                                 scSHC_res.0.2 == 8 ~ "EC-3",
                                 scSHC_res.0.2 == 9 ~ "EC-4",
                                 scSHC_res.0.2 == 7 ~ "EC-5",
                                 scSHC_res.0.2 == 6 ~ "EC-6",
                                 scSHC_res.0.2 == 2 ~ "EC-7",
                                 scSHC_res.0.2 == 3 ~ "EC-8",
                                 scSHC_res.0.2 == 1 ~ "EC-9",
                                 TRUE ~ scSHC_res.0.2),
         scSHC_res.0.2 = NULL)
head(ECs)

# Add new variable to metadata
seuratobj.clusters <- AddMetaData(seuratobj.clusters, metadata = ECs)
head(seuratobj.clusters@meta.data)

# UMAP_Dimplot and Spatial Dimplots of ECs.
# Create colors for expression clusters
ECs.colors <- c("EC-1" = "#ff944c",
                "EC-2" = "#76b74b",
                "EC-3" = "#8762c9",
                "EC-4" = "#f65b9a",
                "EC-5" = "#6d8cce",
                "EC-6" = "#b76512", 
                "EC-7" = "#4db497",
                "EC-8" = "#c26083",
                "EC-9" = "#6c823c")

dim.clusters <- DimPlot(seuratobj.clusters, group.by = "ECs_res.0.2", cols = ECs.colors) +
  ggtitle(label = NULL)




spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "ECs_res.0.2", combine = F) 
spatial.clusters <- lapply(spatial.clusters, function(i){
  i + scale_fill_manual(values = ECs.colors) + labs(fill = "ECs")
})
spatial.clusters <- wrap_plots(spatial.clusters)



# Boxplot by cluster, celltypes proportion
# boxplot.celltypes.clusters <- seuratobj.clusters@meta.data %>%
#   select(B.cells:T.cells,SCT_snn_res.0.3) %>%
#   rename(clusters = SCT_snn_res.0.3) %>%
#   mutate(clusters = as.factor(clusters)) %>%
#   pivot_longer(cols = -clusters,
#                names_to = "cell.type",
#                values_to = "cell.prop") %>%
#   ggplot(., aes(x = cell.type, y = cell.prop)) +
#   geom_boxplot(aes(fill = cell.type)) + 
#   facet_wrap(~clusters) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks = element_blank())

# Barplot cell types proportion by expr. clusters (ECs)
colors.categories <- toupper(c("MYELOID" ="#be7adc", #violet
                               "CAFS" = "#dec36f", #ocre
                               "ENDOTHELIAL" = "#549f42", #green
                               "LYMPHOID" = "#f1703a", #orange
                               "OTHERS" ="#79696B", #grey
                               "TUMOUR" = "#c4534e")) #dark.red

names(colors.categories)
colors.categories <- colors.categories[order(names(colors.categories))]
colors.categories

cell.types.by.EC <- seuratobj.clusters@meta.data %>%
  select(spot.collapse,ECs_res.0.2) %>%
  group_by(ECs_res.0.2, spot.collapse) %>%
  summarise(n = n())

total_counts <- cell.types.by.EC %>%
  group_by(ECs_res.0.2) %>%
  summarise(total = sum(n))

final.df <- cell.types.by.EC %>%
  left_join(total_counts, by = "ECs_res.0.2")

final.df <- final.df %>%
  mutate(relat.prop = n / total)

barplot.celltypes.clusters <- ggplot(final.df, aes(x=ECs_res.0.2, y=relat.prop, fill = spot.collapse)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = colors.categories) +
  ylab(label = "Relative proportion") +
  labs(fill = "Cell type categories") +
  theme_minimal() +
  theme(axis.line.y = element_line(colour = "grey"),
        axis.line.x = element_line(colour = "grey"),
        axis.title.y = element_text(vjust = 1.5),
        axis.title.x = element_blank())

barplot.celltypes.clusters



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
dim.dual <- DimPlot(seuratobj.clusters, group.by = "spot.dual", cols = dual.colors) +
  ggtitle(label = NULL)

spatial.dual <- SpatialDimPlot(seuratobj.clusters, group.by = "spot.dual", cols = dual.colors, combine = F)
spatial.dual <- lapply(spatial.dual, function(i) {
  i + labs(fill = "Tumour-TME")
})
spatial.dual <- wrap_plots(spatial.dual)

# Create a final figure for clustering
panel.a <- dim.clusters + 
  theme(legend.text = element_text(size = 8))
panel.b <- dim.dual +
  theme(legend.text = element_text(size = 8))
panel.c <- barplot.celltypes.clusters
panel.d <- spatial.dual + plot_layout(guides = "collect")
panel.e <- spatial.clusters + plot_layout(guides = "collect")

dim.panel <- plot_grid(panel.a,panel.b, align = "h", labels = c("A","B"))
left.panels <- plot_grid(dim.panel,panel.c, ncol = 1, labels = c("","C"), rel_heights = c(1,1.5))
right.panel <- plot_grid(panel.d,panel.e, labels = c("D","E"), ncol = 1)

figure <- left.panels | right.panel


# Save plots and ggplots
ggsave(filename = "clustree.png",
       plot = clustree.plot,
       path = "./results/plots/clustering/")
ggsave(filename = "clustree.svg",
       plot = clustree.plot,
       path = "./results/plots/clustering/")

ggsave(filename = "sankeydiagram_res_0.2.png",
       plot = sankey.scSHC,
       path = "./results/plots/clustering/")
ggsave(filename = "sankeydiagram_res_0.2.svg",
       plot = sankey.scSHC,
       path = "./results/plots/clustering/")

ggsave(filename = "clustree_sankey.png",
       plot = clustree.sankey,
       path = "./results/plots/clustering/")


ggsave(filename = "dimplot_clusters_res_0.2.png",
       plot = dim.clusters,
       path = "./results/plots/clustering/")
ggsave(filename = "dimplot_clusters_res_0.2.svg",
       plot = dim.clusters,
       path = "./results/plots/clustering/")

ggsave(filename = "spatialdimplot_clusters_res.0.2.png",
       plot = spatial.clusters,
       path = "./results/plots/clustering/")
ggsave(filename = "spatialdimplot_clusters_res.0.2.svg",
       plot = spatial.clusters,
       path = "./results/plots/clustering/")

ggsave(filename = "barplot_clusters_celltypes.png",
       plot = barplot.celltypes.clusters,
       path = "./results/plots/clustering/")
ggsave(filename = "barplot_clusters_celltypes.svg",
       plot = barplot.celltypes.clusters,
       path = "./results/plots/clustering/")

ggsave(filename = "dimplot_Tumor_TME.png",
       plot = dim.dual,
       path = "./results/plots/clustering/")
ggsave(filename = "dimplot_Tumor_TME.svg",
       plot = dim.dual,
       path = "./results/plots/clustering/")

ggsave(filename = "spatialdimplot_tumor_tme.png",
       plot = spatial.dual,
       path = "./results/plots/clustering/")
ggsave(filename = "spatialdimplot_tumor_tme.svg",
       plot = spatial.dual,
       path = "./results/plots/clustering/")

ggsave(filename = "figure_clustering.png",
       plot = figure,
       path = "./results/plots/clustering/")
ggsave(filename = "figure_clustering.svg",
       plot = figure,
       path = "./results/plots/clustering/")

all.plots <- list(clustree.plot,sankey.scSHC,clustree.sankey,dim.clusters,
                  spatial.clusters,barplot.celltypes.clusters,dim.dual,
                  spatial.dual,figure)
save(all.plots, file = paste0("./results/ggplots/clustering.RData"))

# Save data
saveRDS(seuratobj.clusters, file = "./results/analysis/seuratobj.clusters.rds")







