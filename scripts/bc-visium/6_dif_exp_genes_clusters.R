rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("clustree")
library("tidygraph")

# Load seurat object

seuratobj.clusters <- readRDS(file = "./results/analysis/seuratobj.clusters.rds")

##cell.1 <- seuratobj.clusters@meta.data %>%
##  filter(SCT_snn_res.0.2 == 1)
##
##cell.2 <- seuratobj.clusters@meta.data %>%
##  filter((SCT_snn_res.0.2 != 1 & SCT_snn_res.0.2 != 0))
##
##FindMarkers(seuratobj.clusters, cells.1 = cell.1, cells.2 = cell.2)
##
##markers.c1 <- FindMarkers(seuratobj.clusters, ident.1 = 1, ident.2 = c(2,max(levels(seuratobj.clusters@meta.data$SCT_snn_res.0.2))))
##top_n(x = markers.c1,
##      n = 10,
##      wt = avg_log2FC)

#####

all.markers <- FindAllMarkers(seuratobj.clusters)
filter(.data = all.markers, gene == "LINC00052")
top20.genes <- all.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  mutate(rank = dense_rank(desc(avg_log2FC))) %>%
  filter(rank <= 20)

Idents(seuratobj.clusters) <- "SCT_snn_res.0.2"

heat.top20 <- DoHeatmap(seuratobj.clusters, features = top20.genes$gene, group.by = "SCT_snn_res.0.2") +
  theme(axis.text.y = element_text(size = 5))
heat.top20
ggsave(filename = "heatmap.png", plot = heat.top20, path = "./results/plots/", height = 10)

top20.cluster0 <- top20.genes %>%
  filter(cluster == 0)
FeaturePlot(seuratobj.clusters, features = top20.cluster0$gene)

top20.cluster1 <- top20.genes %>%
  filter(cluster == 1)
FeaturePlot(seuratobj.clusters, features = top20.cluster1$gene)

top20.cluster2 <- top20.genes %>%
  filter(cluster == 2) %>%
  arrange(rank)
FeaturePlot(seuratobj.clusters, features = top20.cluster2$gene)



VlnPlot(seuratobj.clusters, features = "percent.mt")
hist(seuratobj.clusters@meta.data$percent.mt)

# Differential between cluster 2 and rest of clusters

markers.cluster2 <- FindMarkers(seuratobj.clusters, ident.1 = 2)
markers.cluster2 <- markers.cluster2 %>%
  filter(p_val_adj < 0.05)
