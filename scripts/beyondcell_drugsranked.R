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
# Read beyondcell object
bc <- readRDS("./results/analysis/beyondcellobject.rds")
head(bc@meta.data)

# Compute drugs ranking by selected clustering
bc.ranked <- bcRanks(bc, idents = "bc_clusters_res.0.3")


list.4squares <- bc4Squares(bc.ranked, idents = "bc_clusters_res.0.3")
list.4squares <- wrap_plots(list.4squares)



bc.clusters <- bcClusters(bc.ranked, UMAP = "Seurat", idents = "bc_clusters_res.0.3", pt.size = 0.5) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.3)")

diff.high <- FindDrugs(bc.ranked, x=c("DIHYDROROTENONE", "ERLOTINIB"))$IDs
diff.high.plots <- bcSignatures(bc.ranked, UMAP = "Seurat", signatures = list(values=diff.high), pt.size = 0.5)
diff.high.plots <- wrap_plots(diff.high.plots) + plot_annotation(title = "Differencial High Sensitive")

diff.low <- FindDrugs(bc.ranked, x = c("MITOXANTRONE", "STAUROSPORINE"))$IDs
diff.low.plots <- bcSignatures(bc.ranked, UMAP = "Seurat", signatures = list(values = diff.low), pt.size = 0.5)
diff.low.plots <- wrap_plots(diff.low.plots)


layout <- "
AAAAAA
AAAAAA
##BB##
##BB##
CCCCCC
CCCCCC
"

diff.high.plots / bc.clusters / diff.low.plots + plot_layout(design = layout)

a <- bc.ranked@expr.matrix
head(a)
((list.4squares[[1]] + bc.clusters) / (diff.low.plots[[1]] +diff.high.plots[[1]])) | 
  bcHistogram(bc.ranked, signatures = "sig-20938", idents = "bc_clusters_res.0.3") |
  bcHistogram(bc.ranked, signatures = "sig-20943", idents = "bc_clusters_res.0.3") 
(list.4squares[[2]] + bc.clusters) / (diff.low.plots[[1]] + plot_spacer())
(list.4squares[[3]] + bc.clusters) / (diff.low.plots[[1]] + diff.high.plots[[2]])
(list.4squares[[4]] + bc.clusters) / (plot_spacer() + diff.high.plots[[1]])
(list.4squares[[5]] + bc.clusters) / (diff.low.plots[[3]] + diff.high.plots[[1]])
(list.4squares[[6]] + bc.clusters) / (diff.low.plots[[1]] + plot_spacer())

bcHistogram(bc.ranked, signatures = "sig-20938", idents = "bc_clusters_res.0.3")
bcHistogram(bc.ranked, signatures = "sig-21167")

diff.drugs <- c(diff.high, diff.low,"sig-21167")
diff.drugs <- c("sig-21167")
diff.drugs.names <- FindDrugs(bc.ranked, x = diff.drugs)$preferred.drug.names
diff.matrix <- bc.ranked@normalized[diff.drugs,]
diff.dataframe <- data.frame(diff.matrix)

diff.drugs.union <- paste0(diff.drugs.names,"(",diff.drugs,")")

library(ComplexHeatmap)

bc.names.clusters <- bc.ranked@meta.data$bc_clusters_res.0.3

orden_clusters <- order(bc.names.clusters)
datos_ordenados <- diff.matrix[, orden_clusters]
clusters_ordenados <- bc.names.clusters[orden_clusters]

dend1 = cluster_between_groups(datos_ordenados, clusters_ordenados)
dend2 = cluster_within_group(datos_ordenados, clusters_ordenados)

cell.colors <- c("B-cells" = "cornflowerblue",
                 "CAFs" = "goldenrod2",
                 "Cancer Epithelial" = "red3",
                 "Endothelial" = "seagreen4",
                 "Myeloid" = "darkorchid",
                 "Normal Epithelial" = "hotpink",
                 "Plasmablasts" = "greenyellow",
                 "PVL" = "darkorange1",
                 "T-cells" = "steelblue4")

Heatmap(
  datos_ordenados,
  cluster_columns = dend2,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     col = list(clusters = c("0" = "cornflowerblue",
                                             "1" = "goldenrod2",
                                             "2" = "red3",
                                             "3" = "seagreen4",
                                             "4" = "darkorchid",
                                             "5" = "darkorange1"))),
  show_column_names = FALSE,
  row_labels = diff.drugs.union
)

bc.ranked@ranks$bc_clusters_res.0.3[bc.ranked@ranks$bc_clusters_res.0.3$rank.3 == 1,]
gs.ssc@genelist[["sig-21305"]]
bcSignatures(bc.ranked, UMAP = "Seurat", signatures = list(values="sig-21167"), pt.size = 1) + bc.clusters
bc.ranked@meta.data

bc.ranked <- bcRanks(bc.ranked, idents = "spot.composition.filter")
bc.ranked@ranks$spot.composition.filter
bc4Squares(bc.ranked, idents = "spot.composition.filter")
wrap_plots(bc4Squares(bc.ranked, idents = "spot.composition.filter"), ncol = 3)

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/beyondcell_drugs_rank/"))
ggsave(filename = "bc4squares.png",
       plot = list.4squares,
       path = "./results/plots/beyondcell_drugs_rank/")