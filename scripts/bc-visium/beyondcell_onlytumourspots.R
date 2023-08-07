rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")


out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

set.seed(1)

# Load beyondcell objects
bc.drugs <- readRDS("./results/analysis/beyondcellobject.rds")
bc.functional <- readRDS("./results/analysis/beyondcell_functional.rds")

# Read additional tiles
collapse.moas <- read_tsv(file = "./data/collapsed_MoAs.tsv")
gs.ssc <- GetCollection(SSc)
gs.functional <- GenerateGenesets(x = "./data/gmts/functional.gmt")

# Create merged bc.object
bc.merged <- bcMerge(bc1 = bc.drugs, bc2 = bc.functional)

# Filter "Pure_Tumour" spots
pure.tumour.spots <- bc.merged@meta.data %>%
  filter(spot.composition.filter == "Pure_Tumour") %>%
  rownames_to_column("spots") %>%
  pull(spots)

# Subset bc object with puretumourspots
bc.puretumour <- bcSubset(bc.merged, cells = pure.tumour.spots)

# Recalculate clusters
res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
bc.pure.clusters <- bcUMAP(bc.puretumour, pc = 20, k.neighbors = 40, res = res)
head(bc.pure.clusters@meta.data)

spatial.tumourspots <- bcClusters(bc.pure.clusters, UMAP = "beyondcell", idents = "spot.composition.filter", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
ggsave(filename = "spatial_onlytumour.png",
       plot = spatial.tumourspots,
       path = "./results/plots/beyondcell_onlyPuretumour/")
spatial.clusters.tumour <- bcClusters(bc.pure.clusters, UMAP = "beyondcell", idents = "bc_clusters_res.0.3", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
ggsave(filename = "spatial_clusters_tumour.png",
       plot = spatial.clusters.tumour,
       path = "./results/plots/beyondcell_onlyPuretumour/")
dimplot.clusters.tumour <- bcClusters(bc.pure.clusters, UMAP = "beyondcell", idents = "bc_clusters_res.0.3", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.3)")
ggsave(filename = "dimplot_clusters_tumour.png",
       plot = dimplot.clusters.tumour,
       path = "./results/plots/beyondcell_onlyPuretumour/")

# Drugs ranking
bc.pure.ranked <- bcRanks(bc.pure.clusters, idents = "bc_clusters_res.0.3")
bc.pure.4squares <- bc4Squares(bc.pure.ranked, idents = "bc_clusters_res.0.3")
bc4squares.plots <- wrap_plots(bc.pure.4squares)
ggsave(filename = "bc4squares.tumour.png",
       plot = bc4squares.plots,
       path = "./results/plots/beyondcell_onlyPuretumour/")

# Select TOP-Differential-Drugs
top.diff <- as.data.frame(bc.pure.ranked@ranks) %>%
  select(starts_with(match = "bc_clusters_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("bc_clusters_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

top.diff.drugs <- top.diff[which(grepl(pattern = "sig\\-", x = top.diff))]
names.drugs <- as.data.frame(FindDrugs(bc.merged, x = top.diff.drugs)) %>%
  select(IDs, preferred.and.sigs)
rownames(names.drugs) <- names.drugs$IDs
sigs.to.moas <- collapse.moas %>%
  filter(IDs %in% top.diff.drugs) %>%
  mutate(collapse = paste(collapsed.MoAs,"-",dual.MoAs)) 
rownames(sigs.to.moas) <- sigs.to.moas$IDs
merge.drugs.names <- merge(names.drugs,sigs.to.moas, by = "IDs")
rownames(merge.drugs.names) <- merge.drugs.names$IDs

# HeatMap Drugs
# Matrix of TOP-Differential Drugs
drugs.matrix <- bc.pure.ranked@normalized[top.diff.drugs,]
dim(drugs.matrix)
## Calculate maximum and minimum for matrix
drugs.max.matrix <- max(apply(drugs.matrix, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))

##1 Extract the number of cluster
bc.names.clusters <- bc.pure.ranked@meta.data$bc_clusters_res.0.3

##2 Order the cluster vector
orden_clusters <- order(bc.names.clusters)
##3 Rearrange the drugs matrix in order by cluster
datos_ordenados_drugs <- drugs.matrix[, orden_clusters]
##4 Rearrange vector of cluster in order
clusters_ordenados <- bc.names.clusters[orden_clusters]

moas.levels <- levels(factor(merge.drugs.names$collapsed.MoAs))
col_vector <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                "#aec7e8", "#ffbb78", "#98df8a", "#ff9896")
names(col_vector) <- moas.levels
col_vector
## Create heatmap
heatmap.drugs <- Heatmap(
  datos_ordenados_drugs,
  #drugs.matrix,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     col = list(clusters = c("0" = "cornflowerblue",
                                                             "1" = "goldenrod2",
                                                             "2" = "red3",
                                                             "3" = "seagreen4",
                                                             "4" = "darkorchid",
                                                             "5" = "darkorange1"))),
  right_annotation = rowAnnotation(MoA = merge.drugs.names$collapsed.MoAs,
                                   #Dual = merge.drugs.names$dual.MoAs,
                                   col = list(MoA = col_vector)),
  show_column_names = FALSE,
  #column_split = categorical.tumor.tme,
  row_names_gp = gpar(fontsize = 6),
  row_labels = merge.drugs.names$preferred.and.sigs,
  row_split = 4,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs

png(filename = "./results/plots/beyondcell_onlyPuretumour/heatmap_onlytumour.png", width = 3000, height = 2000, res = 100)
heatmap.drugs
dev.off()

View(bc.pure.ranked@expr.matrix)

a <- a["ERBB2", pure.tumour.spots]
a <- as.matrix(t(a))
View(a)

datos_ordenados_ERBB2 <- a[,orden_clusters]
erbb2.max.matrix <- max(apply(datos_ordenados_ERBB2, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))
heatmap.ERBB2 <- Heatmap(
  datos_ordenados_ERBB2,
  #drugs.matrix,
  name = "ERBB2",
  cluster_columns = FALSE,
  #top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     #col = list(clusters = c("0" = "cornflowerblue",
                                      #                       "1" = "goldenrod2",
                                       #                      "2" = "red3",
                                        #                     "3" = "seagreen4"))),
  #right_annotation = rowAnnotation(MoA = merge.drugs.names$collapsed.MoAs,
                                   #Dual = merge.drugs.names$dual.MoAs,
                                   #col = list(MoA = col_vector)),
  show_column_names = FALSE,
  #column_split = categorical.tumor.tme,
  row_names_gp = gpar(fontsize = 6),
  #row_labels = merge.drugs.names$preferred.and.sigs,
  #row_split = 4,
  #show_row_dend = F,
  row_title = NULL,
  #col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.ERBB2
