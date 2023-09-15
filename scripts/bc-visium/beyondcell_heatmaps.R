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

bc.recomputed <- readRDS("./results/analysis/beyondcell_recomputed.rds")


# Drugs ranking
# Establish cutoff  of 5% (more restrictive) and calculate for both resolutions
bc.ranked.95 <- bcRanks(bc.recomputed, idents = "bc_clusters_new_renamed_res_0.1", resm.cutoff = c(0.05,0.95))
names(bc.ranked.95@ranks)
bc.ranked.95 <- bcRanks(bc.ranked.95, idents = "bc_clusters_new_renamed_res_0.2", resm.cutoff = c(0.05,0.95))
names(bc.ranked.95@ranks)

# Plots bc4squares
bc4squares.95.res.0.1 <- bc4Squares(bc.ranked.95, idents = "bc_clusters_new_renamed_res_0.1") 
bc4squares.95.res.0.1.plots <- wrap_plots(bc4squares.95.res.0.1) +
  plot_layout(guides = "collect")

bc4squares.95.res.0.2 <- bc4Squares(bc.ranked.95, idents = "bc_clusters_new_renamed_res_0.2") 
bc4squares.95.res.0.2.plots <- wrap_plots(bc4squares.95.res.0.2) +
  plot_layout(guides = "collect")


# Select TOP-Differential-Drugs
top.diff.res.0.1 <- as.data.frame(bc.ranked.95@ranks) %>%
  select(starts_with(match = "bc_clusters_new_renamed_res_0.1.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("bc_clusters_new_renamed_res_0.1.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

top.diff.res.0.2 <- as.data.frame(bc.ranked.95@ranks) %>%
  select(starts_with(match = "bc_clusters_new_renamed_res_0.2.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("bc_clusters_new_renamed_res_0.2.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

# HeatMap Drugs
# Matrix of TOP-Differential Drugs

drugs.matrix.res.0.1 <- bc.ranked.95@normalized[top.diff.res.0.1,]
dim(drugs.matrix.res.0.1)

drugs.matrix.res.0.2 <- bc.ranked.95@normalized[top.diff.res.0.2,]
dim(drugs.matrix.res.0.2)

## Calculate maximum and minimum for matrix

drugs.max.matrix.res.0.1 <- max(apply(drugs.matrix.res.0.1, 1, function(row) max(row)))
drugs.min.matrix.res.0.1 <- min(apply(drugs.matrix.res.0.1, 1, function(row) min(row)))

drugs.max.matrix.res.0.2 <- max(apply(drugs.matrix.res.0.2, 1, function(row) max(row)))
drugs.min.matrix.res.0.2 <- min(apply(drugs.matrix.res.0.2, 1, function(row) min(row)))

##1 Extract the number of cluster
num.bc.clusters.res.0.1 <- bc.ranked.95@meta.data$bc_clusters_new_renamed_res_0.1
num.bc.clusters.res.0.2 <- bc.ranked.95@meta.data$bc_clusters_new_renamed_res_0.2

##2 Order the cluster vector
orden_clusters.res.0.1 <- order(num.bc.clusters.res.0.1)
orden_clusters.res.0.2 <- order(num.bc.clusters.res.0.2)

##3 Rearrange the drugs matrix in order by cluster
datos_ordenados_drugs_res.0.1 <- drugs.matrix.res.0.1[, orden_clusters.res.0.1]
datos_ordenados_drugs_res.0.2 <- drugs.matrix.res.0.2[, orden_clusters.res.0.2]

##4 Rearrange vector of cluster in order
clusters_ordenados.res.0.1 <- num.bc.clusters.res.0.1[orden_clusters.res.0.1]
clusters_ordenados.res.0.2 <- num.bc.clusters.res.0.2[orden_clusters.res.0.2]



# Collapsed moas
collapsed.moas <- read_tsv(file = "./data/selected_breast_signatures - Hoja 3.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$signature_complete

collapsed.moas.res.0.1 <- collapsed.moas[match(rownames(datos_ordenados_drugs_res.0.1),collapsed.moas$signature_complete),]
names.moas.res.0.1 <- levels(factor(collapsed.moas.res.0.1$collapsed.MoAs))
length.moas.res.0.1 <- length(names.moas.res.0.1)
col.moas.res.0.1 <- brewer.pal(n=length.moas.res.0.1, name = "Paired")
names(col.moas.res.0.1) <- names.moas.res.0.1

collapsed.moas.res.0.2 <- collapsed.moas[match(rownames(datos_ordenados_drugs_res.0.2),collapsed.moas$signature_complete),]
names.moas.res.0.2 <- levels(factor(collapsed.moas.res.0.2$collapsed.MoAs))
length.moas.res.0.2 <- length(names.moas.res.0.2)
col.moas.res.0.2 <- brewer.pal(n=length.moas.res.0.2, name = "Paired")
names(col.moas.res.0.2) <- names.moas.res.0.2

# Expression of ERBB2

ERBB2.matrix <- bc.ranked.95@expr.matrix["ERBB2",]

ERBB2.matrix.res.0.1 <- ERBB2.matrix[colnames(datos_ordenados_drugs_res.0.1)]
expr.levels.ERBB2.res.0.1 <- seq(from=min(ERBB2.matrix.res.0.1), to=max(ERBB2.matrix.res.0.1))
col.ERBB2.res.0.1 <- brewer.pal(n = length(expr.levels.ERBB2.res.0.1), name = "YlGnBu")
col.anno.ERBB2.res.0.1 = colorRamp2(breaks = expr.levels.ERBB2.res.0.1, colors = rev(col.ERBB2.res.0.1))

ERBB2.matrix.res.0.2 <- ERBB2.matrix[colnames(datos_ordenados_drugs_res.0.2)]
expr.levels.ERBB2.res.0.2 <- seq(from=min(ERBB2.matrix.res.0.2), to=max(ERBB2.matrix.res.0.2))
col.ERBB2.res.0.2 <- brewer.pal(n = length(expr.levels.ERBB2.res.0.2), name = "YlGnBu")
col.anno.ERBB2.res.0.2 = colorRamp2(breaks = expr.levels.ERBB2.res.0.2, colors = rev(col.ERBB2.res.0.2))

# Expression EGFR

EGFR.matrix <- bc.ranked.95@expr.matrix["EGFR",]

EGFR.matrix.res.0.1 <- EGFR.matrix[colnames(datos_ordenados_drugs_res.0.1)]
expr.levels.EGFR.res.0.1 <- round(seq(from=min(EGFR.matrix.res.0.1), to=max(EGFR.matrix.res.0.1), length.out = 4), digits = 1)
col.EGFR.res.0.1 <- brewer.pal(n = length(expr.levels.EGFR.res.0.1), name = "PuRd")
col.anno.EGFR.res.0.1 = colorRamp2(breaks = expr.levels.EGFR.res.0.1, colors = rev(col.EGFR.res.0.1))

EGFR.matrix.res.0.2 <- EGFR.matrix[colnames(datos_ordenados_drugs_res.0.2)]
expr.levels.EGFR.res.0.2 <- round(seq(from=min(EGFR.matrix.res.0.2), to=max(EGFR.matrix.res.0.2), length.out = 4), digits = 1)
col.EGFR.res.0.2 <- brewer.pal(n = length(expr.levels.EGFR.res.0.2), name = "PuRd")
col.anno.EGFR.res.0.2 = colorRamp2(breaks = expr.levels.EGFR.res.0.2, colors = rev(col.EGFR.res.0.2))

# Create heatmap with annotations
heatmap.drugs.res.0.1 <- Heatmap(
  datos_ordenados_drugs_res.0.1,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados.res.0.1,
                                     ERBB2 = ERBB2.matrix.res.0.1,
                                     EGFR = EGFR.matrix.res.0.1,
                                     col = list(clusters = c("1" = "tomato",
                                                             "2" = "olivedrab",
                                                             "3" = "turquoise2"),
                                                ERBB2 = col.anno.ERBB2.res.0.1,
                                                EGFR = col.anno.EGFR.res.0.1)),
  right_annotation = rowAnnotation(MoA = collapsed.moas.res.0.1$collapsed.MoAs,
                                   col = list(MoA = col.moas.res.0.1)),
  show_column_names = FALSE,
  column_split = clusters_ordenados.res.0.1,
  row_names_gp = gpar(fontsize = 6),
  row_labels = collapsed.moas.res.0.1$preferred.drug.names,
  row_split = 5,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix.res.0.1, 0, drugs.max.matrix.res.0.1), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix.res.0.1, 0, drugs.max.matrix.res.0.1))
)      
heatmap.drugs.res.0.1 <- draw(heatmap.drugs.res.0.1, merge_legend = TRUE)
heatmap.drugs.res.0.1

# Create heatmap with annotations
heatmap.drugs.res.0.2 <- Heatmap(
  datos_ordenados_drugs_res.0.2,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados.res.0.2,
                                     ERBB2 = ERBB2.matrix.res.0.2,
                                     EGFR = EGFR.matrix.res.0.2,
                                     col = list(clusters = c("1" = "tomato",
                                                             "2" = "olivedrab",
                                                             "3" = "turquoise2",
                                                             "4" = "blueviolet"),
                                                ERBB2 = col.anno.ERBB2.res.0.2,
                                                EGFR = col.anno.EGFR.res.0.2)),
  right_annotation = rowAnnotation(MoA = collapsed.moas.res.0.2$collapsed.MoAs,
                                   col = list(MoA = col.moas.res.0.2)),
  show_column_names = FALSE,
  column_split = clusters_ordenados.res.0.2,
  row_names_gp = gpar(fontsize = 6),
  row_labels = collapsed.moas.res.0.2$preferred.drug.names,
  row_split = 5,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix.res.0.2, 0, drugs.max.matrix.res.0.2), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix.res.0.2, 0, drugs.max.matrix.res.0.2))
)      
heatmap.drugs.res.0.2 <- draw(heatmap.drugs.res.0.2, merge_legend = TRUE)
heatmap.drugs.res.0.2

png(filename = "./results/plots/beyondcell_pure_breast/prueba.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
pdf(file = "./results/plots/beyondcell_pure_breast/prueba.pdf")
draw(heatmap.drugs.res.0.2)
dev.off()


# Store Heatmap as an object to work in patchwork. 
w = convertWidth(unit(1, "npc")*(9/10), "inch", valueOnly = TRUE)
h = convertHeight(unit(1, "npc")*(4/5), "inch", valueOnly = TRUE)
grob.res.0.1 <- grid.grabExpr(draw(heatmap.drugs.res.0.1), width = w, height = h)
grob.res.0.2 <- grid.grabExpr(draw(heatmap.drugs.res.0.2), width = w, height = h)


# Plot spatial distribution of clusters and its heatmap
spatial.bc.clusters.new.0.1 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.1", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.new.0.2 <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))

layout <- "
##CCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
##CCCCCCCCCCCCCCCC
"

spatial.heatmap.res.0.1 <- (spatial.bc.clusters.new.0.1[[1]] / spatial.bc.clusters.new.0.1[[2]]) + 
  grob.res.0.1 +
  plot_layout(design = layout)
spatial.heatmap.res.0.2 <- (spatial.bc.clusters.new.0.2[[1]] / spatial.bc.clusters.new.0.2[[2]]) + 
  grob.res.0.2 +
  plot_layout(design = layout)

# Save plots
ggsave(filename = "spatial_and_heatmap_beyondcell_res_01.png",
       plot = spatial.heatmap.res.0.1,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "spatial_and_heatmap_beyondcell_res_02.pdf",
       plot = spatial.heatmap.res.0.2,
       path = "./results/plots/beyondcell_pure_breast/")

ggsave(filename = "spatial_and_heatmap_beyondcell_res_01_nes.png",
       plot = spatial.heatmap.res.0.1,
       width = 47,
       height = 24,
       path = "./results/plots/beyondcell_pure_breast/")

# Save data
saveRDS(bc.ranked.95, file = paste0("./results/analysis/beyondcell_ranked95.rds"))


# For 10% cutoff bcRanks (default)

# cutoff 10%
#bc.ranked <- bcRanks(bc.recomputed, idents = "bc_clusters_new_renamed")

# Select TOP-Differential-Drugs
# For cutoff 10%
#top.diff.10 <- as.data.frame(bc.ranked@ranks) %>%
#  select(starts_with(match = "bc_clusters_new_renamed.group.")) %>%
#  rownames_to_column("signature") %>%
#  pivot_longer(cols = starts_with("bc_clusters_new_renamed.group."), names_to = "cluster", values_to = "group") %>%
#  filter(group != is.na(group),
#         grepl("Differential", group)) %>%
#  pull("signature") %>%
#  unique()

#drugs.matrix.10 <- bc.ranked@normalized[top.diff.10,]

#drugs.max.matrix.10 <- max(apply(drugs.matrix.10, 1, function(row) max(row)))
#drugs.min.matrix.10 <- min(apply(drugs.matrix.10, 1, function(row) min(row)))

#bc.clusters.10 <- bc.ranked@meta.data$bc_clusters_new_renamed
#orden_clusters.10 <- order(bc.clusters.10)
#datos_ordenados_drugs.10 <- drugs.matrix.10[, orden_clusters.10]
#clusters_ordenados.10 <- bc.clusters.10[orden_clusters.10]

## Create heatmap with annotations
#heatmap.drugs.10 <- Heatmap(
#  datos_ordenados_drugs.10,
#  name = "bcScore",
#  cluster_columns = FALSE,
#  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
#                                     #ERBB2 = ERBB2.matrix,
#                                     #EGFR = EGFR.matrix,
#                                     col = list(clusters = c("1" = "tomato",
#                                                             "2" = "olivedrab",
#                                                             "3" = "turquoise2",
#                                                             "4" = "blueviolet"))),
#  #right_annotation = rowAnnotation(MoA = collapsed.moas$collapsed.MoAs,
#  #                                 col = list(MoA = col.moas)),
#  show_column_names = FALSE,
#  column_split = clusters_ordenados,
#  row_names_gp = gpar(fontsize = 6),
#  #row_labels = collapsed.moas$preferred.drug.names,
#  row_split = 5,
#  #show_row_dend = F,
#  row_title = NULL,
#  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
#  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
#)      
#heatmap.drugs.10