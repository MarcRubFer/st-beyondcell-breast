rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")


out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

set.seed(1)

# Read SeuratObjects
seuratobj.distances <- readRDS("./results/analysis/seuratobj.distances.rds")
seuratobj.distances.alt <- readRDS("./results/analysis/seuratobj.distances.alt.rds")

# Read Beyondcell objects (Uncomment the necessary objects)
# For bcScore compute with Ssc signature use these objects:
#bc.recomputed <- readRDS("./results/analysis/beyondcellobject.rds")
#bc.recomputed.alt <- readRDS("./results/analysis/beyondcellobject.alt.rds")

# For bcScore compute with breast-Ssc signature use these others:
bc.recomputed <- readRDS("./results/analysis/beyondcell_allspots_breastsignature.rds")
bc.recomputed.alt <- readRDS("./results/analysis/beyondcell_allspots_breastsignature.alt.rds")

head(bc.recomputed@meta.data)

# Drugs ranking
## Defatul cutoff at 10%
#bc.ranked <- bcRanks(bc.recomputed, idents = "bc_clusters_res.0.3")

## Establish cutoff  of 5% (more restrictive)
bc.ranked <- bcRanks(bc.recomputed, idents = "bc_clusters_res.0.3",  resm.cutoff = c(0.05,0.95))
bc.4squares <- bc4Squares(bc.ranked, idents = "bc_clusters_res.0.3")
bc4squares.plots <- wrap_plots(bc.4squares)

head(bc.ranked@ranks)

# Select TOP-Differential-Drugs
top.diff <- as.data.frame(bc.ranked@ranks) %>%
  select(starts_with(match = "bc_clusters_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("bc_clusters_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

names.drugs <- as.data.frame(FindDrugs(bc.ranked, x = top.diff)) %>%
  select(IDs, preferred.and.sigs)
rownames(names.drugs) <- names.drugs$IDs
head(names.drugs)

# HeatMap Drugs
# Matrix of TOP-Differential Drugs
drugs.matrix <- bc.ranked@normalized[top.diff,]
dim(drugs.matrix)
## Calculate maximum and minimum for matrix
drugs.max.matrix <- max(apply(drugs.matrix, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))

##1 Extract the number of cluster
heatmap.clusters <- bc.ranked@meta.data$bc_clusters_res.0.3
##2 Order the cluster vector
orden_clusters <- order(heatmap.clusters)
##3 Rearrange the drugs matrix in order by cluster
datos_ordenados_drugs <- drugs.matrix[, orden_clusters]
##4 Rearrange vector of cluster in order
clusters_ordenados <- heatmap.clusters[orden_clusters]

# Collapsed moas for breast-SSc signature
collapsed.moas <- read_tsv(file = "./data/selected_breast_signatures - Hoja 3.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$signature_complete

collapsed.moas <- collapsed.moas[match(rownames(datos_ordenados_drugs),collapsed.moas$signature_complete),]
names.moas <- levels(factor(collapsed.moas$collapsed.MoAs))
length.moas <- length(names.moas)
col.moas <- colorRampPalette(brewer.pal(12,name = "Paired"))(length.moas)
names(col.moas) <- names.moas

#Colors for cell types categories
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

## Create heatmap
heatmap.drugs <- Heatmap(
  datos_ordenados_drugs,
  #drugs.matrix2.order,
  #drugs.matrix,
  name = "bcScore",
  cluster_columns = FALSE,
  #top_annotation = HeatmapAnnotation(#df = col.order,
  #                                   TCs = col.order$bc_clusters_res.0.3,
  #                                   cell.type = col.order$spot.collapse),
  #top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
  #                                   cell.type = bc.ranked@meta.data$spot.collapse,
  #                                   col = list(cell.type = colors.categories)),
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados),
  #right_annotation = rowAnnotation(MoA = collapsed.moas$collapsed.MoAs,
  #                                col = list(MoA = col.moas)),
  show_column_names = FALSE,
  column_split = col.order$bc_clusters_res.0.3,
  #column_order = col.order,
  row_names_gp = gpar(fontsize = 6),
  #row_labels = collapsed.moas$preferred.drug.names,
  row_split = 5,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs
heatmap.drugs <- draw(heatmap.drugs, merge_legend = TRUE)
heatmap.drugs
png(filename = "./results/plots/Beyondcell_oct23_breastsig/heatmap_drugs_allspots_breastsig_rank95.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs)
dev.off()


col.order <- bc.ranked@meta.data %>%
  select(bc_clusters_res.0.3, spot.collapse) %>%
  arrange(bc_clusters_res.0.3, spot.collapse)
drugs.matrix2 <- drugs.matrix
drugs.matrix2.order <- drugs.matrix2[,col.order]
