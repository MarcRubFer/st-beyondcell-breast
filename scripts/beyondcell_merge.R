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
# Read beyondcell object
bc <- readRDS("./results/analysis/beyondcellobject.rds")
bc.functional <- readRDS("./results/analysis/beyondcell_functional.rds")

bc.merged <- bcMerge(bc1 = bc, bc2 = bc.functional)
gs.functional <- GenerateGenesets(x = "./data/gmts/functional.gmt")

functional.sigs <- which(!grepl(pattern = "sig\\-", x = rownames(bc.merged@scaled)))
functional.names <- rownames(bc.merged@normalized[functional.sigs,])
functional.matrix <- bc.merged@normalized[functional.sigs,]

functional.matrix.centered <- t(functional.matrix)
functional.matrix.centered <- scale(functional.matrix.centered, scale = FALSE)
functional.matrix.centered <- t(functional.matrix.centered)

row.neg <- apply(functional.matrix, 1, function(row) any(row < 0))
which(row.neg)
length(which(row.neg))

# Heatmap drugs matrix
drugs.sigs <- which(grepl(pattern = "sig\\-", x = rownames(bc.merged@scaled)))
drugs.names <- rownames(bc.merged@normalized[drugs.sigs,])
drugs.matrix <- bc.merged@normalized[drugs.sigs,]
drugs.max.matrix <- apply(drugs.matrix, 1, function(row) max(row))
drugs.max.matrix <- max(drugs.max.matrix)

drugs.min.matrix <- apply(drugs.matrix, 1, function(row) min(row))
drugs.min.matrix <- min(drugs.min.matrix)

bc.names.clusters <- bc.merged@meta.data$bc_clusters_res.0.3

orden_clusters <- order(bc.names.clusters)
datos_ordenados_drugs <- drugs.matrix[, orden_clusters]
clusters_ordenados <- bc.names.clusters[orden_clusters]

dend.drugs = cluster_between_groups(datos_ordenados_drugs, clusters_ordenados)

heatmap.drugs <- Heatmap(
  datos_ordenados_drugs,
  name = "bcScore",
  cluster_columns = dend.drugs,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     col = list(clusters = c("0" = "cornflowerblue",
                                                             "1" = "goldenrod2",
                                                             "2" = "red3",
                                                             "3" = "seagreen4",
                                                             "4" = "darkorchid",
                                                             "5" = "darkorange1"))),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 3),
  row_labels = drugs.names,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs

png(filename = "./results/plots/heatmap_drugs.png", width = 2560, height = 1440)
heatmap.drugs
dev.off()

row.neg <- apply(drugs.matrix, 1, function(row) any(row < 0))
which(row.neg)
length(which(row.neg))

# Maximum of functional matrix
max.matrix <- apply(functional.matrix, 1, function(row) max(row))
max.matrix <- max(max.matrix)

# Minimum of functional matrix
min.matrix <- apply(functional.matrix, 1, function(row) min(row))
min.matrix <- min(min.matrix)

# Heatmap of functional matrix, grouped by therapeutic clusters
bc.names.clusters <- bc.merged@meta.data$bc_clusters_res.0.3

orden_clusters <- order(bc.names.clusters)
datos_ordenados <- functional.matrix[, orden_clusters]
clusters_ordenados <- bc.names.clusters[orden_clusters]

dend1 = cluster_between_groups(datos_ordenados, clusters_ordenados)

heatmap.functional <- Heatmap(
  datos_ordenados,
  name = "bcScore",
  cluster_columns = dend1,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     col = list(clusters = c("0" = "cornflowerblue",
                                                             "1" = "goldenrod2",
                                                             "2" = "red3",
                                                             "3" = "seagreen4",
                                                             "4" = "darkorchid",
                                                             "5" = "darkorange1"))),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  row_labels = functional.names,
  col = colorRamp2(c(min.matrix, 0, max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(min.matrix, 0, max.matrix))
)      
heatmap.functional

plot.IFG <- bcSignatures(bc.merged, UMAP = "Seurat", signatures = list(values="HALLMARK_INTERFERON_GAMMA_RESPONSE"), pt.size = 0.5)
plot.EMT <- bcSignatures(bc.merged, UMAP = "Seurat", signatures = list(values="sig_EMT_HALLMARKS"), pt.size = 0.5)
plot.ductal <- bcSignatures(bc.merged, UMAP = "Seurat", signatures = list(values="SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE"), pt.size = 0.5)


schuetz.genes <- gs.functional@genelist$SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE$up
schuetz.genes.filt <- schuetz.genes[which(schuetz.genes %in% rownames(bc.merged@expr.matrix))]


exprss.matrix.schuetz <- bc.merged@expr.matrix[schuetz.genes.filt,]
exprss.matrix.schuetz <- scale(exprss.matrix.schuetz)
matrix.order <- exprss.matrix.schuetz[, orden_clusters]
dend2 = cluster_between_groups(matrix.order, clusters_ordenados)
heatmap.genes.schuetz <- Heatmap(
  matrix.order,
  cluster_columns = dend2,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     col = list(clusters = c("0" = "cornflowerblue",
                                                             "1" = "goldenrod2",
                                                             "2" = "red3",
                                                             "3" = "seagreen4",
                                                             "4" = "darkorchid",
                                                             "5" = "darkorange1"))),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  row_labels = schuetz.genes.filt
)      

# UMAPs

plot.emt.groger <- bcSignatures(bc.merged, UMAP = "Seurat", signatures = list(values="sig_EMT_GROGER_2012"), pt.size = 0.5)
bcSignatures(bc.merged, UMAP = "Seurat", signatures = list(values="sig-21330"), pt.size = 0.5)
# Ranking
bc.merged@meta.data
bc.ranked <- bcRanks(bc.merged, idents = "bc_clusters_res.0.3")
ranking <- as.data.frame(bc.ranked@ranks)
rownames(ranking)
bc4squares.drugs <- wrap_plots(bc4Squares(bc.ranked, idents = "bc_clusters_res.0.3"))
ggsave(filename = "bc4squares-drugs.png",
       plot = bc4squares.drugs,
       path = "./results/plots/beyondcell_drugs_rank/")

ranking.cluster3 <- ranking %>%
  select(bc_clusters_res.0.3.rank.3:bc_clusters_res.0.3.group.3)


bc.clusters <- as.data.frame(bc.merged@meta.data$bc_clusters_res.0.3)
rownames(bc.clusters) <- rownames(bc.merged@meta.data)
bc.clusters

bc.functional <- bcAddMetadata(bc.functional, metadata = bc.clusters)
bc.functional@meta.data
bc.functional@meta.data$bc.clusters <- bc.functional@meta.data$`bc.merged@meta.data$bc_clusters_res.0.3`
bc.functional@meta.data$`bc.merged@meta.data$bc_clusters_res.0.3` <- NULL

bc.functional.ranked <- bcRanks(bc.functional, idents = "bc.clusters")
wrap_plots(bc4Squares(bc.functional.ranked, idents = "bc.clusters", topnames = "sig_EMT_GROGER_2012"))
bc.functional.ranked@ranks[["sig_EMT_GROGER_2012"]]
bc.functional.ranked@normalized

a <- as.data.frame(bc.functional.ranked@ranks)
a["sig_EMT_GROGER_2012",]

df.groger <- a["sig_EMT_GROGER_2012",]
groger.residuals <- df.groger %>%
  select(starts_with("bc.clusters.residuals.mean")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("bc.clusters.residuals.mean"), names_to = "clusters", values_to = "residuals.mean") %>%
  mutate(clusters = gsub("bc.clusters.residuals.mean", replacement = "cluster", clusters))

ggplot(groger.residuals, aes(x=residuals.mean, y=clusters)) +
  geom_point() +
  geom_text(aes(label = paste(round(residuals.mean, digits = 1),clusters), vjust = -1.5))

bc.functional.ranked@switch.point["sig_EMT_GROGER_2012"]

# Save heatmaps
png(filename = "./results/plots/heatmap_functional.png", width = 2560, height = 1440)
heatmap.functional
dev.off()

png(filename = "./results/plots/heatmap_genes_schuetz.png", width = 2560, height = 1440)
heatmap.genes.schuetz
dev.off()
