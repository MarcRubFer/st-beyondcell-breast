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

#Read mechanisms of actions(MoAs)
collapse.moas <- read_tsv(file = "./data/collapsed_MoAs.tsv")

# Ranking
bc.merged@meta.data
bc.ranked.drugs <- bcRanks(bc.merged, idents = "bc_clusters_res.0.3")
bc4squares.drugs <- wrap_plots(bc4Squares(bc.ranked.drugs, idents = "bc_clusters_res.0.3"))
ggsave(filename = "bc4squares-drugs.png",
       plot = bc4squares.drugs,
       path = "./results/plots/beyondcell_drugs_rank/")

# Select TOP-Differential-Drugs
top.diff <- as.data.frame(bc.ranked.drugs@ranks) %>%
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


# Heatmap drugs matrix
#drugs.sigs <- which(grepl(pattern = "sig\\-", x = rownames(bc.merged@scaled)))
#drugs.names <- rownames(bc.merged@normalized[drugs.sigs,])

# Matrix of TOP-Differential Drugs
drugs.matrix <- bc.merged@normalized[top.diff.drugs,]
dim(drugs.matrix)
## Calculate maximum and minimum for matrix
drugs.max.matrix <- max(apply(drugs.matrix, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))

## Create cluster format for heatmap

##1 Extract the number of cluster
bc.names.clusters <- bc.merged@meta.data$bc_clusters_res.0.3

categorical.tumor.tme.nonorder <- sapply(bc.names.clusters, function(num) {
  if (num %in% c(0, 3, 4)) {
    return("Tumour")
  } else {
    return("TME")
  }
})
##2 Order the cluster vector
orden_clusters <- order(bc.names.clusters)
##3 Rearrange the drugs matrix in order by cluster
datos_ordenados_drugs <- drugs.matrix[, orden_clusters]
##4 Rearrange vector of cluster in order
clusters_ordenados <- bc.names.clusters[orden_clusters]

categorical.tumor.tme <- sapply(clusters_ordenados, function(num) {
  if (num %in% c(0, 3, 4)) {
    return("Tumour")
  } else {
    return("TME")
  }
})


## Create a dendogram only between clusters (0,1,2,...)
dend.drugs = cluster_between_groups(datos_ordenados_drugs, clusters_ordenados)


moas.levels <- levels(factor(merge.drugs.names$collapsed.MoAs))
col_vector <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#d627d9")
names(col_vector) <- moas.levels
col_vector
## Create heatmap
heatmap.drugs <- Heatmap(
  datos_ordenados_drugs,
  #drugs.matrix,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(Tumour_TME = categorical.tumor.tme,
                                     clusters = clusters_ordenados,
                                     col = list(clusters = c("0" = "cornflowerblue",
                                                             "1" = "goldenrod2",
                                                             "2" = "red3",
                                                             "3" = "seagreen4",
                                                             "4" = "darkorchid",
                                                             "5" = "darkorange1"),
                                                Tumour_TME = c("Tumour" = "red",
                                                               "TME" = "green"))),
  right_annotation = rowAnnotation(MoA = merge.drugs.names$collapsed.MoAs,
                                   Dual = merge.drugs.names$dual.MoAs,
                                   col = list(MoA = col_vector)),
  show_column_names = FALSE,
  column_split = categorical.tumor.tme,
  row_names_gp = gpar(fontsize = 6),
  row_labels = merge.drugs.names$preferred.and.sigs,
  row_split = 4,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs

# Matrix only tumor clusters(0,3 and 4)
tumour.spots <- bc.ranked.drugs@meta.data %>%
  filter(bc_clusters_res.0.3 == 0 | bc_clusters_res.0.3 == 3 | bc_clusters_res.0.3 == 4) %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.tumour <- drugs.matrix[,tumour.spots]

heatmap.drugs.tumour <- Heatmap(
  datos_ordenados_drugs,
  #drugs.matrix,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(Tumour_TME = categorical.tumor.tme,
                                     clusters = clusters_ordenados,
                                     col = list(clusters = c("0" = "cornflowerblue",
                                                             "1" = "goldenrod2",
                                                             "2" = "red3",
                                                             "3" = "seagreen4",
                                                             "4" = "darkorchid",
                                                             "5" = "darkorange1"),
                                                Tumour_TME = c("Tumour" = "red",
                                                               "TME" = "green"))),
  right_annotation = rowAnnotation(MoA = merge.drugs.names$collapsed.MoAs,
                                   Dual = merge.drugs.names$dual.MoAs,
                                   col = list(MoA = col_vector)),
  show_column_names = FALSE,
  column_split = categorical.tumor.tme,
  row_names_gp = gpar(fontsize = 6),
  row_labels = merge.drugs.names$preferred.and.sigs,
  row_split = 4,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs


png(filename = "./results/plots/beyondcell_drugs_rank/heatmap_drugs.png", width = 1920, height = 1080)
heatmap.drugs
dev.off()

# Functional Matrix
functional.sigs <- which(!grepl(pattern = "sig\\-", x = rownames(bc.merged@normalized)))
functional.names <- rownames(bc.merged@normalized[functional.sigs,])
functional.matrix <- bc.merged@normalized[functional.sigs,]


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
