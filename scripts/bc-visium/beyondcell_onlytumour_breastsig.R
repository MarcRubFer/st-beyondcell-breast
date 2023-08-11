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

# Read single-cell experiment.
seuratobj.distances <- readRDS("./results/analysis/seuratobj.distances.rds")

# Reformat spot composition in two levels: Pure Tumour and Rest
spot.dual <- seuratobj.distances@meta.data %>%
  select(spot.composition.filter) %>%
  mutate(spot.dual = case_when(spot.composition.filter == "Pure_Tumour" ~ "Pure_Tumour",
                               TRUE ~ "Rest"),
         spot.composition.filter = NULL)

seuratobj.distances <- AddMetaData(seuratobj.distances, metadata = spot.dual)
seuratobj.distances@meta.data

seuratobj.pure <- subset(seuratobj.distances, subset = spot.dual == "Pure_Tumour")
seuratobj.pure@meta.data

# Read additional tiles
gs.breast <- GenerateGenesets(x = "./data/gmts/drug_signatures_classic_nodup.gmt")

# Compute beyondcell score with breast signatures
bc <- bcScore(seuratobj.pure, gs.breast, expr.thres = 0.1) 

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

# Run the bcUMAP function again, specifying the k.params
res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
bc.recomputed <- bcUMAP(bc.recomputed, pc = 20, k.neighbors = 40, res = res)
head(bc.recomputed@meta.data)




# Reanalysis of sc3stability for plotting
clustree.plot <- clustree(bc.recomputed@meta.data, 
                          prefix = "bc_clusters_res.",
                          node_colour = "sc3_stability") 

clustree.graph <- clustree(bc.recomputed@meta.data, 
                           prefix = "bc_clusters_res.",
                           prop_filter = 0.1, 
                           node_colour = "sc3_stability",
                           return = "graph")

max.stability <- clustree.graph %>%
  activate(nodes) %>%
  as.data.frame() %>%
  group_by(bc_clusters_res.) %>%
  summarise(median.stability = median(sc3_stability),
            n.clusters = length(unique(cluster)))

max.stability.plot <- ggplot(data=max.stability, aes(x=bc_clusters_res., y=median.stability, group = 1)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability")

clustree.plot <- clustree.plot + ggtitle(label = "Clustree with 40 k.param")

# Test significance of clusters (sc-SHC package)

## From Script
k.neighbors <- 40
chosen.res <- 0.2
# Round scaled data
scaled <- round(bc.recomputed@scaled * 100, 0)

# Compute Therapeutic Clusters (TCs)

original.clusters <- bc.recomputed@meta.data[colnames(scaled), ] %>%
  pull(all_of("bc_clusters_res.0.2")) %>%
  as.character()
new.clusters <- testClusters(scaled, cluster_ids = original.clusters,
                             batch = NULL, alpha = 0.05, 
                             num_features = 100, num_PCs = 20, 
                             parallel = FALSE)
# Add metadata
TCs <- as.data.frame(new.clusters[[1]])
colnames(TCs) <- "bc_clusters_new"
rownames(TCs) <- colnames(scaled)
TCs <- TCs %>%
  mutate(bc_clusters_new = str_remove(bc_clusters_new, pattern = "^new"),
         bc_clusters_new_renamed = case_when(bc_clusters_new == "2" ~ "3",
                                             bc_clusters_new == "3" ~ "2",
                                             TRUE ~ bc_clusters_new),
         bc_clusters_new = as.factor(bc_clusters_new),
         bc_clusters_new_renamed = as.factor(bc_clusters_new_renamed))
print(head(TCs))


bc.recomputed <- bcAddMetadata(bc.recomputed, TCs)
head(bc.recomputed@meta.data)



# Plot UMAPs Beyondcell and Seurat clustering
bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.2)")
bc.clusters.seurat <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (res 0.2)")
patch.bc.clusters <- bc.clusters | bc.clusters.seurat

bc.clusters.new <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (sc-SHC)")

bc.clusters.seurat.new <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_new", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (sc-SHC)")
patch.bc.clusters.new <- bc.clusters.new | bc.clusters.seurat.new

bc.clusters | bc.clusters.new

# Plot spatial distribution
spatial.bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters <- spatial.bc.clusters +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters")

spatial.bc.clusters.new <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_new_renamed", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.new <- spatial.bc.clusters.new +
  plot_annotation(title = "Spatial Distribution of Beyondcell clusters (sc-SHC)")

# Drugs ranking

bc.ranked <- bcRanks(bc.recomputed, idents = "bc_clusters_new_renamed")
bc.ranked.95 <- bcRanks(bc.recomputed, idents = "bc_clusters_new_renamed", resm.cutoff = c(0.05,0.95))
bc.4squares.95 <- bc4Squares(bc.ranked.95, idents = "bc_clusters_new_renamed")
bc4squares.plots <- wrap_plots(bc.4squares.95)

# Select TOP-Differential-Drugs
top.diff.10 <- as.data.frame(bc.ranked@ranks) %>%
  select(starts_with(match = "bc_clusters_new_renamed.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("bc_clusters_new_renamed.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

top.diff <- as.data.frame(bc.ranked.95@ranks) %>%
  select(starts_with(match = "bc_clusters_new_renamed.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("bc_clusters_new_renamed.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

# HeatMap Drugs
# Matrix of TOP-Differential Drugs
drugs.matrix.10 <- bc.ranked@normalized[top.diff.10,]
drugs.matrix <- bc.ranked.95@normalized[top.diff,]
dim(drugs.matrix)
## Calculate maximum and minimum for matrix
drugs.max.matrix.10 <- max(apply(drugs.matrix.10, 1, function(row) max(row)))
drugs.min.matrix.10 <- min(apply(drugs.matrix.10, 1, function(row) min(row)))

drugs.max.matrix <- max(apply(drugs.matrix, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))

##1 Extract the number of cluster
bc.clusters.10 <- bc.ranked@meta.data$bc_clusters_new_renamed
bc.clusters <- bc.ranked.95@meta.data$bc_clusters_new_renamed
##2 Order the cluster vector
orden_clusters.10 <- order(bc.clusters.10)
orden_clusters <- order(bc.clusters)
##3 Rearrange the drugs matrix in order by cluster
datos_ordenados_drugs.10 <- drugs.matrix.10[, orden_clusters.10]
datos_ordenados_drugs <- drugs.matrix[, orden_clusters]
##4 Rearrange vector of cluster in order
clusters_ordenados.10 <- bc.clusters.10[orden_clusters.10]
clusters_ordenados <- bc.clusters[orden_clusters]

# Create heatmap with annotations
heatmap.drugs.10 <- Heatmap(
  datos_ordenados_drugs.10,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     ERBB2 = ERBB2.matrix,
                                     EGFR = EGFR.matrix,
                                     col = list(clusters = c("1" = "tomato",
                                                             "2" = "olivedrab",
                                                             "3" = "turquoise2",
                                                             "4" = "blueviolet"),
                                                ERBB2 = col.anno.ERBB2)),
  #right_annotation = rowAnnotation(MoA = collapsed.moas$collapsed.MoAs,
  #                                 col = list(MoA = col.moas)),
  show_column_names = FALSE,
  column_split = clusters_ordenados,
  row_names_gp = gpar(fontsize = 6),
  #row_labels = collapsed.moas$preferred.drug.names,
  row_split = 5,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs.10


# Collapsed moas
collapsed.moas <- read_tsv(file = "./data/selected_breast_signatures - Hoja 3.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$signature_complete

collapsed.moas <- collapsed.moas[match(rownames(datos_ordenados_drugs),collapsed.moas$signature_complete),]

names.moas <- levels(factor(collapsed.moas$collapsed.MoAs))
length.moas <- length(names.moas)
col.moas <- brewer.pal(n=length.moas, name = "Paired")
names(col.moas) <- names.moas

# Expression of ERBB2
seuratobj.pure@assays$SCT@scale.data
ERBB2.matrix <- bc.ranked.95@expr.matrix["ERBB2",]
ERBB2.matrix <- ERBB2.matrix[colnames(datos_ordenados_drugs)]

expr.levels.ERBB2 <- seq(from=min(ERBB2.matrix), to=max(ERBB2.matrix))
col.ERBB2 <- brewer.pal(n = length(expr.levels.ERBB2), name = "Spectral")
col.anno.ERBB2 = colorRamp2(breaks = expr.levels.ERBB2, colors = rev(col.ERBB2))

# Expression EGFR

EGFR.matrix <- bc.ranked.95@expr.matrix["EGFR",]
EGFR.matrix <- EGFR.matrix[colnames(datos_ordenados_drugs)]

expr.levels.EGFR <- round(seq(from=min(EGFR.matrix), to=max(EGFR.matrix), length.out = 4), digits = 1)
col.EGFR <- brewer.pal(n = length(expr.levels.EGFR), name = "Spectral")
col.anno.EGFR = colorRamp2(breaks = expr.levels.EGFR, colors = rev(col.EGFR))

# Create heatmap with annotations
heatmap.drugs <- Heatmap(
  datos_ordenados_drugs,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     ERBB2 = ERBB2.matrix,
                                     EGFR = EGFR.matrix,
                                     col = list(clusters = c("1" = "tomato",
                                                             "2" = "olivedrab",
                                                             "3" = "turquoise2",
                                                             "4" = "blueviolet"),
                                                ERBB2 = col.anno.ERBB2,
                                                EGFR = col.anno.EGFR)),
  right_annotation = rowAnnotation(MoA = collapsed.moas$collapsed.MoAs,
                                   col = list(MoA = col.moas)),
  show_column_names = FALSE,
  column_split = clusters_ordenados,
  row_names_gp = gpar(fontsize = 6),
  row_labels = collapsed.moas$preferred.drug.names,
  row_split = 5,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs
heatmap.drugs <- draw(heatmap.drugs, merge_legend = TRUE)


# Store Heatmap as an objecto to work in patchwork. 
w = convertWidth(unit(1, "npc")*(9/10), "inch", valueOnly = TRUE)
h = convertHeight(unit(1, "npc")*(4/5), "inch", valueOnly = TRUE)
grob <- grid.grabExpr(draw(heatmap.drugs), width = w, height = h)

layout <- "
##CCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
##CCCCCCCCCCCCCCCC
"

(spatial.bc.clusters.new[[1]] / spatial.bc.clusters.new[[2]]) + 
  grob +
  plot_layout(design = layout)

