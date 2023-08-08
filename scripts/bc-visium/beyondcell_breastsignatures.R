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

# Read single-cell experiment.
seuratobj.distances <- readRDS("./results/analysis/seuratobj.distances.rds")

spot.dual <- seuratobj.distances@meta.data %>%
  select(spot.composition.filter) %>%
  mutate(spot.dual = case_when(spot.composition.filter == "Pure_Tumour" ~ "Pure_Tumour",
                               TRUE ~ "Rest"),
         spot.composition.filter = NULL)

seuratobj.distances <- AddMetaData(seuratobj.distances, metadata = spot.dual)
seuratobj.distances@meta.data
# Read additional tiles
gs.breast <- GenerateGenesets(x = "./data/gmts/drug_signatures_classic_nodup.gmt")
#gs.ssc <- GetCollection(SSc)
#gs.functional <- GenerateGenesets(x = "./data/gmts/functional.gmt")


# Compute with breast signatures
bc <- bcScore(seuratobj.distances, gs.breast, expr.thres = 0.1) 

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

bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "spot.dual", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Beyondcell (res 0.3)")
bc.clusters.seurat <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "spot.dual", pt.size = 1) +
  ggtitle("BeyondCell Clusters - UMAP Seurat (res 0.3)")
patch.bc.clusters <- bc.clusters | bc.clusters.seurat

# Drugs ranking
bc.ranked <- bcRanks(bc.recomputed, idents = "spot.dual")
bc.4squares <- bc4Squares(bc.ranked, idents = "spot.dual")
bc4squares.plots <- wrap_plots(bc.4squares)


# Select TOP-Differential-Drugs
top.diff <- as.data.frame(bc.ranked@ranks) %>%
  select(starts_with(match = "spot.dual.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("spot.dual.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

names.drugs <- as.data.frame(FindDrugs(bc.ranked, x = top.diff)) %>%
  select(IDs, preferred.and.sigs)

  
  # HeatMap Drugs
  # Matrix of TOP-Differential Drugs
drugs.matrix <- bc.ranked@normalized[top.diff,]
dim(drugs.matrix)
## Calculate maximum and minimum for matrix
drugs.max.matrix <- max(apply(drugs.matrix, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))

##1 Extract the number of cluster
bc.clusters <- bc.ranked@meta.data$spot.dual
##2 Order the cluster vector
orden_clusters <- order(bc.clusters)
##3 Rearrange the drugs matrix in order by cluster
datos_ordenados_drugs <- drugs.matrix[, orden_clusters]
##4 Rearrange vector of cluster in order
clusters_ordenados <- bc.clusters[orden_clusters]

## Create heatmap
heatmap.drugs <- Heatmap(
  datos_ordenados_drugs,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     col = list(clusters = c("Pure_Tumour" = "cornflowerblue",
                                                             "Rest" = "red3"))),
  #right_annotation = rowAnnotation(MoA = merge.drugs.names$collapsed.MoAs,
  #                                 #Dual = merge.drugs.names$dual.MoAs,
  #                                 col = list(MoA = col_vector)),
  show_column_names = FALSE,
  column_split = clusters_ordenados,
  row_names_gp = gpar(fontsize = 6),
  #row_labels = merge.drugs.names$preferred.and.sigs,
  row_split = 7,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs

# Subset signatures which have differences between groups
moas.selected <- read_tsv(file = "./data/selected_breast_signatures.tsv")
rownames(moas.selected) <- moas.selected$signature
signatures.selected <- moas.selected %>%
  select(which(!grepl(pattern = "gene", x = names(moas.selected)))) %>%
  pull(signature)

datos_seleccionados <- datos_ordenados_drugs[signatures.selected,]
heatmap.subset <- Heatmap(
  datos_seleccionados,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(clusters = clusters_ordenados,
                                     col = list(clusters = c("Pure_Tumour" = "cornflowerblue",
                                                             "Rest" = "red3"))),
  right_annotation = rowAnnotation(MoA_PRISM = moas.selected$PRISM_MoA,
                                   MoA_CTRP = moas.selected$CTRP_MoA,
                                   MoA_GDSC = moas.selected$GDSC_MoA),
  show_column_names = FALSE,
  column_split = clusters_ordenados,
  row_names_gp = gpar(fontsize = 6),
  #row_labels = merge.drugs.names$preferred.and.sigs,
  row_split = 7,
  #show_row_dend = F,
  row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.subset

# Save plots
ggsave(filename = "bc4squares.png",
       plot = bc4squares.plots,
       path = "./results/plots/beyondcell_breastsignatures/")


