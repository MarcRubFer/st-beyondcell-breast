rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("viridis")
library("RColorBrewer")

library("fgsea")
library("ggseabubble")
library("ggtree")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# stablish seed
set.seed(1)

# Read Reactome and functional GMT 
reactome.functional.gmt <- readGMT(x= "./data/gmts/reactome_functional.gmt")

seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")
bc.ranked.95 <- readRDS(file = "./results/analysis/beyondcell_ranked95.rds")

# Read table data for gsea
table.gsea.res.0.1 <- read_tsv(file = "./data/dea_tables/dea.gsea.res_0.1.tsv")
table.gsea.res.0.2 <- read_tsv(file = "./data/dea_tables/dea.gsea.res_0.2.tsv")

# Extract number of therapeutic clusters
num.tcs.res.0.1 <- table.gsea.res.0.1 %>%
  select(TC) %>%
  pull(TC) %>%
  unique()

num.tcs.res.0.2 <- table.gsea.res.0.2 %>%
  select(TC) %>%
  pull(TC) %>%
  unique()

# Get, as named vectors, the logFC from each TC (names = gene symbol)
list.vector.gsea.res.0.1 <- lapply(num.tcs.res.0.1, function(tc) {
  names.genes <- table.gsea.res.0.1 %>%
    filter(TC == tc) %>%
    pull(gene)
  logfc <- table.gsea.res.0.1 %>%
    filter(TC == tc) %>%
    pull(avg_log2FC)
  names(logfc) <- names.genes
  logfc <- logfc[order(logfc, decreasing = TRUE)]
  return(logfc)
})

list.vector.gsea.res.0.2 <- lapply(num.tcs.res.0.2, function(tc) {
  names.genes <- table.gsea.res.0.2 %>%
    filter(TC == tc) %>%
    pull(gene)
  logfc <- table.gsea.res.0.2 %>%
    filter(TC == tc) %>%
    pull(avg_log2FC)
  names(logfc) <- names.genes
  logfc <- logfc[order(logfc, decreasing = TRUE)]
  return(logfc)
})

# Compute fastGSEA (fgsea) for each TC
#results.gsea <- lapply(X = list.vector.gsea, FUN = function(vector.tc) {
#  fgseaRes <- fgsea(pathways = reactome.functional.gmt, 
#                    stats    = vector.tc,
#                    minSize  = 15,
#                    maxSize  = 500)
#  fgseaRes <- as.data.frame(fgseaRes)
#  return(fgseaRes)
#})

results.gsea.res.0.1 <- lapply(X = num.tcs.res.0.1, FUN = function(tc) {
  vector.tc <- list.vector.gsea.res.0.1[[tc]]
  fgseaRes <- fgsea(pathways = reactome.functional.gmt, 
                    stats    = vector.tc,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 5000)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$TC <- tc
  return(fgseaRes)
})

results.gsea.res.0.2 <- lapply(X = num.tcs.res.0.2, FUN = function(tc) {
  vector.tc <- list.vector.gsea.res.0.2[[tc]]
  fgseaRes <- fgsea(pathways = reactome.functional.gmt, 
                    stats    = vector.tc,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 100000)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$TC <- tc
  return(fgseaRes)
})

# Format table
gsea.table.res.0.1 <- results.gsea.res.0.1 %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC)

gsea.table.res.0.2 <- results.gsea.res.0.2 %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC)

# Get the top 20 signatures
chosen.fun.sigs.res.0.1 <- gsea.table.res.0.1 %>%
  slice_max(order_by = abs(NES), n = 50, by = COMPARISON) %>%
  #top_n(n = 20, wt = abs(NES)) %>%
  pull(NAME) %>%
  unique()

table.bub.heatmap.res.0.1 <- gsea.table.res.0.1 %>%
  filter(!is.na(NES) & NAME %in% chosen.fun.sigs.res.0.1) %>%
  mutate(COMPARISON = factor(COMPARISON))


chosen.fun.sigs.res.0.2 <- gsea.table.res.0.2 %>%
  slice_max(order_by = abs(NES), n = 50, by = COMPARISON) %>%
  #top_n(n = 20, wt = abs(NES)) %>%
  pull(NAME) %>%
  unique()

table.bub.heatmap.res.0.2 <- gsea.table.res.0.2 %>%
  filter(!is.na(NES) & NAME %in% chosen.fun.sigs.res.0.2) %>%
  mutate(COMPARISON = factor(COMPARISON))

# Plot bubbleheatmap
bub.heatmap.res.0.1 <- ggbubbleHeatmap(table.bub.heatmap.res.0.1, cluster.cols = TRUE, n.perm = 5000) 
bub.heatmap.res.0.1[[2]] <- bub.heatmap.res.0.1[[2]] +
  theme(axis.text.y = element_text(size = 5))
bub.heatmap.res.0.1 

bub.heatmap.res.0.2 <- ggbubbleHeatmap(table.bub.heatmap.res.0.2, cluster.cols = TRUE, n.perm = 100000) +
  theme(axis.text.y = element_text(size = 3))
bub.heatmap.res.0.2[[2]] <- bub.heatmap.res.0.2[[2]] +
  theme(axis.text.y = element_text(size = 5))
bub.heatmap.res.0.2

spatial.bc.clusters.new.res.0.1 <- bcClusters(bc.ranked.95, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.1", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))
spatial.bc.clusters.new.res.0.2 <- bcClusters(bc.ranked.95, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.2", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))

layout <- "
####CCCCCCCCCCCCCCCC
AAAACCCCCCCCCCCCCCCC
AAAACCCCCCCCCCCCCCCC
BBBBCCCCCCCCCCCCCCCC
BBBBCCCCCCCCCCCCCCCC
####CCCCCCCCCCCCCCCC
"

spatial.heatmap.res.0.1 <- (spatial.bc.clusters.new.res.0.1[[1]] / spatial.bc.clusters.new.res.0.1[[2]]) + 
  bub.heatmap.res.0.1 +
  plot_layout(design = layout)

spatial.heatmap.res.0.2 <- (spatial.bc.clusters.new.res.0.2[[1]] / spatial.bc.clusters.new.res.0.2[[2]]) + 
  bub.heatmap.res.0.2 +
  plot_layout(design = layout)

# Save bubble heatmaps
ggsave(filename = "GSEA_spatial_bubbleheatmap_res_0.1.png",
       plot = spatial.heatmap.res.0.1,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "GSEA_spatial_bubbleheatmap_res_0.2.png",
       plot = spatial.heatmap.res.0.2,
       path = "./results/plots/beyondcell_pure_breast/")

############################################################################################################3

gsea.table.res.0.1.sub <- gsea.table.res.0.1 %>%
  filter(!grepl(pattern = "^REACTOME", x = NAME))
fig.sub <- ggbubbleHeatmap(gsea.table.res.0.1.sub, cluster.cols = TRUE, n.perm = 100000) 


# GSEA only breast sigs
breast.functional.gmt <- readGMT(x= "./data/gmts/breast_functional.gmt")
results.gsea.breast <- lapply(X = num.tcs, FUN = function(tc) {
  vector.tc <- list.vector.gsea[[tc]]
  fgseaRes <- fgsea(pathways = breast.functional.gmt, 
                    stats    = vector.tc,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 100000)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$TC <- tc
  return(fgseaRes)
})
# Format table
gsea.table.breast <- results.gsea.breast %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC) %>%
  filter(!is.na(NES))
fig.breast <- ggbubbleHeatmap(gsea.table.breast, cluster.cols = TRUE, n.perm = 100000)



spatial.bc.clusters.new.res.0.1 <- bcClusters(bc.ranked.95, UMAP = "beyondcell", idents = "bc_clusters_new_renamed_res_0.1", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))

layout <- "
##CCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
##CCCCCCCCCCCCCCCC
"

spatial.heatmap <- (spatial.bc.clusters.new.res.0.1[[1]] / spatial.bc.clusters.new.res.0.1[[2]]) + 
  bub.heatmap +
  plot_layout(design = layout)

spatial.sub <- (spatial.bc.clusters.new[[1]] / spatial.bc.clusters.new[[2]]) | fig.sub 
spatial.breast <- (spatial.bc.clusters.new[[1]] / spatial.bc.clusters.new[[2]]) | fig.breast 

# Save gsea results as tables

# Tables global GSEA
for (index in 1:length(num.tcs)) {
  leading_edge_lists <- results.gsea2[[index]]$leadingEdge
  leading_edge_strings <- lapply(leading_edge_lists, paste, collapse = "|")
  leading.edge.df <- as.data.frame(unlist(leading_edge_strings))
  df_tmp <- results.gsea2[[index]] %>%
    select(-leadingEdge) %>%
    mutate(leadingEdge = leading.edge.df$`unlist(leading_edge_strings)`)
  file_name <- paste0("./data/dea_tables/fgsea.results.TC.", index, ".tsv")
  write.table(df_tmp, file = file_name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# Tables breast signatures GSEA
for (index in 1:length(num.tcs)) {
  leading_edge_lists <- results.gsea.breast[[index]]$leadingEdge
  leading_edge_strings <- lapply(leading_edge_lists, paste, collapse = "|")
  leading.edge.df <- as.data.frame(unlist(leading_edge_strings))
  df_tmp <- results.gsea.breast[[index]] %>%
    select(-leadingEdge) %>%
    mutate(leadingEdge = leading.edge.df$`unlist(leading_edge_strings)`)
  file_name <- paste0("./data/dea_tables/fgsea.breast.results.TC.", index, ".tsv")
  write.table(df_tmp, file = file_name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}


