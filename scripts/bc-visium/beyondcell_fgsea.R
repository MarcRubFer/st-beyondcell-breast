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
table.gsea <- read_tsv(file = "./data/dea_tables/dea.gsea.tsv")

# Extract number of therapeutic clusters
num.tcs <- table.gsea %>%
  select(TC) %>%
  pull(TC) %>%
  unique()

# Get, as named vectors, the logFC from each TC (names = gene symbol)
list.vector.gsea <- lapply(num.tcs, function(tc) {
  names.genes <- table.gsea %>%
    filter(TC == tc) %>%
    pull(gene)
  logfc <- table.gsea %>%
    filter(TC == tc) %>%
    pull(avg_log2FC)
  names(logfc) <- names.genes
  logfc <- logfc[order(logfc, decreasing = TRUE)]
  return(logfc)
})

positive <- sapply(num.tcs, function(x) {length(which(list.vector.gsea[[x]] > 0))})
negative <- sapply(num.tcs, function(x) {length(which(list.vector.gsea[[x]] < 0))})

distrib.fc <- cbind(num.tcs,positive,negative)

# Compute fastGSEA (fgsea) for each TC
#results.gsea <- lapply(X = list.vector.gsea, FUN = function(vector.tc) {
#  fgseaRes <- fgsea(pathways = reactome.functional.gmt, 
#                    stats    = vector.tc,
#                    minSize  = 15,
#                    maxSize  = 500)
#  fgseaRes <- as.data.frame(fgseaRes)
#  return(fgseaRes)
#})

results.gsea2 <- lapply(X = num.tcs, FUN = function(tc) {
  vector.tc <- list.vector.gsea[[tc]]
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
gsea.table <- results.gsea2 %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC)

# Get the top 20 signatures
chosen.fun.sigs <- gsea.table %>%
  slice_max(order_by = abs(NES), n = 50, by = COMPARISON) %>%
  #top_n(n = 20, wt = abs(NES)) %>%
  pull(NAME) %>%
  unique()

gsea.table <- gsea.table %>%
  filter(!is.na(NES) & NAME %in% chosen.fun.sigs) %>%
  mutate(COMPARISON = factor(COMPARISON))

# Plot bubbleheatmap
fig2B <- ggbubbleHeatmap(gsea.table, cluster.cols = TRUE, n.perm = 100000) +
  theme(axis.text.y = element_text(size = 6))
fig2B[[2]] <- fig2B[[2]] +
  theme(axis.text.y = element_text(size = 6))
fig2B 

gsea.table.sub <- gsea.table %>%
  filter(!grepl(pattern = "^REACTOME", x = NAME))
fig.sub <- ggbubbleHeatmap(gsea.table.sub, cluster.cols = TRUE, n.perm = 100000) 


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



spatial.bc.clusters.new <- bcClusters(bc.ranked.95, UMAP = "beyondcell", idents = "bc_clusters_new_renamed", pt.size = 1.5, spatial = TRUE, mfrow = c(1,2))

layout <- "
##CCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
AACCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
BBCCCCCCCCCCCCCCCC
##CCCCCCCCCCCCCCCC
"

spatial.heatmap <- (spatial.bc.clusters.new[[1]] / spatial.bc.clusters.new[[2]]) + 
  fig2B +
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

# Save bubble heatmaps
ggsave(filename = "global.bubbleheatmap.png",
       plot = spatial.heatmap,
       path = "./results/plots/beyondcell_pure_breast/")
ggsave(filename = "breast.sigs.bubbleheatmap.png",
       plot = spatial.breast,
       path = "./results/plots/beyondcell_pure_breast/")
