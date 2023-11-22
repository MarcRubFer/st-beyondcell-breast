rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("cowplot")
library("viridis")
library("RColorBrewer")
library("fgsea")
library("ggseabubble")
library("ggtree")

# stablish seed
set.seed(1)

out.dir <- "./results/tables/"
dir.create(path = out.dir, recursive = TRUE)

# Subset seurat object for Tumour TCs cells
seurat.tumour <- readRDS(file = "./results/analysis/seuratobj.Tumour-TCs.rds")


# Load gmts
reactome.gmt <- readGMT(x= "./data/gmts/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
hallmarks.extended.gmt <- readGMT(x = "./data/gmts/hallmarks.schuetz.brcaness.tis.gmt")

DefaultAssay(seurat.tumour) <- "Spatial"
Idents(seurat.tumour) <- "TCs_res.0.3"

TCs.Tumour <- levels(seurat.tumour@meta.data$TCs_res.0.3)
tumour.dea.gsea <- lapply(TCs.Tumour, FUN = function(i) {
  markers <- FindMarkers(seurat.tumour, 
                         ident.1 = i, 
                         min.pct = 0, 
                         logfc.threshold = 0, 
                         test.use = "wilcox")
  markers <- markers %>%
    rownames_to_column("gene") %>%
    mutate(TC = i)
  return(markers)
}) %>%
  bind_rows()

write.table(tumour.dea.gsea, file = "./results/tables/tumour_TC_DEA.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

tumour.dea.gsea <- read_tsv(file = "./results/tables/tumour_TC_DEA.tsv")

# Fast-GSEA: 
# Process data for fgsea function
num.tcs <- tumour.dea.gsea %>%
  select(TC) %>%
  pull(TC) %>%
  unique()

list.vector.gsea <- lapply(num.tcs, function(tc) {
  names.genes <- tumour.dea.gsea %>%
    filter(TC == tc) %>%
    pull(gene)
  logfc <- tumour.dea.gsea %>%
    filter(TC == tc) %>%
    pull(avg_log2FC)
  names(logfc) <- names.genes
  logfc <- logfc[order(logfc, decreasing = TRUE)]
  return(logfc)
})

# fGSEA for Reactome GMT
results.gsea <- lapply(X = seq_along(num.tcs), FUN = function(tc) {
  vector.tc <- list.vector.gsea[[tc]]
  fgseaRes <- fgsea(pathways = reactome.gmt, 
                    stats    = vector.tc,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 100000)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$TC <- tc
  return(fgseaRes)
})

# Process fgsea results for Bubbleheatmap
gsea.table <- results.gsea %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC)

chosen.fun.sigs <- gsea.table %>%
  slice_max(order_by = NES, n = 30, by = COMPARISON) %>%
  pull(NAME) %>%
  unique()

table.bub.heatmap <- gsea.table %>%
  filter(!is.na(NES) & NAME %in% chosen.fun.sigs) %>%
  mutate(COMPARISON = COMPARISON + 2,
         COMPARISON = factor(COMPARISON),
         COMPARISON = paste0("TC-",COMPARISON))

# BubbleHeatmap plotting
reactome.plot <- ggbubbleHeatmap(df = table.bub.heatmap,
                                 cluster.cols = F, 
                                 n.perm = 10000)
reactome.plot[[2]] <- reactome.plot[[2]] +
  theme(axis.text.y = element_text(size = 5,
                                   face = "bold"),
        axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 0,
                                   face = "bold",
                                   size = 10))
reactome.plot

out.dir.plots <- "./results/plots/TC_Tumour_analysis"
dir.create(path = out.dir.plots, recursive = TRUE)
ggsave(filename = "reactome_bubbleheatmap.png",
       plot = reactome.plot,
       path = out.dir.plots)
##############################################################################
# fGSEA for Hallmarks GMT
results.gsea.hallmarks <- lapply(X = seq_along(num.tcs), FUN = function(tc) {
  vector.tc <- list.vector.gsea[[tc]]
  fgseaRes <- fgsea(pathways = hallmarks.extended.gmt, 
                    stats    = vector.tc,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 100000)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$TC <- tc
  return(fgseaRes)
})

# Process fgsea results for Bubbleheatmap
gsea.table.hallmarks <- results.gsea.hallmarks %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC)

chosen.fun.sigs.hallmarks <- gsea.table.hallmarks %>%
  slice_max(order_by = NES, n = 50, by = COMPARISON) %>%
  pull(NAME) %>%
  unique()

table.bub.heatmap.hallmarks <- gsea.table.hallmarks %>%
  filter(!is.na(NES) & NAME %in% chosen.fun.sigs.hallmarks) %>%
  mutate(COMPARISON = COMPARISON + 2,
         COMPARISON = factor(COMPARISON),
         COMPARISON = paste0("TC-",COMPARISON))

# BubbleHeatmap plotting
hallmarks.plot <- ggbubbleHeatmap(df = table.bub.heatmap.hallmarks,
                                 cluster.cols = F, 
                                 n.perm = 100000)
hallmarks.plot[[2]] <- hallmarks.plot[[2]] +
  theme(axis.text.y = element_text(size = 5,
                                   face = "bold"),
        axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 0,
                                   face = "bold",
                                   size = 10))
hallmarks.plot

out.dir.plots <- "./results/plots/TC_Tumour_analysis"
dir.create(path = out.dir.plots, recursive = TRUE)
ggsave(filename = "hallmarks_bubbleheatmap.png",
       plot = hallmarks.plot,
       path = out.dir.plots)

