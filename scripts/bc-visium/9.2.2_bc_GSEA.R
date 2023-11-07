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


out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# stablish seed
set.seed(1)

# Read Reactome and functional GMT 
reactome.functional.gmt <- readGMT(x= "./data/gmts/reactome_functional.gmt")
reactome.brcaness.tis.gmt <- readGMT(x= "./data/gmts/reactome_brcaness_tis.gmt")
hallmarks.gmt <- readGMT(x="./data/gmts/hallmarks.all.v2023.1.Hs.symbols.gmt")
hallmarks.extended.gmt <- readGMT(x = "./data/gmts/hallmarks.schuetz.brcaness.tis.gmt")

seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")
bc.recomputed <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")

# Read table data for gsea

table.gsea <- read_tsv(file = "./data/dea_tables/dea.gsea.tsv")

# Generate data for BubbleHeatmap from Hallmarks GMT
num.tcs <- table.gsea %>%
  select(TC) %>%
  pull(TC) %>%
  unique()

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

results.gsea <- lapply(X = seq_along(num.tcs), FUN = function(tc) {
  vector.tc <- list.vector.gsea[[tc]]
  fgseaRes <- fgsea(pathways = hallmarks.extended.gmt, 
                    stats    = vector.tc,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 1000)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$TC <- tc
  return(fgseaRes)
})

gsea.table <- results.gsea %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC)
chosen.fun.sigs <- gsea.table %>%
  slice_max(order_by = NES, n = 50, by = COMPARISON) %>%
  pull(NAME) %>%
  unique()
table.bub.heatmap <- gsea.table %>%
  filter(!is.na(NES) & NAME %in% chosen.fun.sigs) %>%
  mutate(COMPARISON = factor(COMPARISON),
         COMPARISON = paste0("TC-",COMPARISON))

hallmarks.plot <- ggbubbleHeatmap(df = table.bub.heatmap,
                                  cluster.cols = F, 
                                  n.perm = 1000)
hallmarks.plot[[2]] <- hallmarks.plot[[2]] +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   face = "bold",
                                   size = 10))
hallmarks.plot


# Generate data for BubbleHeatmap from REACTOME GMT
results.gsea.reactome <- lapply(X = seq_along(num.tcs), FUN = function(tc) {
  vector.tc <- list.vector.gsea[[tc]]
  fgseaRes <- fgsea(pathways = reactome.functional.gmt, 
                    stats    = vector.tc,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 500)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$TC <- tc
  return(fgseaRes)
})

gsea.table.reactome <- results.gsea.reactome %>%
  bind_rows() %>%
  rename(NAME = pathway,
         FDR.q.val = padj, COMPARISON = TC)
chosen.fun.sigs.reactome <- gsea.table.reactome %>%
  slice_max(order_by = NES, n = 25, by = COMPARISON) %>%
  pull(NAME) %>%
  unique()
table.bub.heatmap.reactome <- gsea.table.reactome %>%
  filter(!is.na(NES) & NAME %in% chosen.fun.sigs.reactome) %>%
  mutate(COMPARISON = factor(COMPARISON))

reactome.plot <- ggbubbleHeatmap(df = table.bub.heatmap.reactome,
                                  cluster.cols = F, 
                                  n.perm = 500)
reactome.plot[[2]] <- reactome.plot[[2]] +
  theme(axis.text.y = element_text(size = 5,
                                   face = "bold"),
        axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 0,
                                   face = "bold",
                                   size = 10))
reactome.plot

# Patchwork spatial TC clustering and heatmaps
TC.colors <- c("TC-1" = "#00b2d7",
               "TC-2" = "#e3c26a",
               "TC-3" = "#903ca2",
               "TC-4" = "#3f8741",
               "TC-5" = "#ff7b00",
               "TC-6" = "#cb5c42")

spatial.TC.clusters <- SpatialDimPlot(seuratobj.tcs, group.by = "TCs_res.0.3", cols = TC.colors, ncol = 1)
spatial.TC.clusters <- spatial.TC.clusters + plot_layout(guides = "collect")
spatial.TC.clusters <- spatial.TC.clusters & theme(legend.position = "bottom",
                                                   legend.title = element_blank())

figure <- plot_grid(spatial.TC.clusters,NULL,hallmarks.plot, labels = c("A","","B"), rel_widths = c(0.15,0.05,1), label_y = 0.87, nrow = 1)
figure


ggsave(filename = "halmarks_bubbleheatmap.png",
       plot = hallmarks.plot,
       path = "./results/plots/Beyondcell_oct23_GSEA/")
ggsave(filename = "reactome_bubbleheatmap.png",
       plot = reactome.plot,
       path = "./results/plots/Beyondcell_oct23_GSEA/")
ggsave(filename = "figue_spatial_gseabubbleheatmap.png",
       plot = figure,
       path = "./results/plots/Beyondcell_oct23_GSEA/")




