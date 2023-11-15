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

# Read SeuratObjects
seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")

reactome.gmt <- readGMT(x= "./data/gmts/c2.cp.reactome.v2023.1.Hs.symbols.gmt")

# stablish seed
set.seed(1)

DefaultAssay(seuratobj.tcs) <- "Spatial"
Idents(seuratobj.tcs) <- "TCs_res.0.3"
TME.dea.gsea <- FindMarkers(seuratobj.tcs, 
                            ident.1 = "TC-1",
                            ident.2 = "TC-2",
                            min.pct = 0, 
                            logfc.threshold = 0, 
                            test.use = "wilcox")

names.vector.fgsea <- TME.dea.gsea %>%
  rownames_to_column("gene") %>%
  pull(gene)

log2fc.vector <- TME.dea.gsea %>%
  select(avg_log2FC) %>%
  pull()

names(log2fc.vector) <- names.vector.fgsea
fgseaRes <- fgsea(pathways = reactome.gmt, 
                  stats    = log2fc.vector,
                  minSize  = 15,
                  maxSize  = 500,
                  nPermSimple = 1000)

topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=20), pathway]
plotGseaTable(pathways = reactome.gmt[topPathwaysUp], 
              stats = log2fc.vector, 
              fgseaRes = fgseaRes, 
              gseaParam=0.5)

obj@reductions$pca.rev@feature.loadings

E <- seuratobj.tcs@reductions$pca@feature.loadings
seuratobj.tcs@assays$SCT@data

set.seed(1)
gesecaRes <- geseca(reactome.gmt, E, minSize = 15, maxSize = 500, center = FALSE)
head(gesecaRes, 10)
topPathways <- gesecaRes[, pathway] |> head(4)
titles <- sub("HALLMARK_", "", topPathways)
ps <- plotCoregulationProfileSpatial(reactome.gmt[topPathways], 
                                     seuratobj.tcs,
                                     title=titles)
cowplot::plot_grid(plotlist=ps, ncol=2)

plotGesecaTable(gesecaRes = gesecaRes, pathways = reactome.gmt[topPathways], E = E)



x <- GetAssay(seuratobj.tcs, assay=DefaultAssay(object))
E <- x@scale.data

res <- object
pathway <- reactome.gmt[[2]]
pathway <- intersect(unique(pathway), rownames(E))

score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
score <- scale(score, center=TRUE, scale=T)
