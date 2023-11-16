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

reactome.gmt <- readGMT(x= "./data/gmts/c2.cp.reactome.v2023.2.Hs.symbols.gmt")

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

write.table(x = TME.dea.gsea,
            file = "./results/tables/TC1_TC2_findMarkers_results.tsv",
            sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

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

write.table(x = fgseaRes,
            file = "./results/tables/TC1_TC2_fGSEA_results.tsv",
            sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

fgseaRes <- as.data.frame(fgseaRes)

barplot.data <- fgseaRes %>%
  filter(padj < 0.05) %>%
  mutate(UP_DOWN = case_when(NES > 0 ~ "UP",
                             TRUE ~ "DOWN")) %>%
  group_by(UP_DOWN) %>%
  arrange(desc(NES)) %>%
  slice_max(NES, n = 50) 

barplot.data.export <- barplot.data %>%
  select(-leadingEdge)

write.table(x = barplot.data.export,
            file = "./results/tables/TC1_TC2_barplot_data.tsv",
            sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

barplot.data.annotated <- read_tsv(file = "./data/tsv/TC1_TC2_barplot_data_reannotated.tsv")

barplot.data <- barplot.data %>%
  mutate(REACTOME_METAPROGRAM = barplot.data.annotated$REACTOME_METAPROGRAM)

barplot <- barplot.data %>%
  ggplot(aes(y = reorder(as.factor(pathway), REACTOME_METAPROGRAM), x = NES)) +
    geom_bar(aes(fill = REACTOME_METAPROGRAM),stat = "identity") 
 
barplot <- barplot.data %>%
  ggplot(aes(y = REACTOME_METAPROGRAM)) +
  geom_bar(aes(fill=REACTOME_METAPROGRAM)) +
  scale_x_continuous(breaks = c(1:10)) 
barplot  
  
  
topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=20), pathway]
plotGseaTable(pathways = reactome.gmt[topPathwaysUp], 
              stats = log2fc.vector, 
              fgseaRes = fgseaRes, 
              gseaParam=0.5)

