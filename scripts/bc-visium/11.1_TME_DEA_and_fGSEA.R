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

# stablish seed
set.seed(1)

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Read SeuratObjects
#seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")

seurat.TME <- readRDS(file = "./results/analysis/seuratobj.TME-TCs.rds")
# Load gmts
reactome.gmt <- readGMT(x= "./data/gmts/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
hallmarks.extended.gmt <- readGMT(x = "./data/gmts/hallmarks.schuetz.brcaness.tis.gmt")


# Change to Spatial assay and establish TCs as idents
DefaultAssay(seurat.TME) <- "Spatial"
Idents(seurat.TME) <- "TCs_res.0.3"

# Check if results were obtained previously
path.to.results <- "./results/tables/TC1_TC2_findMarkers_results.tsv"

if (file.exists(path.to.results)) {
  # Load file from path
  TME.dea.gsea <- read_tsv(file = "./results/tables/TC1_TC2_findMarkers_results.tsv")
} else {
  # Obtain differential expressed markers/genes between TC-1 and TC-2
  TME.dea.gsea <- FindMarkers(seurat.TME, 
                            ident.1 = "TC-1",
                            ident.2 = "TC-2",
                            min.pct = 0, 
                            logfc.threshold = 0, 
                            test.use = "wilcox")
  # Move gene names (rownames) to a new column
  TME.dea.gsea <- TME.dea.gsea %>%
    rownames_to_column("gene")
  
  # Export results as tsv table
  write.table(x = TME.dea.gsea,
            file = "./results/tables/TC1_TC2_findMarkers_results.tsv",
            sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
  }

# Fast-GSEA: 
## For fgsea we need a named vector with the selected parameter (avg_log2FC)
## Extract gene names
names.vector.fgsea <- TME.dea.gsea %>%
  pull(gene)

## Extract avg_log2FC values
log2fc.vector <- TME.dea.gsea %>%
  select(avg_log2FC) %>%
  pull()

## Naming the vector
names(log2fc.vector) <- names.vector.fgsea

## Calculate GSEA results
fgseaRes <- fgsea(pathways = reactome.gmt, 
                  stats    = log2fc.vector,
                  minSize  = 15,
                  maxSize  = 500,
                  nPermSimple = 1000)

## TODO: create a table which include leadingedge list. At the moment, we save 
## fgseaRes without this column

fgseaRes.toexport <- fgseaRes %>%
  select(-leadingEdge)

write.table(x = fgseaRes.toexport,
            file = "./results/tables/TC1_TC2_fGSEA_results.tsv",
            sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


fgseaRes <- read_tsv("./results/tables/TC1_TC2_fGSEA_results.tsv")

fgseaRes <- as.data.frame(fgseaRes)

# In this case we have only two TCs, so we cannot use ggbubbleheatmap (multiple 
# comparations). For this results we decided to use a standard barplot


barplot.data <- fgseaRes %>%
  filter(padj < 0.05) %>%
  mutate(UP_DOWN = case_when(NES > 0 ~ "UP",
                             TRUE ~ "DOWN")) %>%
  group_by(UP_DOWN) %>%
  arrange(desc(NES)) %>%
  slice_max(NES, n = 50) 

barplot <- barplot.data %>%
  ggplot(aes(y = reorder(as.factor(pathway), NES), x = NES)) +
  geom_bar(aes(fill = "firebrick2"), stat = "identity") 


barplot.data.export <- barplot.data %>%
  select(-leadingEdge)

write.table(x = barplot.data.export,
            file = "./results/tables/TC1_TC2_barplot_data.tsv",
            sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)



barplot.data.annotated <- read_tsv(file = "./data/tsv/TC1_TC2_barplot_data_reannotated.tsv")

barplot.data2 <- barplot.data %>%
  mutate(REACTOME_METAPROGRAM = barplot.data.annotated$REACTOME_METAPROGRAM)

barplot2 <- barplot.data2 %>%
  ggplot(aes(y = reorder(as.factor(pathway), NES), x = NES)) +
    geom_bar(aes(fill = REACTOME_METAPROGRAM),stat = "identity") 
 
barplot3 <- barplot.data2 %>%
  ggplot(aes(y = REACTOME_METAPROGRAM)) +
  geom_bar(aes(fill=REACTOME_METAPROGRAM)) +
  scale_x_continuous(breaks = c(1:10)) 
barplot3

ggsave(filename = "TME_barplot_fGSEA_pathways.png",
       plot = barplot,
       path = "./results/plots/TC_TME_analysis/")
ggsave(filename = "TME_barplot_fGSEA_pathways_annotated.png",
       plot = barplot2,
       path = "./results/plots/TC_TME_analysis/")
ggsave(filename = "TME_barplot_fGSEA_metapathways.png",
       plot = barplot3,
       path = "./results/plots/TC_TME_analysis/")


###############################################################################
log2fc.vector.ordered <- sort(log2fc.vector, decreasing = T)
## Calculate GSEA results
fgseaRes.hallmarks <- fgsea(pathways = hallmarks.extended.gmt, 
                  stats    = log2fc.vector.ordered,
                  minSize  = 15,
                  maxSize  = 500,
                  nPermSimple = 1000)

## TODO: create a table which include leadingedge list. At the moment, we save 
## fgseaRes without this column

fgseaRes.hallmarks.toexport <- fgseaRes.hallmarks %>%
  select(-leadingEdge)

write.table(x = fgseaRes.hallmarks.toexport,
            file = "./results/tables/TC1_TC2_fGSEA_results_HALLMARKS.tsv",
            sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


fgseaRes.hallmarks <- read_tsv("./results/tables/TC1_TC2_fGSEA_results_HALLMARKS.tsv")

fgseaRes.hallmarks <- as.data.frame(fgseaRes.hallmarks)

barplot.hallmarks.data <- fgseaRes.hallmarks %>%
  filter(padj < 0.05) %>%
  mutate(UP_DOWN = case_when(NES > 0 ~ "UP",
                             TRUE ~ "DOWN")) %>%
  group_by(UP_DOWN) %>%
  arrange(desc(NES)) %>%
  slice_max(NES, n = 50) 

barplot.hallmarks <- barplot.hallmarks.data %>%
  ggplot(aes(y = reorder(as.factor(pathway), NES), x = NES)) +
  geom_bar(aes(fill = "firebrick2"), stat = "identity") 

ggsave(filename = "TME_barplot_HALLMARKS_fGSEA_pathways.png",
       plot = barplot.hallmarks,
       path = "./results/plots/TC_TME_analysis/")

myc.v1.founders <- readGMT(x = "./data/gmts/HALLMARK_MYC_TARGETS_V1_FOUNDERS.v2023.2.Hs.gmt")
fgseaRes.myc.founders <- fgsea(pathways = myc.v1.founders, 
                            stats    = log2fc.vector,
                            minSize  = 15,
                            maxSize  = 500,
                            nPermSimple = 1000)
data.myc.founders <- fgseaRes.myc.founders %>%
  filter(padj < 0.05) %>%
  mutate(UP_DOWN = case_when(NES > 0 ~ "UP",
                             TRUE ~ "DOWN")) %>%
  group_by(UP_DOWN) %>%
  arrange(desc(NES)) %>%
  slice_max(NES, n = 50) 
barplot.myc <- data.myc.founders %>%
  ggplot(aes(y = reorder(as.factor(pathway), NES), x = NES)) +
  geom_bar(aes(fill = "firebrick2"), stat = "identity")
