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
library("ComplexHeatmap")
library("circlize")


# stablish seed
set.seed(1)

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")
bc.recomputed <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")


reactome.gmt <- readGMT(x= "./data/gmts/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
reactome.bc.gmt <- GenerateGenesets(x= "./data/gmts/reactome_bc_modificado.gmt")
fgseaRes.TME <- read_tsv("./results/tables/TC1_TC2_fGSEA_results.tsv")

# Analysis of TME TCs (TCs 1 and 2)
TCs.TME <- levels(bc.recomputed@meta.data$TCs_res.0.3)[1:2]
cells.TCs.TME <- bc.recomputed@meta.data %>%
  filter(TCs_res.0.3 %in% TCs.TME) %>%
  rownames_to_column(var = "spots") %>%
  pull(spots)

seurat.TME <- subset(seuratobj.tcs, cells = cells.TCs.TME)
dim(seurat.TME)

bc.TME.reactome <- bcScore(sc = seurat.TME,
                           gs = reactome.bc.gmt,
                           expr.thres = 0.1)

# Number of NAs
n.NA <- data.frame(nNAs = colSums(is.na(bc.TME.reactome@normalized)),
                   row.names = colnames(bc.TME.reactome@normalized))
bc.TME.reactome <- bcAddMetadata(bc.TME.reactome, n.NA)
bc.TME.reactome@meta.data

# Filter out spots with a high percentage of NAs
bc.TME.reactome <- bcSubset(bc.TME.reactome, nan.cells = 0.95)

# Replace NAs by 0s
bc.TME.reactome@normalized[is.na(bc.TME.reactome@normalized)] <- 0
bc.TME.reactome <- bcRecompute(bc.TME.reactome, slot = "normalized")

bc.TME.reactome <- bcUMAP(bc.TME.reactome, 
                        pc = 10, 
                        k.neighbors = 20, 
                        res = 0.2)
# Extract name of top 50 pathways
top.50 <- fgseaRes %>%
  filter(padj < 0.05) %>%
  mutate(UP_DOWN = case_when(NES > 0 ~ "UP",
                             TRUE ~ "DOWN")) %>%
  group_by(UP_DOWN) %>%
  arrange(desc(NES)) %>%
  slice_max(NES, n = 50) %>%
  pull(pathway)

#Colors for cell types categories
colors.categories <- toupper(c(#"#4fafe3",
  "MYELOID" ="#be7adc", #violet
  "CAFS" = "#dec36f", #ocre
  "ENDOTHELIAL" = "#549f42", #green
  "LYMPHOID" = "#f1703a", #orange
  "OTHERS" ="#79696B", #grey
  "TUMOUR" = "#c4534e")) #dark.red
names(colors.categories)
colors.categories <- colors.categories[order(names(colors.categories))]
colors.categories

TC.colors <- c("TC-1" = "#00b2d7",
               "TC-2" = "#e5c22f",
               "TC-3" = "#903ca2",
               "TC-4" = "#3f8741",
               "TC-5" = "#ff7b00",
               "TC-6" = "#cb5c42")

# Extract normalized data for top50 pathway
pathway.matrix <- bc.TME.reactome@normalized[top.50,]
pathway.matrix.scaled <- bc.TME.reactome@scaled[top.50,]
dim(pathway.matrix)

col.order <- bc.TME.reactome@meta.data %>%
  select(TCs_res.0.3, spot.collapse) %>%
  arrange(TCs_res.0.3, spot.collapse)

col.order.spots <- col.order  %>%
  rownames_to_column("spots") %>%
  pull(spots)

pathway.matrix <- pathway.matrix[,col.order.spots]
pathway.matrix.scaled <- pathway.matrix.scaled[ ,col.order.spots]


## Create heatmap
heatmap.pathway <- Heatmap(
  #pathway.matrix,
  pathway.matrix.scaled,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order$TCs_res.0.3,
                                     "Cell type" = col.order$spot.collapse,
                                     #"Cell type" = anno_block(gp = gpar(fill = colors.categories),
                                     #                         labels = levels(as.factor(col.order$spot.collapse)),
                                     #                         labels_gp = gpar(fontsize = 5)),
                                     col = list("TCs" = TC.colors,
                                                "Cell type" = colors.categories)),
  #right_annotation = rowAnnotation(MoA = collapsed.moas$collapsed.MoAs,
  #                                 col = list(MoA = col.moas)),
  #left_annotation = rowAnnotation(labels = c(1:15)),
  #row_dend_reorder = TRUE,
  show_column_names = FALSE,
  column_split = col.order$TCs_res.0.3,
  column_title_gp = gpar(fontsize = 10),
  #column_order = col.order,
  row_names_gp = gpar(fontsize = 6),
  #row_labels = toupper(collapsed.moas$preferred.drug.names),
  #row_km = 6,
  #row_order = 1:16,
  row_split = 8,
  #show_row_dend = F,
  #row_title = NULL,
  #col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.pathway
