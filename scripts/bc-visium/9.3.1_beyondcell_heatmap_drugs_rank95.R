rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")


out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

set.seed(1)

# Read SeuratObjects
seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")
bc.recomputed <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")

head(bc.recomputed@meta.data)

gs.breast <- GenerateGenesets(x = "./data/gmts/drug_signatures_classic_nodup.gmt")

# Drugs ranking
## Establish cutoff  of 5% (more restrictive)
bc.ranked <- bcRanks(bc.recomputed, idents = "TCs_res.0.3",  resm.cutoff = c(0.05,0.95))
bc.4squares <- bc4Squares(bc.ranked, idents = "TCs_res.0.3")
bc4squares.plots <- wrap_plots(bc.4squares)

head(bc.ranked@ranks)

# Select TOP-Differential-Drugs
top.diff <- as.data.frame(bc.ranked@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

df.top.diff <- as.data.frame(top.diff) 
df.top.diff <- df.top.diff %>%
  separate(col = top.diff, 
           into = c("drug","data.base","db.id"), 
           sep = "_",
           remove = F)
head(df.top.diff)
# Export top differential drugs
write.table(x = df.top.diff, 
            file = "./results/tables/top.differential.drugs.tsv", 
            sep = "\t",
            row.names = F)

names.drugs <- as.data.frame(FindDrugs(bc.ranked, x = top.diff)) %>%
  select(IDs, preferred.and.sigs)
rownames(names.drugs) <- names.drugs$IDs
head(names.drugs)

# HeatMap Drugs
# Matrix of TOP-Differential Drugs
dim(bc.ranked@normalized)
drugs.matrix <- bc.ranked@normalized[top.diff,]
dim(drugs.matrix)

col.order <- bc.ranked@meta.data %>%
  select(TCs_res.0.3, spot.collapse) %>%
  arrange(TCs_res.0.3, spot.collapse)

col.order.spots <- col.order  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.ordered <- drugs.matrix[,col.order.spots]
drugs.matrix <- drugs.matrix.ordered

## Calculate maximum and minimum for matrix
drugs.max.matrix <- max(apply(drugs.matrix, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))

# Collapsed moas for breast-SSc signature
## Import manual collapsed moas drugs
collapsed.moas <- read_tsv(file = "./data/tsv/collapsed.moas.top.differential.drugs - top.differential.drugs.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$top.diff

#collapsed.moas <- collapsed.moas[match(rownames(datos_ordenados_drugs),collapsed.moas$signature_complete),]
names.moas <- levels(factor(collapsed.moas$collapsed.MoAs))
length.moas <- length(names.moas)
#col.moas <- colorRampPalette(brewer.pal(12,name = "Paired"))(length.moas)
#names(col.moas) <- names.moas

col.moas <- c("#db4470",
              "#5cc151",
              "#9b58cf",
              "#9bb932",
              "#5c6fda",
              "#cca939",
              "#dd6fd1",
              "#508f36",
              "#c83f97",
              "#63c385",
              "#8850a1",
              "#da8b2f",
              "#5b70b5",
              "#c45926",
              "#43c4c4",
              "#d5433c",
              "#61a2da",
              "#926a2e",
              "#c190d8",
              "#317341",
              "#e286a5",
              "#4a9f7c",
              "#9c4c77",
              "#abb061",
              "#ae5050",
              "#697329",
              "#df936d")
names(col.moas) <- names.moas

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
## Create heatmap
heatmap.drugs <- Heatmap(
  t(scale(t(drugs.matrix))),
  #drugs.matrix,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order$TCs_res.0.3,
                                     "Cell type" = col.order$spot.collapse,
                                     col = list("TCs" = TC.colors,
                                                "Cell type" = colors.categories)),
  right_annotation = rowAnnotation(MoA = collapsed.moas$collapsed.MoAs,
                                  col = list(MoA = col.moas)),
  #left_annotation = rowAnnotation(labels = c(1:15)),
  #row_dend_reorder = TRUE,
  show_column_names = FALSE,
  #column_split = col.order$TCs_res.0.3,
  #column_order = col.order,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas$preferred.drug.names),
  #row_km = 16,
  #row_order = 1:16,
  row_split = 16,
  #show_row_dend = F,
  #row_title = NULL,
  #col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs
heatmap.drugs <- draw(heatmap.drugs, merge_legend = TRUE)
heatmap.drugs
png(filename = "./results/plots/Beyondcell_oct23_DrugRank/heatmap_drugs_allspots_breastsig_rank95.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs)
dev.off()

# Analysis of Tumour TCs (TCs 3 to 6)
TCs.tumours <- levels(bc.recomputed@meta.data$TCs_res.0.3)[3:6]
cells.TCs.tumours <- bc.recomputed@meta.data %>%
  filter(TCs_res.0.3 %in% TCs.tumours) %>%
  rownames_to_column(var = "spots") %>%
  pull(spots)

bc.TCs.tumour <- bcSubset(bc.recomputed, cells = cells.TCs.tumours)
bc.ranked.tumour <- bcRanks(bc.TCs.tumour, idents = "TCs_res.0.3",  resm.cutoff = c(0.05,0.95))
top.diff.tumour <- as.data.frame(bc.ranked.tumour@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

# HeatMap Drugs
# Matrix of TOP-Differential Drugs
dim(bc.ranked.tumour@normalized)
drugs.matrix.tumour <- bc.ranked.tumour@normalized[top.diff.tumour,]
dim(drugs.matrix.tumour)

col.order.tumour <- bc.ranked.tumour@meta.data %>%
  select(TCs_res.0.3, spot.collapse) %>%
  arrange(TCs_res.0.3, spot.collapse)

col.order.spots <- col.order.tumour  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.ordered <- drugs.matrix.tumour[,col.order.spots]
drugs.matrix.tumour <- drugs.matrix.ordered


## Calculate maximum and minimum for matrix
drugs.max.matrix <- max(apply(drugs.matrix.tumour, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix.tumour, 1, function(row) min(row)))

tmp <- t(scale(t(drugs.matrix.tumour)))
tmp2 <- scale(tmp)
heatmap.drugs.tumour <- Heatmap(
  #drugs.matrix.tumour,
  #tmp,
  tmp2,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order.tumour$TCs_res.0.3,
                                     "Cell type" = col.order.tumour$spot.collapse,
                                     col = list("TCs" = TC.colors[3:6],
                                                "Cell type" = colors.categories)),
  #right_annotation = rowAnnotation(MoA = collapsed.moas$collapsed.MoAs,
  #                                 col = list(MoA = col.moas)),
  #left_annotation = rowAnnotation(labels = c(1:15)),
  #row_dend_reorder = TRUE,
  show_column_names = FALSE,
  column_split = col.order.tumour$TCs_res.0.3,
  #column_order = col.order,
  row_names_gp = gpar(fontsize = 6),
  #row_labels = toupper(collapsed.moas$preferred.drug.names),
  #row_km = 16,
  #row_order = 1:16,
  cluster_rows = T,
  row_split = 6,
  #show_row_dend = F,
  #row_title = NULL,
  #col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs.tumour
heatmap.drugs.tumour <- draw(heatmap.drugs.tumour, merge_legend = TRUE)
heatmap.drugs.tumour
png(filename = "./results/plots/Beyondcell_oct23_DrugRank/heatmap_drugs_allspots_breastsig_rank95_TCtumour.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs.tumour)
dev.off()
