rm(list = ls())

library("beyondcell")
library("Seurat")
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

saveRDS(object = bc.ranked, file = "./results/analysis/bc_ranked_all.rds")

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

# Scaled matrix
dim(bc.ranked@scaled)
length(bc.ranked@switch.point)

drugs.matrix <- bc.ranked@scaled[top.diff,]
sp <- bc.ranked@switch.point[top.diff]

drugs.matrix <- apply(X = drugs.matrix, MARGIN = 2, FUN = function(x) {x-sp})

library(scales)
drugs.matrix <- t(apply(X = drugs.matrix, MARGIN = 1, FUN = function(x){
  if (all(x >= 0)){
    scaled.row <- scales::rescale(x, to = c(0,1))
  } else if (all(x <= 0)) {
    scaled.row <- scales::rescale(x, to = c(-1,0))
  } else {
    max.positive <- max(x[x > 0])
    min.negative <- min(x[x < 0])
    row.positive <- x[x >= 0]
    row.negative <- x[x <= 0]
    
    row.positive.scaled <- rescale(row.positive, to = c(0, 1), from = c(min(row.positive), max.positive))
    row.negative.scaled <- rescale(row.negative, to = c(-1, 0), from = c(max(row.negative), min.negative))
    
    scaled.row <- x
    scaled.row[x >= 0] <- row.positive.scaled
    scaled.row[x <= 0] <- row.negative.scaled
  } 
  return(scaled.row)
}))
drugs.matrix <- t(apply(drugs.matrix, 1, function(x){round(x,2)}))

col.order <- bc.ranked@meta.data %>%
  select(TCs_res.0.3, spot.collapse) %>%
  arrange(TCs_res.0.3, spot.collapse)

col.order.spots <- col.order  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix <- drugs.matrix[,col.order.spots]

## Calculate maximum and minimum for matrix
drugs.max.matrix <- max(apply(drugs.matrix, 1, function(row) max(row)))
drugs.min.matrix <- min(apply(drugs.matrix, 1, function(row) min(row)))

# Collapsed moas for breast-SSc signature
## Import manual collapsed moas drugs
collapsed.moas <- read_tsv(file = "./data/tsv/collapsed.moas.top.differential.drugs - top.differential.drugs.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$top.diff

names.moas <- levels(factor(collapsed.moas$collapsed.MoAs))
length.moas <- length(names.moas)

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

# High-Low sensitivity top diff drugs
top.diff.df <- as.data.frame(bc.ranked@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  mutate(cluster = gsub(pattern = "TCs_res.0.3.group.TC.", replacement = "", x = cluster),
         group = gsub(pattern = "TOP-Differential-", replacement = "", x = group))

# TODO: Generate a loop (lapply or similar) for this data extraction
TC1.sensitivity <- top.diff.df %>%
  filter(cluster == 1) %>%
  mutate(cluster = NULL) %>%
  rename(TC1.sensitivity = group,
         top.diff = signature)
TC2.sensitivity <- top.diff.df %>%
  filter(cluster == 2) %>%
  mutate(cluster = NULL) %>%
  rename(TC2.sensitivity = group,
         top.diff = signature)
TC3.sensitivity <- top.diff.df %>%
  filter(cluster == 3) %>%
  mutate(cluster = NULL) %>%
  rename(TC3.sensitivity = group,
         top.diff = signature)
TC4.sensitivity <- top.diff.df %>%
  filter(cluster == 4) %>%
  mutate(cluster = NULL) %>%
  rename(TC4.sensitivity = group,
         top.diff = signature)
TC5.sensitivity <- top.diff.df %>%
  filter(cluster == 5) %>%
  mutate(cluster = NULL) %>%
  rename(TC5.sensitivity = group,
         top.diff = signature)
TC6.sensitivity <- top.diff.df %>%
  filter(cluster == 6) %>%
  mutate(cluster = NULL) %>%
  rename(TC6.sensitivity = group,
         top.diff = signature)

collapsed.moas <- collapsed.moas %>%
  left_join(TC1.sensitivity, by = join_by(top.diff)) %>%
  left_join(TC2.sensitivity, by = join_by(top.diff)) %>%
  left_join(TC3.sensitivity, by = join_by(top.diff)) %>%
  left_join(TC4.sensitivity, by = join_by(top.diff)) %>%
  left_join(TC5.sensitivity, by = join_by(top.diff)) %>%
  left_join(TC6.sensitivity, by = join_by(top.diff))

head(collapsed.moas)
col.sensitivity <- c("HighSensitivity" = "yellow",
                     "LowSensitivity" = "purple")
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

#scaled.matrix <- t(scale(t(drugs.matrix)))
## Create heatmap
heatmap.drugs <- Heatmap(
  drugs.matrix,
  #scaled.matrix,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order$TCs_res.0.3,
                                     "Cell type" = col.order$spot.collapse,
                                     col = list("TCs" = TC.colors,
                                                "Cell type" = colors.categories)),
  right_annotation = rowAnnotation(TC1.sens = collapsed.moas$TC1.sensitivity,
                                   TC2.sens = collapsed.moas$TC2.sensitivity,
                                   TC3.sens = collapsed.moas$TC3.sensitivity,
                                   TC4.sens = collapsed.moas$TC4.sensitivity,
                                   TC5.sens = collapsed.moas$TC5.sensitivity,
                                   TC6.sens = collapsed.moas$TC6.sensitivity,
                                   MoA = collapsed.moas$collapsed.MoAs,
                                   col = list(TC1.sens = col.sensitivity,
                                              TC2.sens = col.sensitivity,
                                              TC3.sens = col.sensitivity,
                                              TC4.sens = col.sensitivity,
                                              TC5.sens = col.sensitivity,
                                              TC6.sens = col.sensitivity,
                                              MoA = col.moas)),
  #left_annotation = rowAnnotation(labels = c(1:15)),
  #row_dend_reorder = TRUE,
  show_column_names = FALSE,
  column_split = col.order$TCs_res.0.3,
  #column_order = col.order,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas$preferred.drug.names),
  #row_km = 6,
  #row_order = collapsed.moas$collapsed.MoAs,
  row_split = 7,
  #show_row_dend = F,
  #row_title = NULL,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
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

###############################################################################
# HEatmap including Cancer Epithelial percentage

col.order2 <- bc.ranked@meta.data %>%
  select(TCs_res.0.3, Cancer.Epithelial,B.cells,T.cells,spot.collapse) %>%
  mutate(Lymphoid = round(B.cells + T.cells, 2)) %>%
  arrange(TCs_res.0.3,Cancer.Epithelial)

col.order.spots2 <- col.order2  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix2 <- drugs.matrix[,col.order.spots2]

# Scale color for cancer epithelial
n.cancer = length(col.order2$Cancer.Epithelial)
min_cancer = round(min(col.order2$Cancer.Epithelial), 2)
max_cancer = round(max(col.order2$Cancer.Epithelial), 2)
scale.cancer = circlize::colorRamp2(seq(min_cancer, max_cancer, length = n.cancer), hcl.colors(n.cancer,"viridis"))

## Scale color for Lymphoid
#n.lymphoid = length(col.order2$Lymphoid)
#min_lymphoid = round(min(col.order2$Lymphoid), 2)
#max_lymphoid = round(max(col.order2$Lymphoid), 2)
#scale.lymphoid = circlize::colorRamp2(seq(min_lymphoid, max_lymphoid, length = n.lymphoid), hcl.colors(n.lymphoid,"viridis"))
#
## Scale color for B cells
#n.bcells = length(col.order2$B.cells)
#min_bcells = round(min(col.order2$B.cells), 2)
#max_bcells = round(max(col.order2$B.cells), 2)
#scale.bcells = circlize::colorRamp2(seq(min_bcells, max_bcells, length = n.bcells), hcl.colors(n.bcells,"viridis"))
#
## Scale color for T cells
#n.tcells = length(col.order2$T.cells)
#min_tcells = round(min(col.order2$T.cells), 2)
#max_tcells = round(max(col.order2$T.cells), 2)
#scale.tcells = circlize::colorRamp2(seq(min_tcells, max_tcells, length = n.tcells), hcl.colors(n.tcells,"viridis"))

# Create Heatmap
heatmap.drugs.cancer <- Heatmap(
  drugs.matrix2,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order2$TCs_res.0.3,
                                     "Cancer Epithelial" = col.order2$Cancer.Epithelial,
                                     #"Lymphoid" = col.order2$Lymphoid,
                                     #"B cells" = col.order2$B.cells,
                                     #"T cells" = col.order2$T.cells,
                                     "Cell type" = col.order2$spot.collapse,
                                     col = list("TCs" = TC.colors,
                                                "Cancer Epithelial" = scale.cancer,
                                                #"Lymphoid" = scale.lymphoid,
                                                #"B cells" = scale.bcells,
                                                #"T cells" = scale.tcells,
                                                "Cell type" = colors.categories)),
  show_column_names = FALSE,
  column_split = col.order2$TCs_res.0.3,
  right_annotation = rowAnnotation(TC1.sens = collapsed.moas$TC1.sensitivity,
                                   TC2.sens = collapsed.moas$TC2.sensitivity,
                                   TC3.sens = collapsed.moas$TC3.sensitivity,
                                   TC4.sens = collapsed.moas$TC4.sensitivity,
                                   TC5.sens = collapsed.moas$TC5.sensitivity,
                                   TC6.sens = collapsed.moas$TC6.sensitivity,
                                   MoA = collapsed.moas$collapsed.MoAs,
                                   col = list(TC1.sens = col.sensitivity,
                                              TC2.sens = col.sensitivity,
                                              TC3.sens = col.sensitivity,
                                              TC4.sens = col.sensitivity,
                                              TC5.sens = col.sensitivity,
                                              TC6.sens = col.sensitivity,
                                              MoA = col.moas)),
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas$preferred.drug.names),
  row_split = 5,
  col = colorRamp2(c(drugs.min.matrix, 0, drugs.max.matrix), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.min.matrix, 0, drugs.max.matrix))
)      
heatmap.drugs.cancer
heatmap.drugs.cancer <- draw(heatmap.drugs.cancer, merge_legend = TRUE)
heatmap.drugs.cancer

dir.create(path = "./results/plots/Beyondcell_oct23_DrugRank")
png(filename = "./results/plots/Beyondcell_oct23_DrugRank/Heatmap_DrugRank95_2.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs.cancer)
dev.off()
svg(filename = "./results/plots/Beyondcell_oct23_DrugRank/Heatmap_DrugRank95_2.svg",
    width = 48,
    height = 24)
draw(heatmap.drugs.cancer)
dev.off()
