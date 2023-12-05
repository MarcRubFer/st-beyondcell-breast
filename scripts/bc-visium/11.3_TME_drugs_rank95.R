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

# Other variables
TC.colors <- c("TC-1" = "#00b2d7",
               "TC-2" = "#e5c22f",
               "TC-3" = "#903ca2",
               "TC-4" = "#3f8741",
               "TC-5" = "#ff7b00",
               "TC-6" = "#cb5c42")

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

# Load bc object TME -TCs

bc.TCs.TME <- readRDS(file = "./results/analysis/beyondcell.TME-TCs.rds")

# Drugs ranking
bc.ranked.TME <- bcRanks(bc.TCs.TME, idents = "TCs_res.0.3",  resm.cutoff = c(0.05,0.95))

saveRDS(object = bc.ranked.TME, file = "./results/analysis/bc_ranked_TME.rds")

# Plots 4 squares
bc4squares.TME <- bc4Squares(bc.ranked.TME, idents = "TCs_res.0.3")
bc4squares.TME.plots <- wrap_plots(bc4squares.TME)


###############################################################################

top.diff.TME <- as.data.frame(bc.ranked.TME@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

#Export differential drugs for Tumour TCs
df.top.diff.TME <- as.data.frame(top.diff.TME) 
df.top.diff.TME <- df.top.diff.TME %>%
  separate(col = top.diff.TME, 
           into = c("drug","data.base","db.id"), 
           sep = "_",
           remove = F)
head(df.top.diff.TME)
write.table(x = df.top.diff.TME, 
            file = "./results/tables/top.differential.drugs.TCs.TME.tsv", 
            sep = "\t",
            row.names = F)

# HeatMap Drugs
# Matrix of TOP-Differential Drugs
dim(bc.ranked.TME@normalized)
drugs.matrix.TME <- bc.ranked.TME@normalized[top.diff.TME,]
dim(drugs.matrix.TME)

# Update: for maximize differences in heatmap we will use scaled data. But 
# we need to transform data to mantein switch point info (SP)
dim(bc.ranked.TME@scaled)
length(bc.ranked.TME@switch.point)

drugs.matrix.TME <- bc.ranked.TME@scaled[top.diff.TME,]
sp <- bc.ranked.TME@switch.point[top.diff.TME]

drugs.matrix.TME <- apply(X = drugs.matrix.TME, MARGIN = 2, FUN = function(x) {x-sp})

library(scales)
drugs.matrix.TME <- t(apply(X = drugs.matrix.TME, MARGIN = 1, FUN = function(x){
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
drugs.matrix.TME <- t(apply(drugs.matrix.TME, 1, function(x){round(x,2)}))

col.order.TME <- bc.ranked.TME@meta.data %>%
  select(TCs_res.0.3, spot.collapse, Cancer.Epithelial) %>%
  arrange(TCs_res.0.3, spot.collapse,Cancer.Epithelial)

col.order.spots.tme <- col.order.TME  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.TME <- drugs.matrix.TME[,col.order.spots.tme]


# col.order without spot.collapse
drugs.matrix.TME <- bc.ranked.TME@normalized[top.diff.TME,]
col.order.TME2 <- bc.ranked.TME@meta.data %>%
  select(TCs_res.0.3, Cancer.Epithelial,spot.collapse) %>%
  arrange(TCs_res.0.3,Cancer.Epithelial)

col.order.spots.tme2 <- col.order.TME2  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.TME2 <- drugs.matrix.TME[,col.order.spots.tme2]

# Collapsed moas for breast-SSc signature 
## Import manual collapsed moas TME drugs
collapsed.moas.TME <- read_tsv(file = "./data/tsv/top.differential.drugs.TCs.TME.tsv")
collapsed.moas.TME <- as.data.frame(collapsed.moas.TME)
rownames(collapsed.moas.TME) <- collapsed.moas.TME$top.diff

names.moas.TME <- levels(factor(collapsed.moas.TME$collapsed.MoAs))
length.moas.TME <- length(names.moas.TME)

cols.drugs.TME <- c("#cd9046",
                    "#9d5ccc",
                    "#5eb04d",
                    "#c84eb0",
                    "#b4b540",
                    "#676cc6",
                    "#ce5235",
                    "#4baf90",
                    "#d24376",
                    "#737b36",
                    "#c177ae",
                    "#6b9ad5",
                    "#b96260")
names(cols.drugs.TME) <- names.moas.TME

#High-LowSensitibity Differential by TCs
top.diff.TME.df <- as.data.frame(bc.ranked.TME@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  mutate(cluster = gsub(pattern = "TCs_res.0.3.group.TC.", replacement = "", x = cluster),
         group = gsub(pattern = "TOP-Differential-", replacement = "", x = group))

TC1.sensitivity <- top.diff.TME.df%>%
  filter(cluster == 1) %>%
  mutate(cluster = NULL) %>%
  rename(TC1.sensitivity = group,
         top.diff = signature)
TC2.sensitivity <- top.diff.TME.df %>%
  filter(cluster == 2) %>%
  mutate(cluster = NULL) %>%
  rename(TC2.sensitivity = group,
         top.diff = signature)

collapsed.moas.TME <- collapsed.moas.TME %>%
  left_join(TC1.sensitivity, by = join_by(top.diff)) %>%
  left_join(TC2.sensitivity,by = join_by(top.diff)) 

head(collapsed.moas.TME)
col.sensitivity <- c("HighSensitivity" = "yellow",
                     "LowSensitivity" = "purple")

## Calculate maximum and minimum for matrix
drugs.matrix.TME.max <- max(apply(drugs.matrix.TME, 1, function(row) max(row)))
drugs.matrix.TME.min <- min(apply(drugs.matrix.TME, 1, function(row) min(row)))

heatmap.drugs.TME <- Heatmap(
  drugs.matrix.TME,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order.TME$TCs_res.0.3,
                                     "Cell type" = col.order.TME$spot.collapse,
                                     col = list("TCs" = TC.colors[1:2],
                                                "Cell type" = colors.categories)),
  right_annotation = rowAnnotation(MoA = collapsed.moas.TME$collapsed.MoAs,
                                   TC1.sens = collapsed.moas.TME$TC1.sensitivity,
                                   TC4.sens = collapsed.moas.TME$TC2.sensitivity,
                                   col = list(MoA = cols.drugs.TME,
                                              TC1.sens = col.sensitivity,
                                              TC2.sens = col.sensitivity)),
  show_column_names = FALSE,
  column_split = col.order.TME$TCs_res.0.3,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas.TME$preferred.drug.names),
  cluster_rows = T,
  row_split = 6,
  col = colorRamp2(c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max))
)      
heatmap.drugs.TME
heatmap.drugs.TME <- draw(heatmap.drugs.TME, merge_legend = TRUE)
heatmap.drugs.TME
png(filename = "./results/plots/TC_TME_analysis/heatmap_TME_TCs.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs.TME)
dev.off()

# HEatmap including Cancer Epithelial percentage

# col.order without spot.collapse
#drugs.matrix.TME <- bc.ranked.TME@normalized[top.diff.TME,]
col.order.TME2 <- bc.ranked.TME@meta.data %>%
  select(TCs_res.0.3, Cancer.Epithelial, B.cells, T.cells, spot.collapse) %>%
  mutate(Lymphoid = round(B.cells + T.cells, 2)) %>%
  #arrange(TCs_res.0.3,Lymphoid)
  arrange(TCs_res.0.3,Cancer.Epithelial)

col.order.spots.tme2 <- col.order.TME2  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.TME2 <- drugs.matrix.TME[,col.order.spots.tme2]

# Scale color for cancer epithelial
n.cancer = length(col.order.TME2$Cancer.Epithelial)
min_cancer = round(min(col.order.TME2$Cancer.Epithelial), 2)
max_cancer = round(max(col.order.TME2$Cancer.Epithelial), 2)
scale.cancer = circlize::colorRamp2(seq(min_cancer, max_cancer, length = n.cancer), hcl.colors(n.cancer,"viridis"))

# Scale color for Lymphoid
n.lymphoid = length(col.order.TME2$Lymphoid)
min_lymphoid = round(min(col.order.TME2$Lymphoid), 2)
max_lymphoid = round(max(col.order.TME2$Lymphoid), 2)
scale.lymphoid = circlize::colorRamp2(seq(min_lymphoid, max_lymphoid, length = n.lymphoid), hcl.colors(n.lymphoid,"viridis"))

# Scale color for B cells
n.bcells = length(col.order.TME2$B.cells)
min_bcells = round(min(col.order.TME2$B.cells), 2)
max_bcells = round(max(col.order.TME2$B.cells), 2)
scale.bcells = circlize::colorRamp2(seq(min_bcells, max_bcells, length = n.bcells), hcl.colors(n.bcells,"viridis"))

# Scale color for T cells
n.tcells = length(col.order.TME2$T.cells)
min_tcells = round(min(col.order.TME2$T.cells), 2)
max_tcells = round(max(col.order.TME2$T.cells), 2)
scale.tcells = circlize::colorRamp2(seq(min_tcells, max_tcells, length = n.tcells), hcl.colors(n.tcells,"viridis"))

# Draw Heatmap
heatmap.drugs.TME.cancerepith <- Heatmap(
  drugs.matrix.TME2,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order.TME2$TCs_res.0.3,
                                     "Cancer Epith" = col.order.TME2$Cancer.Epithelial,
                                     "Lymphoid" = col.order.TME2$Lymphoid,
                                     "B cells" = col.order.TME2$B.cells,
                                     "T cells" = col.order.TME2$T.cells,
                                     "Cell type" = col.order.TME2$spot.collapse,
                                     col = list("TCs" = TC.colors[1:2],
                                                "Cell type" = colors.categories,
                                                "Cancer Epith" = scale.cancer,
                                                "Lymphoid" = scale.lymphoid,
                                                "B cells" = scale.bcells,
                                                "T cells" = scale.tcells)),
  right_annotation = rowAnnotation(TC1.sens = collapsed.moas.TME$TC1.sensitivity,
                                   TC2.sens = collapsed.moas.TME$TC2.sensitivity,
                                   MoA = collapsed.moas.TME$collapsed.MoAs,
                                   col = list(MoA = cols.drugs.TME,
                                              TC1.sens = col.sensitivity,
                                              TC2.sens = col.sensitivity)),
  show_column_names = FALSE,
  column_split = col.order.TME2$TCs_res.0.3,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas.TME$preferred.drug.names),
  cluster_rows = T,
  row_split =  2,
  col = colorRamp2(c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max))
)      
heatmap.drugs.TME.cancerepith
heatmap.drugs.TME.cancerepith <- draw(heatmap.drugs.TME.cancerepith, merge_legend = TRUE)
heatmap.drugs.TME.cancerepith
png(filename = "./results/plots/TC_TME_analysis/heatmap_TME_TCs_CancerEpith.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs.TME.cancerepith)
dev.off()

