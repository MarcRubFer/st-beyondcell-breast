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

bc.TCs.Tumour <- readRDS(file = "./results/analysis/beyondcell.Tumour-TCs.rds")

# Drugs ranking
bc.ranked.Tumour <- bcRanks(bc.TCs.Tumour, idents = "TCs_res.0.3",  resm.cutoff = c(0.05,0.95))

saveRDS(object = bc.ranked.Tumour, file = "./results/analysis/bc_ranked_Tumour.rds")

# Plots 4 squares
bc4squares.Tumour <- bc4Squares(bc.ranked.Tumour, idents = "TCs_res.0.3")
bc4squares.Tumour.plots <- wrap_plots(bc4squares.Tumour)


###############################################################################

top.diff.Tumour <- as.data.frame(bc.ranked.Tumour@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

#Export differential drugs for Tumour TCs
df.top.diff.Tumour <- as.data.frame(top.diff.Tumour) 
df.top.diff.Tumour <- df.top.diff.Tumour %>%
  separate(col = top.diff.Tumour, 
           into = c("drug","data.base","db.id"), 
           sep = "_",
           remove = F)
head(df.top.diff.Tumour)
write.table(x = df.top.diff.Tumour, 
            file = "./results/tables/top.differential.drugs.TCs.Tumour.tsv", 
            sep = "\t",
            row.names = F)

# HeatMap Drugs
# Matrix of TOP-Differential Drugs
dim(bc.ranked.Tumour@normalized)
drugs.matrix.Tumour <- bc.ranked.Tumour@normalized[top.diff.Tumour,]
dim(drugs.matrix.Tumour)

# Update: for maximize differences in heatmap we will use scaled data. But 
# we need to transform data to mantein switch point info (SP)
dim(bc.ranked.Tumour@scaled)
length(bc.ranked.Tumour@switch.point)

drugs.matrix.Tumour <- bc.ranked.Tumour@scaled[top.diff.Tumour,]
sp <- bc.ranked.Tumour@switch.point[top.diff.Tumour]

drugs.matrix.Tumour <- apply(X = drugs.matrix.Tumour, MARGIN = 2, FUN = function(x) {x-sp})

library(scales)
drugs.matrix.Tumour <- t(apply(X = drugs.matrix.Tumour, MARGIN = 1, FUN = function(x){
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
drugs.matrix.Tumour <- t(apply(drugs.matrix.Tumour, 1, function(x){round(x,2)}))

col.order.Tumour <- bc.ranked.Tumour@meta.data %>%
  select(TCs_res.0.3, spot.collapse, Cancer.Epithelial) %>%
  arrange(TCs_res.0.3, spot.collapse,Cancer.Epithelial)

col.order.spots.tumour <- col.order.Tumour  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.Tumour <- drugs.matrix.Tumour[,col.order.spots.tumour]




# Collapsed moas for breast-SSc signature 
## Import manual collapsed moas TME drugs
collapsed.moas.Tumour <- read_tsv(file = "./data/tsv/top.differential.drugs.TCs.tumour.tsv")
collapsed.moas.Tumour <- as.data.frame(collapsed.moas.Tumour)
rownames(collapsed.moas.Tumour) <- collapsed.moas.Tumour$top.diff

names.moas.Tumour <- levels(factor(collapsed.moas.Tumour$collapsed.MoAs))
length.moas.Tumour <- length(names.moas.Tumour)

cols.drugs.Tumour <- c("#cd9046",
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
names(cols.drugs.Tumour) <- names.moas.Tumour

#High-LowSensitibity Differential by TCs
top.diff.Tumour.df <- as.data.frame(bc.ranked.Tumour@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  mutate(cluster = gsub(pattern = "TCs_res.0.3.group.TC.", replacement = "", x = cluster),
         group = gsub(pattern = "TOP-Differential-", replacement = "", x = group))

TC3.sensitivity <- top.diff.Tumour.df %>%
  filter(cluster == 3) %>%
  mutate(cluster = NULL) %>%
  rename(TC3.sensitivity = group,
         top.diff = signature)
TC4.sensitivity <- top.diff.Tumour.df %>%
  filter(cluster == 4) %>%
  mutate(cluster = NULL) %>%
  rename(TC4.sensitivity = group,
         top.diff = signature)
TC5.sensitivity <- top.diff.Tumour.df %>%
  filter(cluster == 5) %>%
  mutate(cluster = NULL) %>%
  rename(TC5.sensitivity = group,
         top.diff = signature)
TC6.sensitivity <- top.diff.Tumour.df %>%
  filter(cluster == 6) %>%
  mutate(cluster = NULL) %>%
  rename(TC6.sensitivity = group,
         top.diff = signature)


collapsed.moas.Tumour <- collapsed.moas.Tumour %>%
  left_join(TC3.sensitivity, by = join_by(top.diff)) %>%
  left_join(TC4.sensitivity,by = join_by(top.diff)) %>%
  left_join(TC5.sensitivity,by = join_by(top.diff)) %>%
  left_join(TC6.sensitivity,by = join_by(top.diff))

head(collapsed.moas.Tumour)
col.sensitivity <- c("HighSensitivity" = "yellow",
                     "LowSensitivity" = "purple")


## Calculate maximum and minimum for matrix
drugs.matrix.Tumour.max <- max(apply(drugs.matrix.Tumour, 1, function(row) max(row)))
drugs.matrix.Tumour.min <- min(apply(drugs.matrix.Tumour, 1, function(row) min(row)))

heatmap.drugs.Tumour <- Heatmap(
  drugs.matrix.Tumour,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order.Tumour$TCs_res.0.3,
                                     "Cell type" = col.order.Tumour$spot.collapse,
                                     col = list("TCs" = TC.colors,
                                                "Cell type" = colors.categories)),
  right_annotation = rowAnnotation(MoA = collapsed.moas.Tumour$collapsed.MoAs,
                                   TC3.sens = collapsed.moas.Tumour$TC3.sensitivity,
                                   TC4.sens = collapsed.moas.Tumour$TC4.sensitivity,
                                   TC5.sens = collapsed.moas.Tumour$TC5.sensitivity,
                                   TC6.sens = collapsed.moas.Tumour$TC6.sensitivity,
                                   col = list(MoA = cols.drugs.Tumour,
                                              TC3.sens = col.sensitivity,
                                              TC4.sens = col.sensitivity,
                                              TC5.sens = col.sensitivity,
                                              TC6.sens = col.sensitivity)),
  show_column_names = FALSE,
  column_split = col.order.Tumour$TCs_res.0.3,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas.Tumour$preferred.drug.names),
  cluster_rows = T,
  row_split = 6,
  col = colorRamp2(c(drugs.matrix.Tumour.min, 0, drugs.matrix.Tumour.max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.matrix.Tumour.min, 0, drugs.matrix.Tumour.max))
)      
heatmap.drugs.Tumour
heatmap.drugs.Tumour <- draw(heatmap.drugs.Tumour, merge_legend = TRUE)
heatmap.drugs.Tumour

png(filename = "./results/plots/TC_TME_analysis/heatmap_TME_TCs.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs.TME)
dev.off()

# HEatmap including Cancer Epithelial percentage

# col.order without spot.collapse
#drugs.matrix.Tumour2 <- bc.ranked.Tumour@normalized[top.diff.Tumour,]
col.order.Tumour2 <- bc.ranked.Tumour@meta.data %>%
  select(TCs_res.0.3, Cancer.Epithelial,B.cells,T.cells,spot.collapse) %>%
  mutate(Lymphoid = round(B.cells + T.cells, 2)) %>%
  arrange(TCs_res.0.3,Cancer.Epithelial)

col.order.spots.tumour2 <- col.order.Tumour2  %>%
  rownames_to_column("spots") %>%
  pull(spots)

drugs.matrix.Tumour2 <- drugs.matrix.Tumour[,col.order.spots.tumour2]

# Scale color for cancer epithelial
n = length(col.order.Tumour2$Cancer.Epithelial)
min_v = round(min(col.order.Tumour2$Cancer.Epithelial), 2)
max_v = round(max(col.order.Tumour2$Cancer.Epithelial), 2)
scale.cancer = circlize::colorRamp2(seq(min_v, max_v, length = n), hcl.colors(n,"viridis"))

# Scale color for Lymphoid
n.lymphoid = length(col.order.Tumour2$Lymphoid)
min_lymphoid = round(min(col.order.Tumour2$Lymphoid), 2)
max_lymphoid = round(max(col.order.Tumour2$Lymphoid), 2)
scale.lymphoid = circlize::colorRamp2(seq(min_lymphoid, max_lymphoid, length = n.lymphoid), hcl.colors(n.lymphoid,"viridis"))

# Scale color for B cells
n.bcells = length(col.order.Tumour2$B.cells)
min_bcells = round(min(col.order.Tumour2$B.cells), 2)
max_bcells = round(max(col.order.Tumour2$B.cells), 2)
scale.bcells = circlize::colorRamp2(seq(min_bcells, max_bcells, length = n.bcells), hcl.colors(n.bcells,"viridis"))

# Scale color for T cells
n.tcells = length(col.order.Tumour2$T.cells)
min_tcells = round(min(col.order.Tumour2$T.cells), 2)
max_tcells = round(max(col.order.Tumour2$T.cells), 2)
scale.tcells = circlize::colorRamp2(seq(min_tcells, max_tcells, length = n.tcells), hcl.colors(n.tcells,"viridis"))

# Draw Heatmap
heatmap.drugs.Tumour.cancerepith <- Heatmap(
  drugs.matrix.Tumour2,
  name = "bcScore",
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation("TCs" = col.order.Tumour2$TCs_res.0.3,
                                     "Cancer Epith" = col.order.Tumour2$Cancer.Epithelial,
                                     "Lymphoid" = col.order.Tumour2$Lymphoid,
                                     "B cells" = col.order.Tumour2$B.cells,
                                     "T cells" = col.order.Tumour2$T.cells,
                                     "Cell type" = col.order.Tumour2$spot.collapse,
                                     col = list("TCs" = TC.colors,
                                                "Cancer Epith" = scale.cancer,
                                                "Lymphoid" = scale.lymphoid,
                                                "B cells" = scale.bcells,
                                                "T cells" = scale.tcells,
                                                "Cell type" = colors.categories)),
  right_annotation = rowAnnotation(TC3.sens = collapsed.moas.Tumour$TC3.sensitivity,
                                   TC4.sens = collapsed.moas.Tumour$TC4.sensitivity,
                                   TC5.sens = collapsed.moas.Tumour$TC5.sensitivity,
                                   TC6.sens = collapsed.moas.Tumour$TC6.sensitivity,
                                   MoA = collapsed.moas.Tumour$collapsed.MoAs,
                                   col = list(MoA = cols.drugs.Tumour,
                                              TC3.sens = col.sensitivity,
                                              TC4.sens = col.sensitivity,
                                              TC5.sens = col.sensitivity,
                                              TC6.sens = col.sensitivity)),
  show_column_names = FALSE,
  column_split = col.order.Tumour2$TCs_res.0.3,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas.Tumour$preferred.drug.names),
  cluster_rows = T,
  row_split = 4,
  col = colorRamp2(c(drugs.matrix.Tumour.min, 0, drugs.matrix.Tumour.max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.matrix.Tumour.min, 0, drugs.matrix.Tumour.max))
)      
heatmap.drugs.Tumour.cancerepith
heatmap.drugs.Tumour.cancerepith <- draw(heatmap.drugs.Tumour.cancerepith, merge_legend = TRUE)
heatmap.drugs.Tumour.cancerepith
png(filename = "./results/plots/TC_Tumour_analysis/heatmap_Tumour_TCs_CancerEpith_scaled.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs.Tumour.cancerepith)
dev.off()
