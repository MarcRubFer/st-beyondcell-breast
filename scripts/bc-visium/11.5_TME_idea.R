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

# Load beyondcell ranked TME-TCs
seurat.TME <- readRDS("./results/analysis/seuratobj.TME-TCs.rds")
bc.ranked.TME <- readRDS("./results/analysis/bc_ranked_TME.rds")

# Extract high/low-differential signatures
top.diff.TME <- as.data.frame(bc.ranked.TME@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

# Re-scale matrix
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

# Order matrix by Cancer epithelial deconvoluted proportion
col.order.TME2 <- bc.ranked.TME@meta.data %>%
  select(TCs_res.0.3, Cancer.Epithelial, B.cells, T.cells, spot.collapse) %>%
  mutate(Lymphoid = round(B.cells + T.cells, 2)) %>%
  #arrange(TCs_res.0.3,Lymphoid)
  arrange(TCs_res.0.3,Cancer.Epithelial)

# Extract ordered spots
col.order.spots.tme2 <- col.order.TME2  %>%
  rownames_to_column("spots") %>%
  pull(spots)

# Order matrix by ordered spots
drugs.matrix.TME2 <- drugs.matrix.TME[,col.order.spots.tme2]
dim(drugs.matrix.TME2)

write.table(drugs.matrix.TME2, file = "./results/tables/drugs_matrix_TME_CancerEp_ordered.tsv", sep = "\t")


metadata <- bc.ranked.TME@meta.data %>%
  select(spot.collapse, TCs_res.0.3) %>%
  rownames_to_column("spot")
head(metadata)

enrichment <- as.data.frame(t(drugs.matrix.TME2)) %>%
  rownames_to_column("spot")
head(enrichment)

data.raw <- right_join(metadata,enrichment, by = "spot")
head(data.raw)

data <- data.raw %>%
  pivot_longer(cols = -c(spot,spot.collapse,TCs_res.0.3), names_to = "signature", values_to = "enrichment")
head(data)

sigs.UP <- data %>%
  filter(TCs_res.0.3 == "TC-2",
         enrichment > 0) %>%
  select(signature) %>%
  unique() %>%
  pull()

sigs.DOWN <- data %>%
  filter(TCs_res.0.3 == "TC-2",
         enrichment < 0) %>%
  select(signature) %>%
  unique() %>%
  pull()

violin.UP <- ggplot(subset(data, signature %in% sigs.UP), aes(x = TCs_res.0.3, y = enrichment)) + 
  geom_jitter(alpha = 0.3) + 
  geom_violin(aes(fill = spot.collapse)) +
  scale_fill_manual(values = colors.categories) +
  ggtitle("Sensible signatures") +
  theme_bw()


violin.DOWN <- ggplot(subset(data, signature %in% sigs.DOWN), aes(x = TCs_res.0.3, y = enrichment))  + 
  geom_jitter(alpha = 0.3) + 
  geom_violin(aes(fill = spot.collapse)) +
  scale_fill_manual(values = colors.categories) +
  ggtitle("Unsensible signatures") +
  theme_bw()

(violin.UP | violin.DOWN) + plot_layout(guides = "collect")


TC1.Down.in.UP <- data %>%
  filter(signature %in% sigs.UP,
         TCs_res.0.3 == "TC-1",
         enrichment < 0) %>%
  select(spot) %>%
  unique() %>%
  pull()

TC1.UP.in.DOWN <- data %>%
  filter(signature %in% sigs.DOWN,
         TCs_res.0.3 == "TC-1",
         enrichment > 0) %>%
  select(spot) %>%
  unique() %>%
  pull()

categories.down <- ggplot(subset(data.raw, spot %in% TC1.Down.in.UP), aes(x = spot.collapse)) +
  geom_bar(aes(fill = spot.collapse)) +
  ylim(NA,200) +
  scale_fill_manual(values = colors.categories) +
  ggtitle(label = "Categories in spots down in sigs UP") + 
  theme_bw()
categories.up <- ggplot(subset(data.raw, spot %in% TC1.UP.in.DOWN), aes(x = spot.collapse)) +
  geom_bar(aes(fill = spot.collapse)) +
  ylim(NA,200) +
  scale_fill_manual(values = colors.categories) +
  ggtitle(label = "Categories in spots up in sigs DOWN") + 
  theme_bw()

categories.figures <- (categories.down | categories.up) + plot_layout(guides = "collect")
categories.figures
ggsave(filename = "categories_spots.png",
       plot = categories.figures,
       path = "./results/plots/TC_TME_analysis/")

df.cat.global <- as.data.frame(table(seurat.TME@meta.data$spot.collapse)) %>%
  filter(Var1 != "TUMOUR")
table.cat.global
union.spots <- union(TC1.Down.in.UP, TC1.UP.in.DOWN)
subset.union <- data.raw %>%
  filter(spot %in% union.spots)
table.subset <- as.data.frame(table(subset.union$spot.collapse))
merge.df <- merge(df.cat.global,table.subset, by = "Var1")
merge.df <- merge.df %>%
  rename(category = Var1,
         n.total = Freq.x,
         n.split = Freq.y) %>%
  mutate(prop = round(n.split / n.total,2),
         perc = round(prop * 100, 2))

write.table(merge.df,
            file = "./results/tables/proportions_categories_inspot.tsv",
            sep = "\t")



seurat.TME <- readRDS("./results/analysis/seuratobj.TME-TCs.rds")
SpatialDimPlot(seurat.TME, cells.highlight = list(TC1.Down.in.UP, TC1.UP.in.DOWN), cols.highlight = c("blue","red", "grey50"))

coords.DOWN.in.UP <- seurat.TME@meta.data %>%
  rownames_to_column("spot") %>%
  filter(spot %in% TC1.Down.in.UP)
coords.DOWN.in.UP.plot <- ggplot(coords.DOWN.in.UP, aes(x = corr_y, y = corr_x)) +
  geom_point(aes(col = spot.collapse)) +
  scale_y_reverse() +
  ggtitle(label = "Spots DOWN-enrich in UP-sigs", subtitle = "At least in one sig")

coords.UP.in.DOWN <- seurat.TME@meta.data %>%
  rownames_to_column("spot") %>%
  filter(spot %in% TC1.UP.in.DOWN)
coords.UP.in.DOWN.plot <- ggplot(coords.UP.in.DOWN, aes(x = corr_y, y = corr_x)) +
  geom_point(aes(col = spot.collapse)) +
  scale_y_reverse() +
  ggtitle(label = "Spots UP-enrich in DOWN-sigs", subtitle = "At least in one sig")

intersect.spots <- intersect(TC1.Down.in.UP,TC1.UP.in.DOWN)
coords.intersect <- seurat.TME@meta.data %>%
  rownames_to_column("spot") %>%
  filter(spot %in% intersect.spots)
coords.intersect.plot <- ggplot(coords.intersect, aes(x = corr_y, y = corr_x)) +
  geom_point(aes(col = spot.collapse)) +
  scale_y_reverse() +
  ggtitle(label = "Intersect spots UP_DOWN", subtitle = "At least in one sig")

panel <- (coords.DOWN.in.UP.plot / plot_spacer() / coords.UP.in.DOWN.plot) + plot_layout(heights = c(1,0.3,1))
figure.spot <- (panel | coords.intersect.plot | SpatialDimPlot(seurat.TME,group.by = "spot.collapse", ncol = 1)) + plot_layout(widths = c(2,4,1))

ggsave(filename = "Coords_spots_TC1.png",
       plot = figure.spot,
       path = "./results/plots/TC_TME_analysis/")

cell.type.prop <- coords.intersect %>%
  select(spot, B.cells:T.cells)
head(cell.type.prop)
cell.type.prop <- cell.type.prop %>%
  pivot_longer(cols = -spot, names_to = "cell.type", values_to = "proportion")
head(cell.type.prop)
cell.type.prop.plot <- ggplot(cell.type.prop, aes(x = cell.type, y = proportion)) +
  geom_boxplot()
ggsave(filename = "cell_type_proportion.png",
       plot = cell.type.prop.plot,
       path = "./results/plots/TC_TME_analysis/")


TCs.recluster <- seurat.TME@meta.data %>%
  select(corr_x, corr_y,TCs_res.0.3) %>%
  mutate(TCs_recat = case_when(TCs_res.0.3 == "TC-2" ~ "TC-2",
                               corr_x %in% seq(4500,6100) & corr_y %in% seq(0,2500) & rownames(.) %in% intersect.spots ~ "TC-1.2",
                               corr_x %in% seq(4000,6000) & corr_y > 10000 & rownames(.) %in% intersect.spots ~ "TC-1.3",
                               TRUE ~ "TC-1.1"))

a <- ggplot(TCs.recluster, aes(x = corr_y, y = corr_x)) +
  geom_point(aes(col = TCs_recat)) +
  scale_y_reverse()

coords.intersect.plot | a


TLS.genes <- c("CCL19", "CCL21","CXCL13", "CCR7", "SELL","LAMP3","CXCR4","CD86","BCL6")
col.order.TLS <- TCs.recluster %>%
  arrange(TCs_recat) 
col.order.TLS.spots <- col.order.TLS%>%
  rownames_to_column("spots") %>%
  pull(spots)

TLS.norm.data <- TLS.norm.data[,col.order.TLS.spots]

# Scale color for cancer epithelial
n.cancer = length(col.order.TLS$TCs_recat)
min_cancer = round(min(TLS.norm.data), 2)
max_cancer = round(max(TLS.norm.data), 2)
scale.cancer = circlize::colorRamp2(seq(min_cancer, max_cancer, length = n.cancer), rev(hcl.colors(n.cancer,"OrRd")), reverse = F)
Heatmap(matrix = TLS.norm.data,
        col = scale.cancer,
        cluster_columns = F,
        show_column_names = F,
        top_annotation = HeatmapAnnotation("TC" = col.order.TLS$TCs_recat),
        column_split = col.order.TLS$TCs_recat,
        column_gap = unit(c(4), "mm"))
