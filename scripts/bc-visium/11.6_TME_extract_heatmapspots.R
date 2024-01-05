rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
library("scales")


out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

set.seed(1)

# Functions
scale.heatmap <- function(data, cell.type, palette = "viridis") {
  n.total = length(data[[cell.type]])
  min.scale = round(min(data[[cell.type]]), 2)
  max.scale = round(max(data[[cell.type]]), 2)
  scale = circlize::colorRamp2(seq(min.scale, max.scale, length = n.total), hcl.colors(n.total, palette))
  return(scale)
}

# Load matrix
drugs.matrix.import <- read_tsv(file = "./results/tables/drugs_matrix_TME_CancerEp_ordered.tsv")
drugs.matrix.import <- as.data.frame(drugs.matrix.import) 
rownames(drugs.matrix.import) <- drugs.matrix.import$sigs_diff

drugs.matrix.TME2 <- drugs.matrix.import %>%
  select(-sigs_diff)
drugs.matrix.TME2 <- as.matrix(drugs.matrix.TME2)
dim(drugs.matrix.TME2)

## Calculate maximum and minimum for matrix
drugs.matrix.TME.max <- max(apply(drugs.matrix.TME2, 1, function(row) max(row)))
drugs.matrix.TME.min <- min(apply(drugs.matrix.TME2, 1, function(row) min(row)))

# Load collapsed moas-sensitivity
collapsed.moas.TME <- read_tsv("./results/tables/top.differential.drugs.TCs.TME_sensitivity.tsv")

# Load column ordered data.frame
col.order.TME2 <- read_tsv("./results/tables/dataframe_column_ordered_TME.tsv")

# Scales
scale.cancer <- scale.heatmap(data = col.order.TME2, cell.type = "Cancer.Epithelial")
scale.bcells <- scale.heatmap(data = col.order.TME2, cell.type = "B.cells")
scale.tcells <- scale.heatmap(data = col.order.TME2, cell.type = "T.cells")
scale.cafs <- scale.heatmap(data = col.order.TME2, cell.type = "CAFs")
scale.endothelial <- scale.heatmap(data = col.order.TME2, cell.type = "Endothelial")
scale.myeloid <- scale.heatmap(data = col.order.TME2, cell.type = "Myeloid")

# Colors
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

col.sensitivity <- c("HighSensitivity" = "yellow",
                     "LowSensitivity" = "purple")

# Draw Heatmap
heatmap.drugs.TME.cancerepith2 <- Heatmap(
  drugs.matrix.TME2,
  name = "bcScore",
  cluster_columns = T,
  top_annotation = HeatmapAnnotation("TCs" = col.order.TME2$TCs_res.0.3,
                                     "Cancer Epith" = col.order.TME2$Cancer.Epithelial,
                                     #"Lymphoid" = col.order.TME2$Lymphoid,
                                     "B cells" = col.order.TME2$B.cells,
                                     "T cells" = col.order.TME2$T.cells,
                                     "CAFs" = col.order.TME2$CAFs,
                                     "Endothelial" = col.order.TME2$Endothelial,
                                     "Myeloid" = col.order.TME2$Myeloid,
                                     "Cell type" = col.order.TME2$spot.collapse,
                                     col = list("TCs" = TC.colors[1:2],
                                                "Cell type" = colors.categories,
                                                "Cancer Epith" = scale.cancer,
                                                #"Lymphoid" = scale.lymphoid,
                                                "B cells" = scale.bcells,
                                                "T cells" = scale.tcells,
                                                "CAFs" = scale.cafs,
                                                "Endothelial" = scale.endothelial,
                                                "Myeloid" = scale.myeloid)),
  right_annotation = rowAnnotation(TC1.sens = collapsed.moas.TME$TC1.sensitivity,
                                   TC2.sens = collapsed.moas.TME$TC2.sensitivity,
                                   MoA = collapsed.moas.TME$collapsed.MoAs,
                                   col = list(MoA = cols.drugs.TME,
                                              TC1.sens = col.sensitivity,
                                              TC2.sens = col.sensitivity)),
  show_column_names = FALSE,
  #column_split = col.order.TME2$TCs_res.0.3,
  column_km = 4,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas.TME$preferred.drug.names),
  cluster_rows = T,
  column_gap = unit(5,"mm"),
  row_split =  2,
  col = colorRamp2(c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max))
)      
heatmap.drugs.TME.cancerepith2
png(filename = "./results/plots/TC_TME_analysis/heatmap_TME_TCs_cluster_columns.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(heatmap.drugs.TME.cancerepith2)
dev.off()
svg(filename = "./results/plots/TC_TME_analysis/heatmap_TME_TCs_cluster_columns.svg",
    width = 19,
    height = 10)
draw(heatmap.drugs.TME.cancerepith2)
dev.off()

ht = draw(heatmap.drugs.TME.cancerepith2) 
c.order.list <- column_order(ht)
c.dend <- column_dend(ht)
dev.off()
lapply(c.order.list, function(x) length(x))

# Extract spotID for each cluster and create a data frame
c.order.list <- c.order.list[order(names(c.order.list))]
cluster.colum <- lapply(seq_along(c.order.list), FUN = function(index) {
  cluster <- c.order.list[[index]]
  spot <- colnames(drugs.matrix.TME2)[cluster]
  spot <- as.data.frame(spot)
  spot <- spot %>%
    mutate(cluster = paste0("cluster",index))
})
cluster.colum <- do.call(rbind, cluster.colum)

cluster2 <- cluster.colum %>%
  filter(cluster == "cluster2") %>%
  select(spot) %>%
  pull()

bc.ranked.TME <- readRDS("./results/analysis/bc_ranked_TME.rds")
cluster2.position <- bcClusters(bc.ranked.TME, 
                                idents = "spot.collapse", 
                                spatial = T, mfrow = c(1,2), 
                                cells.highlight = list("cluster2" = cluster2))
ggsave(filename = "cluster2_position.svg",
       plot = cluster2.position,
       path = "./results/plots/TC_TME_analysis/")


bcells.prop <- bcClusters(bc.ranked.TME, 
                          idents = "B.cells", 
                          spatial = T, mfrow = c(2,1),
                          factor.col = F)

tcells.prop <- bcClusters(bc.ranked.TME, 
                          idents = "T.cells", 
                          spatial = T, mfrow = c(2,1),
                          factor.col = F)

patch.proportions <- bcells.prop | tcells.prop
ggsave(filename = "patch_proportions.svg",
       plot = patch.proportions,
       path = "./results/plots/TC_TME_analysis/")

# Analysis of tertiary-lymphocyte structure (TLS) signature

seurat.tcs <- readRDS("./results/analysis/backup/seuratobj.therapeutic.clusters.rds")

# Signature from Cabrita et al. 2020 (DOI: https://doi.org/10.1038/s41586-019-1914-8)
TLS.genes <- c("CCL19", "CCL21","CXCL13", "CCR7", "SELL","LAMP3","CXCR4","CD86","BCL6")
TLS.expression <- SpatialFeaturePlot(seurat.tcs, features = TLS.genes, ncol = 8, combine = F)


colorscale <- function(n = 100, rev = TRUE) {
  colors <- RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  if (rev) {
    colors <- colors |> 
      rev()
  }
  colors <- colors |>
    grDevices::colorRampPalette()
  return(colors(n = n))
}

TLS.expression.scaled <- lapply(X = seq_along(TLS.expression), FUN = function(index){
  TLS.expression[[index]] +
    scale_fill_gradientn(colours = colorscale(n = 100),
                         limits = c(0,3))
})
TLS.expression.wrapped <- wrap_plots(TLS.expression.scaled, ncol = 2)

ggsave(filename = "patch_TLS_expression.svg",
       plot = TLS.expression.wrapped,
       path = "./results/plots/TC_TME_analysis/",
       width = 12.1,
       height = 23.5)



# Calculate beyondcell score for TLS signature
gs.TLS <- GenerateGenesets(x = "./data/gmts/TLS_sig.gmt")
bc.TLS <- bcScore(seurat.tcs,
                  gs = gs.TLS,
                  expr.thres = 0.1)
# Number of NAs
n.NA.TLS <- data.frame(nNAs = colSums(is.na(bc.TLS@normalized)),
                   row.names = colnames(bc.TLS@normalized))
bc.TLS <- bcAddMetadata(bc.TLS, n.NA.TLS)
bc.TLS@meta.data
# Filter out spots with a high percentage of NAs
bc.TLS <- bcSubset(bc.TLS, nan.cells = 0.95)
# Replace NAs by 0s
bc.TLS@normalized[is.na(bc.TLS@normalized)] <- 0
bc.TLS.recomputed <- bcRecompute(bc.TLS, slot = "normalized")

# Plot signature in all spots
spatial.bcScore.TLS <- bcSignatures(bc = bc.TLS.recomputed, UMAP = "Seurat", spatial = T, mfrow = c(1,2), signatures = list(values = "TLS_CABRITA"))
ggsave(filename = "bcScore_TLS_allspots.svg",
       spatial.bcScore.TLS,
       path = "./results/plots/TC_TME_analysis/")

##
TLS.spots <- as.data.frame(t(bc.TLS.recomputed@scaled)) %>%
  mutate(TLS.spots = case_when(TLS_CABRITA >= 0.4 ~ "TLS",
                               TRUE ~ NA))

head(TLS.spots)
head(bc.TLS.recomputed@meta.data)
bc.TLS.spot <- bcAddMetadata(bc.TLS.recomputed, metadata = TLS.spots)
head(bc.TLS.spot@meta.data)

bcClusters(bc = bc.TLS.spot,
           idents = "TLS.spots",
           spatial = T,
           #factor.col = F,
           mfrow = c(1,2))


seurat.tcs <- readRDS("./results/analysis/seuratobj.therapeutic.clusters.rds")



seurat.TLS <- AddMetaData(seurat.tcs, metadata = TLS.spots)
head(seurat.TLS@meta.data)
SpatialDimPlot(seurat.TLS, group.by = "TLS.spots")

library("semla")
# Semla object
semlaobj <- UpdateSeuratForSemla(seurat.TLS, 
                                 image_type = "tissue_lowres",
                                 verbose = T)

semlaobj <- LoadImages(semlaobj)
semlaobj <- DisconnectRegions(semlaobj, column_name = "TLS.spots", selected_groups = "TLS")
semlaobj@tools
MapLabels(semlaobj, 
          column_name = "TLS_split", 
          override_plot_dims = TRUE, 
          image_use = "raw", 
          drop_na = TRUE, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

semlaobj$TLS_split2 <- case_when(semlaobj$TLS_split == "S1_region1" | semlaobj$TLS_split == "S2_region1" ~ "TLS",
                                 TRUE ~ NA)

# Radial distances
semlaobj <- RadialDistance(semlaobj, 
                           column_name = "TLS_split2",
                           selected_groups = "TLS")

MapLabels(semlaobj, 
          column_name = "TLS_split2", 
          override_plot_dims = TRUE, 
          image_use = "raw", 
          drop_na = TRUE, 
          pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  #scale_fill_manual(values = dual.colors) &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
semlaobj@meta.data

MapFeatures(semlaobj, 
            features = "r_dist_TLS", 
            center_zero = TRUE, 
            pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

semlaobj$r_dist_TLS_scaled <- sign(semlaobj$r_dist_TLS)*sqrt(abs(semlaobj$r_dist_TLS))
r_dist_TLS_scaled <- MapFeatures(semlaobj, 
            features = "r_dist_TLS_scaled", 
            center_zero = TRUE, 
            pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
ggsave(filename = "r_dist_TLS_scaled.svg",
       plot = r_dist_TLS_scaled,
       path = "./results/plots/TC_TME_analysis/")


bc.allspots <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")
# Transposed enriched matrix
enrich.matrix <- t(bc.allspots@normalized)

# Vector of radial distance
radial.dist.tumour <- semlaobj@meta.data$r_dist_TLS

# Realizar la prueba de correlación de Pearson entre cada columna de datos y el vector de distancia
results.cor.pearson <- apply(enrich.matrix, 2, function(col) cor.test(col, radial.dist.tumour, method = "pearson"))

results.df <- data.frame(
  corr = sapply(results.cor.pearson, function(res) res$estimate),
  p.value = sapply(results.cor.pearson, function(res) res$p.value)
)

collapsed.moas <- read_tsv(file = "./data/tsv/collapsed.moas.top.differential.drugs - top.differential.drugs.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$top.diff

subset.moas <- collapsed.moas %>%
  select(top.diff, preferred.drug.names, collapsed.MoAs)

top.diff <- collapsed.moas$top.diff

results.top.diff <- results.df[top.diff,]
results.top.diff <- results.top.diff %>%
  rownames_to_column("top.diff") %>%
  mutate(top.diff = gsub("\\.cor$", "", top.diff)) 

results.top.diff <- right_join(y=results.top.diff, x=subset.moas, by = "top.diff")

results.top.diff.filtered <- results.top.diff %>%
  filter(p.value < 0.05) %>%
  arrange(corr)

corr.plot <- ggplot(results.top.diff.filtered, aes(x=corr, y=reorder(preferred.drug.names, corr))) +
  geom_bar(aes(fill = corr), stat = "identity") +
  scale_fill_gradient2(limits = c(-1,1)) +
  xlim(-1,1) +
  labs(fill = "Pearson's correlation") +
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank())
ggsave(filename = "corr.plot.svg",
       plot = corr.plot,
       path = "./results/plots/TC_TME_analysis/")

YM155 <- bcSignatures(bc.allspots, spatial = T, mfrow = c(2,2), signatures = list(values = c("YM-155_PRISM_K76703230")))
ggsave(filename = "YM155.svg",
       plot = YM155,
       path = "./results/plots/TC_TME_analysis/")
AM580 <- bcSignatures(bc.allspots, spatial = T, mfrow = c(2,2), signatures = list(values = c("AM-580_PRISM_K06854232")))
ggsave(filename = "AM580.svg",
       plot = AM580,
       path = "./results/plots/TC_TME_analysis/")
