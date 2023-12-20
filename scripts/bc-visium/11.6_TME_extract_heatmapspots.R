
drugs.matrix.import <- read_tsv(file = "./results/tables/drugs_matrix_TME_CancerEp_ordered.tsv")
drugs.matrix.import <- as.data.frame(drugs.matrix.import) 
rownames(drugs.matrix.import) <- drugs.matrix.import$sigs_diff

drugs.matrix.TME2 <- drugs.matrix.import %>%
  select(-sigs_diff)
drugs.matrix.TME2 <- as.matrix(drugs.matrix.TME2)
dim(drugs.matrix.TME2)
# Draw Heatmap
heatmap.drugs.TME.cancerepith2 <- Heatmap(
  drugs.matrix.TME2,
  name = "bcScore",
  cluster_columns = T,
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


bcClusters(bc.ranked.TME, 
           idents = "spot.collapse", 
           spatial = T, mfrow = c(1,2), 
           cells.highlight = list("cluster2" = cluster2))
bcClusters(bc.ranked.TME, 
           idents = "spot.collapse", 
           spatial = T, mfrow = c(1,2))


seurat.TME <- readRDS("./results/analysis/seuratobj.TME-TCs.rds")

TLS.genes <- c("CCL19", "CCL21","CXCL13", "CCR7", "SELL","LAMP3","CXCR4","CD86","BCL6")

rownames(cluster.colum) <- cluster.colum$spot

col.order.TLS <- cluster.colum %>%
  arrange(desc(cluster))
col.order.TLS.spots <- col.order.TLS %>%
  select(spot) %>%
  pull()


TLS.norm.data <- seurat.TME@assays$SCT@data[TLS.genes,]
TLS.norm.data <- as.matrix(TLS.norm.data)
TLS.norm.data <- TLS.norm.data[,col.order.TLS.spots]
# Scale color for heatmap
n.cluster = length(col.order.TLS$cluster)
min_cluster = round(min(TLS.norm.data), 2)
max_cluster = round(max(TLS.norm.data), 2)
scale.expression = circlize::colorRamp2(seq(min_cluster, max_cluster, length = n.cluster), (hcl.colors(n.cancer,"PuRd")))
scale.expression.rev = circlize::colorRamp2(seq(min_cluster, max_cluster, length = n.cluster), rev(hcl.colors(n.cancer,"PuRd")))
col_fun = colorRamp2(c(min_cluster,max_cluster), c("purple", "yellow"))

ht.TLS <- Heatmap(matrix = TLS.norm.data,
        #col = scale.expression,
        col = scale.expression.rev,
        #col = col_fun,
        show_column_names = F,
        cluster_columns = F,
        column_gap = unit(3, "mm"),
        column_split = col.order.TLS$cluster,
        top_annotation = HeatmapAnnotation("cluster" = col.order.TLS$cluster,
                                           show_annotation_name = T,
                                           annotation_name_gp = gpar(fontsize = 5)),
        column_names_gp = gpar(fontsize = 3))
ht.TLS

seurat.TME <- readRDS("./results/analysis/seuratobj.TME-TCs.rds")
gs.TLS <- GenerateGenesets(x = "./data/gmts/TLS_sig.gmt")

bc.TLS <- bcScore(seurat.TME,
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


##
data <- as.data.frame(t(bc.TLS.recomputed@normalized)) %>%
  rownames_to_column("spot")

spot.data <- data %>%
  select(spot) %>%
  pull()
cluster.selec <- cluster.colum %>%
  filter(spot %in% spot.data)

rownames(cluster.selec) <- cluster.selec$spot
cluster.add <- cluster.selec %>%
  mutate(spot = NULL)

bc.TLS.recomputed <- bcAddMetadata(bc.TLS.recomputed, metadata = cluster.add)
bc.TLS.recomputed@meta.data

spatial.bcScore.TLS <- bcSignatures(bc = bc.TLS.recomputed, UMAP = "Seurat", spatial = T, mfrow = c(1,2), signatures = list(values = "TLS_CABRITA"))

spatial.Bcells <- SpatialFeaturePlot(seurat.TME, features = "B.cells", ncol = 1)
spatial.Tcells <- SpatialFeaturePlot(seurat.TME, features = "T.cells", ncol = 1)

spatial.bcScore.TLS | spatial.Bcells | spatial.Tcells

data.scales <- as.data.frame(t(bc.TLS.recomputed@scaled)) %>%
  mutate(TLS.spots = case_when(TLS_CABRITA >= 0.4 ~ "TLS",
                               TRUE ~ NA))

head(data.scales)

bc.TLS.spot <- bcAddMetadata(bc.TLS.recomputed, metadata = data.scales)
head(bc.TLS.spot@meta.data)

bcClusters(bc = bc.TLS.spot,
           idents = "TLS.spots",
           spatial = T,
           #factor.col = F,
           mfrow = c(1,2))

data.merged <- left_join(data, cluster.selec, by = "spot")
data.boxplot <- ggplot(data.merged, aes(x=cluster, y=TLS_CABRITA)) +
  geom_boxplot()

library(ggstatsplot)

stats.plot <-ggbetweenstats(data = data.merged,
               x = cluster,
               y = TLS_CABRITA)
extract_stats(stats.plot)

data.boxplot | stats.plot

seurat.tcs <- readRDS("./results/analysis/seuratobj.therapeutic.clusters.rds")

TLS.spots <- data.scales %>%
  select(TLS.spots)

seurat.TLS <- AddMetaData(seurat.tcs, metadata = TLS.spots)
head(seurat.TLS@meta.data)
SpatialDimPlot(seurat.TLS, group.by = "TLS.spots")

library(semla)
# Semla object
semlaobj <- UpdateSeuratForSemla(seurat.TLS, 
                                 image_type = "tissue_lowres",
                                 verbose = T)

# Radial distances
semlaobj <- RadialDistance(semlaobj, 
                           column_name = "TLS.spots",
                           selected_groups = "TLS")

dual.colors <- c(
  "TUMOUR" = "#c4534e",
  "TME" = "#098cf0")
MapLabels(semlaobj, 
          column_name = "TLS.spots", 
          override_plot_dims = TRUE, 
          image_use = NULL, 
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
MapFeatures(semlaobj, 
            features = "r_dist_TLS_scaled", 
            center_zero = TRUE, 
            pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

bc.allspots <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")
# Transposed enriched matrix
enrich.matrix <- t(bc.allspots@normalized)

# 
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

ggplot(results.top.diff.filtered, aes(x=corr, y=reorder(preferred.drug.names, corr))) +
  geom_bar(aes(fill = corr), stat = "identity") +
  scale_fill_gradient2(limits = c(-1,1)) +
  xlim(-1,1) +
  labs(fill = "Pearson's correlation") +
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank())