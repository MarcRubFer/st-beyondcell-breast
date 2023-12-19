
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

cluster1 <- cluster.colum %>%
  filter(cluster == "cluster1") %>%
  select(spot) %>%
  pull()

cluster2 <- cluster.colum %>%
  filter(cluster == "cluster2") %>%
  select(spot) %>%
  pull()

cluster3 <- cluster.colum %>%
  filter(cluster == "cluster3") %>%
  select(spot) %>%
  pull()

bc.ranked.TME <- bcUMAP(bc.ranked.TME, 
                        pc = 10, 
                        k.neighbors = 20, 
                        res = 0.3)

cols.highlight <- c("cluster1" = "#4f5be2",
                    "cluster2" = "#ff802b",
                    "non.selected" = "#019453")

bcClusters(bc.ranked.TME, 
           idents = "spot.collapse", 
           spatial = T, mfrow = c(1,2), 
           cells.highlight = list("cluster1" = cluster1,
                                  "cluster2" = cluster2,
                                  "cluster3" = cluster3))

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
spatial.bcScore.TLS <- bcSignatures(bc = bc.TLS.recomputed, UMAP = "Seurat", spatial = T, mfrow = c(2,1), signatures = list(values = "TLS_CABRITA"))

spatial.Bcells <- SpatialFeaturePlot(seurat.TME, features = "B.cells", ncol = 1)
spatial.Tcells <- SpatialFeaturePlot(seurat.TME, features = "T.cells", ncol = 1)

spatial.bcScore.TLS | spatial.Bcells | spatial.Tcells

data.merged <- left_join(data, cluster.selec, by = "spot")
data.boxplot <- ggplot(data.merged, aes(x=cluster, y=TLS_CABRITA)) +
  geom_boxplot()

library(ggstatsplot)

stats.plot <-ggbetweenstats(data = data.merged,
               x = cluster,
               y = TLS_CABRITA)
extract_stats(stats.plot)

data.boxplot | stats.plot
