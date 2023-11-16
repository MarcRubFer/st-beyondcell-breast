reactome.gmt <- readGMT(x= "./data/gmts/reactome_functional.gmt")

obj <- seuratobj.tcs
obj <- SCTransform(obj, assay = "Spatial", verbose = T, variable.features.n = 10000)

obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,
              rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)
E <- obj@reductions$pca.rev@feature.loadings
set.seed(1)
reactome.gmt

#pctVar -- percent of explained variance along gene set;
gesecaRes <- geseca(reactome.gmt, E, minSize = 15, maxSize = 500, center = FALSE)
head(gesecaRes, 10)
topPathways <- gesecaRes[, pathway] |> head(50)
titles <- sub("HALLMARK_", "", topPathways)
ps <- plotCoregulationProfileSpatial(reactome.gmt[topPathways], 
                                     object = obj,
                                     title=titles)
cowplot::plot_grid(plotlist=ps, ncol=3)


x <- GetAssay(obj, assay=DefaultAssay(obj))
E <- x@scale.data
prefix=""

res <- obj
pathway <- reactome.gmt[[5]]
pathway <- intersect(unique(pathway), rownames(E))

split <- E[pathway, 1:5, drop=FALSE]
colSums(split)
score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
score <- scale(score, center=TRUE, scale=F)

score.p <- scale(score[1], center=TRUE, scale=F)
#res@meta.data[[paste0(prefix, names(reactome.gmt)[2])]] <- score
#head(res@meta.data)
pathways = reactome.gmt[topPathways]
for (i in seq_along(pathways)) {
  pathway <- pathways[[i]]
  pathway <- intersect(unique(pathway), rownames(E))
  
  score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
  score <- scale(score, center=TRUE, scale=T)
  res@meta.data[[paste0(prefix, names(pathways)[i])]] <- score
}

analysis <- res@meta.data %>%
  select(TCs_res.0.3, all_of(starts_with("REACTOME")))

analysis.to.matrix <- analysis %>%
  select(-TCs_res.0.3) 

matrix <- as.matrix.data.frame(analysis.to.matrix)
head(matrix)
tmatrix <- t(matrix)
head(tmatrix)

library(ComplexHeatmap)

col.order <- res@meta.data %>%
  select(TCs_res.0.3, spot.collapse) %>%
  arrange(TCs_res.0.3, spot.collapse)

col.order.spots <- col.order  %>%
  rownames_to_column("spots") %>%
  pull(spots)

tmatrix <- tmatrix[,col.order.spots]
ht <- Heatmap(
  matrix = tmatrix,
  top_annotation = HeatmapAnnotation("TCs" = col.order$TCs_res.0.3,
                                     "Cell type" = col.order$spot.collapse,
                                     col = list("TCs" = TC.colors,
                                                "Cell type" = colors.categories)),
  show_column_names = F,
  cluster_columns = F,
  column_split = col.order$TCs_res.0.3,
  row_names_gp = gpar(fontsize = 6),
  #row_km = 4,
  row_split = 7
)
ht

plotCoregulationProfileSpatial(reactome.gmt["REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"], 
                               object = obj,
                               title=NULL)
a <- SpatialFeaturePlot(res, features = "REACTOME_CELL_CYCLE")

p <- Seurat::SpatialFeaturePlot(res, features = "REACTOME_CELL_CYCLE",
                                combine = FALSE, image.alpha = 0)[[1]]
p$scales$scales[p$scales$find("fill")] <- NULL

# suppress message of replacing existing color palette
colors=c("darkblue", "lightgrey", "darkred")
guide="colourbar"
suppressMessages({
  p2 <- p +
    scale_fill_gradientn(limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3),
                         oob=scales::squish,
                         colors=colors,
                         guide = guide,
                         name = "z-score"
    ) + theme(legend.position = theme_get()$legend.position)
})                   
p2
p

a[[1]] | p | p2
