reactome.gmt <- readGMT(x= "./data/gmts/reactome_functional.gmt")
hallmarks.gmt <- readGMT(x="./data/gmts/hallmarks.all.v2023.1.Hs.symbols.gmt")


obj <- seuratobj.tcs
obj <- SCTransform(obj, assay = "Spatial", new.assay.name = "SCT_GESECA", verbose = T, variable.features.n = 10000)

obj <- RunPCA(obj, assay = "SCT_GESECA", verbose = FALSE,
              rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)
pca.rev <- obj@reductions$pca.rev@feature.loadings
set.seed(1)
reactome.gmt

#pctVar -- percent of explained variance along gene set;
gesecaRes <- geseca(reactome.gmt, pca.rev, minSize = 15, maxSize = 500, center = FALSE)
head(gesecaRes, 10)
topPathways <- gesecaRes[, pathway] |> head(50)

geseca.hallmarks <- geseca(hallmarks.gmt, pca.rev, minSize = 15, maxSize = 500, center = FALSE)
head(geseca.hallmarks, 10)
topPathways.hallmarks <- geseca.hallmarks[, pathway] |> head(50)

x <- GetAssay(obj, assay=DefaultAssay(obj))
E <- x@scale.data
prefix=""

res <- obj
#pathway <- reactome.gmt[[5]]
#pathway <- intersect(unique(pathway), rownames(E))
#
#split <- E[pathway, 1:5, drop=FALSE]
#colSums(split)
#score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
#score <- scale(score, center=TRUE, scale=F)
#
#score.p <- scale(score[1], center=TRUE, scale=F)
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
hallmarks <- hallmarks.gmt[topPathways.hallmarks]
for (i in seq_along(hallmarks)) {
  pathway <- hallmarks[[i]]
  pathway <- intersect(unique(pathway), rownames(E))
  
  score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
  score <- scale(score, center=TRUE, scale=T)
  res@meta.data[[paste0(prefix, names(hallmarks)[i])]] <- score
}
analysis <- res@meta.data %>%
  select(TCs_res.0.3, all_of(starts_with("REACTOME")))

analysis.to.matrix <- analysis %>%
  select(-TCs_res.0.3) 

matrix <- as.matrix.data.frame(analysis.to.matrix)
head(matrix)
tmatrix <- t(matrix)
head(tmatrix)

hallmarks.matrix <- res@meta.data %>%
  select(all_of(starts_with("HALLMARK")))
hallmarks.matrix <- t(as.matrix.data.frame(hallmarks.matrix))

rm(res)
library(ComplexHeatmap)

col.order <- res@meta.data %>%
  select(TCs_res.0.3, spot.collapse) %>%
  arrange(TCs_res.0.3, spot.collapse)

col.order.spots <- col.order  %>%
  rownames_to_column("spots") %>%
  pull(spots)

tmatrix <- tmatrix[,col.order.spots]
hallmarks.matrix <- hallmarks.matrix[ ,col.order.spots]

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

ht.hallmarks <- Heatmap(
  matrix = hallmarks.matrix,
  top_annotation = HeatmapAnnotation("TCs" = col.order$TCs_res.0.3,
                                     "Cell type" = col.order$spot.collapse,
                                     col = list("TCs" = TC.colors,
                                                "Cell type" = colors.categories)),
  show_column_names = F,
  cluster_columns = F,
  column_split = col.order$TCs_res.0.3,
  row_names_gp = gpar(fontsize = 6),
  #row_km = 4,
  row_split = 5
)
ht.hallmarks

plotCoregulationProfileSpatial(reactome.gmt["REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"], 
                               object = obj,
                               title=NULL)
tc.spatial <- SpatialDimPlot(res, group.by = "TCs_res.0.3")
a <- SpatialFeaturePlot(res, features = "REACTOME_COMPLEMENT_CASCADE")

p <- Seurat::SpatialFeaturePlot(res, features = "REACTOME_COMPLEMENT_CASCADE",
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

a[[1]] | p | p2 | tc.spatial[[1]]

tumour.tc2.spots <- res@meta.data %>%
  filter(TCs_res.0.3 == "TC-2" & spot.collapse == "TUMOUR") %>%
  rownames_to_column("spots") %>%
  pull(spots)
