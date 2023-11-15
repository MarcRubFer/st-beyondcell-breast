obj <- seuratobj.tcs
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE, variable.features.n = 10000)

obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,
              rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)
E <- obj@reductions$pca.rev@feature.loadings
set.seed(1)
reactome.gmt

#pctVar -- percent of explained variance along gene set;
gesecaRes <- geseca(reactome.gmt, E, minSize = 15, maxSize = 500, center = FALSE)
head(gesecaRes, 10)
topPathways <- gesecaRes[, pathway] |> head(10)
titles <- sub("HALLMARK_", "", topPathways)
ps <- plotCoregulationProfileSpatial(reactome.gmt[topPathways], 
                                     object = obj,
                                     title=titles)
cowplot::plot_grid(plotlist=ps, ncol=3)


x <- GetAssay(obj, assay=DefaultAssay(obj))
E <- x@scale.data
prefix=""

res <- obj
pathway <- reactome.gmt[[2]]
pathway <- intersect(unique(pathway), rownames(E))

score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
score <- scale(score, center=TRUE, scale=F)
res@meta.data[[paste0(prefix, names(reactome.gmt)[2])]] <- score
head(res@meta.data)

as.data.frame(res@meta.data$REACTOME_ABACAVIR_ADME)
