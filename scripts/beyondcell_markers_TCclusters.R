library("ape")
library("patchwork")

head(bc.recomputed@meta.data)
head(seuratobj.distances@meta.data)
spatial.markers <- FindSpatiallyVariableFeatures(seuratobj.distances, selection.method = "markvariogram")

metadata <- bc.recomputed@meta.data %>%
  select(bc_clusters_res.0.3)

seurat.prueba <- AddMetaData(seuratobj.distances, metadata = metadata)
head(seurat.prueba@meta.data)

Idents(seurat.prueba) <- "bc_clusters_res.0.3"
puretumour.de.markers <- FindMarkers(seurat.prueba, ident.1 = "3", ident.2 = "4")
head(puretumour.de.markers)

ggplot(data = puretumour.de.markers, aes(x=avg_log2FC)) +
  geom_histogram() 
"ERBB2" %in% rownames(puretumour.de.markers)
"ERBB2" %in% rownames(puretumour.filtered)

puretumour.filtered <- puretumour.de.markers %>%
  filter(p_val_adj < 0.05,
         (avg_log2FC < - 1 | avg_log2FC > 1))

puretumour.filtered

SpatialFeaturePlot(object = seurat.prueba, features = "FN1", alpha = c(0.1, 1), ncol = 2)

puretumour.de.markers["ERBB2",]

############################
all.markers <- FindAllMarkers(seurat.prueba)
all.markers.filtered <- all.markers %>%
  filter(p_val_adj < 0.05,
         avg_log2FC < -0.58 | avg_log2FC > 0.58) 



heat.top20 <- DoHeatmap(seurat.prueba, features = top20.genes$gene, group.by = "bc_clusters_res.0.3") +
  theme(axis.text.y = element_text(size = 5))
top20.cluster3 <- top20.genes %>%
  filter(cluster == 3)
SpatialFeaturePlot(object = seurat.prueba, features = top20.cluster3$gene[1:4], alpha = c(0.1, 1), ncol = 4)
SpatialFeaturePlot(object = seurat.prueba, features = top20.cluster3$gene[5:8], alpha = c(0.1, 1), ncol = 4)
SpatialFeaturePlot(object = seurat.prueba, features = top20.cluster3$gene[9:12], alpha = c(0.1, 1), ncol = 4)

top20.cluster4 <- top20.genes %>%
  filter(cluster == 4)
SpatialFeaturePlot(object = seurat.prueba, features = "CRISP3", alpha = c(0.1, 1), ncol = 4)


############

seurat.spatial <- FindSpatiallyVariableFeatures(seurat.prueba, 
                                                assay = "SCT", 
                                                features = VariableFeatures(seurat.prueba)[1:2000],
                                                selection.method = "moransi")
a <- FindSpatiallyVariableFeatures(seurat.prueba, 
                                   assay = "SCT", 
                                   features = VariableFeatures(seurat.prueba)[1:100],
                                   selection.method = "moransi")
a@assays$SCT@meta.features

morans.data <- seurat.spatial@assays$SCT@meta.features
morans.data <- morans.data %>%
  filter(MoransI_p.value < 0.05)


SpatialFeaturePlot(seurat.spatial, features = "DCD", ncol = 2)
SpatialFeaturePlot(seurat.spatial, features = "AC079296.1", ncol = 2)
SpatialFeaturePlot(seurat.spatial, features = "PTPRZ1", ncol = 2)

head(seurat.spatial@assays$SCT@var.features)



top.features <- head(SpatiallyVariableFeatures(seurat.spatial, selection.method = "moransi"), 6)
l <- SpatialFeaturePlot(seurat.spatial, features = top.features, ncol = 3, alpha = c(0.1, 1), combine = F)
class(l)
l[1:2]

spatial.genes <- SpatiallyVariableFeatures(seurat.spatial, selection.method = "moransi")
length(spatial.genes)
top20.genes.cluster4 <- top20.cluster4$gene

top20.genes.cluster4 %in% spatial.genes
l[[1]] | l[[2]]


SpatialFeaturePlot(seurat.spatial, features = spatial.genes[2000], ncol = 3, alpha = c(0.1, 1))

#####

rownames(puretumour.de.markers)
