rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE, showWarnings = FALSE)

cell.colors <- c("B-cells" = "cornflowerblue",
                 "CAFs" = "goldenrod2",
                 "Cancer Epithelial" = "red3",
                 "Endothelial" = "seagreen4",
                 "Myeloid" = "darkorchid",
                 "Normal Epithelial" = "hotpink",
                 "Plasmablasts" = "greenyellow",
                 "PVL" = "darkorange1",
                 "T-cells" = "steelblue4")

cell.cycle.colors <- c("G1" = "tomato",
                       "G2M" = "royalblue",
                       "S" = "green3")

## CODE ##
set.seed(1)

# Load seurat reference HER2 patients
seuratobj.refHER2.filtered <- readRDS("./results/analysis/reference/seuratobj.refHER2.filtered.rds")
head(seuratobj.refHER2.filtered)

## Systematic analysis
patients.id <- levels(factor(seuratobj.refHER2.filtered$Patient))

# Create a list of seurat object with each patient ID
list.seurat <- lapply(patients.id, function(x) {
  object <- subset(seuratobj.refHER2.filtered,
                   subset = (orig.ident == x))
  return(object)
})

# Name list entries
names(list.seurat) <- patients.id


# Number of cells by cell type

barplots.list <- lapply(patients.id, function (x) {
  seurat_obj <- list.seurat[[x]]
  celltype_data <- seurat_obj@meta.data$celltype_major
  plot <- data.frame(cell.type = celltype_data) %>%
    count(cell.type) %>%
    ggplot(aes(x=cell.type, y=n)) +
    geom_col(aes(fill = cell.type)) +
    scale_fill_manual(values = cell.colors) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(paste("Cell types in ",x)) 
  return(plot)
})
names(barplots.list) <- patients.id

patchwork.barplots <- wrap_plots(barplots.list) 
ggsave(filename = "patchwork_barplots_patients.png",
       plot = patchwork.barplots,
       path = "./results/plots/reference/beyondcell/")

# SCTransform of seurat objects
## Note: min_cells used as in Spatial analysis

list.seurat.norm <- lapply(list.seurat, FUN = function(x) {
  obj.norm <- SCTransform(x, 
                          assay = "RNA", 
                          return.only.var.genes = TRUE, 
                          verbose = TRUE,
                          min_cells = 100)
  obj.norm <- RunPCA(obj.norm, 
                     assay = "SCT", 
                     npcs = 50, 
                     features = VariableFeatures(obj.norm))
})

# Dimplots seurat objects only normalized

list.dimplot.norm <- lapply(patients.id, FUN = function(x) {
  seurat_obj <- list.seurat.norm[[x]]
  dimplot <- DimPlot(seurat_obj, group.by = "celltype_major", cols = cell.colors) +
    ggtitle(paste("DimPlot normalize patient", x))
})
names(list.dimplot.norm) <- patients.id

patchwork.dimplots <- wrap_plots(list.dimplot.norm)
ggsave(filename = "patchwork_dimplots_patients.png",
       plot = patchwork.dimplots,
       path = "./results/plots/reference/beyondcell/")

list.bar_dim.patients <- lapply(patients.id, function(x) {
  plot <- barplots.list[[x]] | list.dimplot.norm[[x]]
  return(plot)
})
patchwork.bar.dim.patients <- wrap_plots(list.bar_dim.patients, ncol = 2)

ggsave(filename = "patchwork_dim_bar_patients.png",
       plot = patchwork.bar.dim.patients,
       path = "./results/plots/reference/beyondcell/")

# Cell cycle effect
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes

# Run PCA without cell cycle regression
list.seurat.non.phase <- lapply(list.seurat.norm, FUN = function(x) {
  obj.phase <- CellCycleScoring(x,
                                g2m.features = g2m.genes, 
                                s.features = s.genes, 
                                assay = "SCT")
  obj.phase <- RunPCA(obj.phase, 
                      assay = "SCT", 
                      npcs = 50, 
                      features = VariableFeatures(obj.phase))
})


# PCA plot without cell cycle regression
list.cell.cycle.without <- lapply(patients.id, function(x) {
  seurat_obj <- list.seurat.non.phase[[x]]
  Idents(seurat_obj) <- "Phase"
  cell.cycle.without <- DimPlot(seurat_obj, cols = cell.cycle.colors) + ggtitle(x)
  cell.cycle.phase.without <- DimPlot(seurat_obj, split.by = "Phase", cols = cell.cycle.colors)
  patch.cell.cycle.without <- cell.cycle.without + cell.cycle.phase.without 
  return(patch.cell.cycle.without)
})


list.cc.wo.plots <- wrap_plots(list.cell.cycle.without, ncol = 2)  +
  plot_annotation(title = "No cell cycle regression plots")

ggsave(filename = "patchwork_dimplots_cellcycle_without_patients.png",
       plot = list.cc.wo.plots,
       path = "./results/plots/reference/beyondcell/")



# Clustering

# Elbow
list.elbows <- lapply(list.seurat.norm, function(x) {
  elbow <- ElbowPlot(x, 
                     ndims = 50, 
                     reduction = "pca")
  return(elbow)
})
wrap_plots(list.elbows)

# Calculate SC3 kparameters
list.sc3kparam <- lapply(patients.id, function(patient) {
  object <- list.seurat.norm[[patient]]
  k.param <- c(10,20,30,40,50)
  l <- lapply(X = k.param, FUN = function(x){
    res <- c(0.1,0.2,0.3,0.4,0.5)
    a <- FindNeighbors(object, 
                       reduction = "pca", 
                       dims = 1:20, 
                       k.param = x)
    a <- FindClusters(a, 
                      resolution = res, 
                      verbose = FALSE)
    # Run UMAP
    a <- RunUMAP(a, 
                 reduction = "pca", 
                 dims = 1:20, 
                 n.components = 2)
    clustree.graph <- clustree(a, 
                               prefix = "SCT_snn_res.",
                               prop_filter = 0.1, 
                               node_colour = "sc3_stability",
                               return = "graph")
    max.stability <- clustree.graph %>%
      activate(nodes) %>%
      as.data.frame() %>%
      group_by(SCT_snn_res.) %>%
      summarise(median.stability = median(sc3_stability),
                n.clusters = length(unique(cluster))) %>%
      mutate(k.param = x)
    return(max.stability)
  })
  
  l <- l %>%
    bind_rows() %>%
    mutate(k.param = as.factor(k.param))
  
  clustersc3.kparam <- ggplot(data=l, aes(x=SCT_snn_res., y=median.stability, group = k.param, color = k.param)) + 
    geom_line()+
    geom_point() +
    labs(x = "Resolution",
         y = "Median of SC3 stability") +
    ggtitle(paste("SC3 stability in ",patient))
  return(clustersc3.kparam)
})

patchwork.sc3kparam <- wrap_plots(list.sc3kparam)
ggsave(filename = "patchwork_sc3kparam_patients.png",
       plot = patchwork.sc3kparam,
       path = "./results/plots/reference/beyondcell/")

# Recompute Clustering at selected K.param
# (in this case k.param = 20)

list.seurat.clusters <- lapply(list.seurat.norm, function(object) {
  seurat.clusters <- FindNeighbors(object, 
                                        reduction = "pca", 
                                        dims = 1:20, 
                                        k.param = 20)
  res <- c(0.1,0.2,0.3,0.4,0.5)
  seurat.clusters <- FindClusters(seurat.clusters, 
                                       resolution = res, 
                                       verbose = FALSE)
  
  # ReRun UMAP
  seurat.clusters <- RunUMAP(seurat.clusters, 
                                  reduction = "pca", 
                                  dims = 1:20, 
                                  n.components = 2)
  return(seurat.clusters)
})

# Clustree plots
list.clustree.plots <- lapply(patients.id, function(patient) {
  object <- list.seurat.clusters[[patient]]
  clustree.plot <- clustree(object, 
                            prefix = "SCT_snn_res.",
                            node_colour = "sc3_stability") +
    guides(colour = "none", size = "none", count = "none") 
  
  clustree.graph <- clustree(object, 
                             prefix = "SCT_snn_res.",
                             prop_filter = 0.1, 
                             node_colour = "sc3_stability",
                             return = "graph")
  
  max.stability <- clustree.graph %>%
    activate(nodes) %>%
    as.data.frame() %>%
    group_by(SCT_snn_res.) %>%
    summarise(median.stability = median(sc3_stability),
              n.clusters = length(unique(cluster)))
  
  max.stability.plot <- ggplot(data=max.stability, aes(x=SCT_snn_res., y=median.stability, group = 1)) + 
    geom_line()+
    geom_point() +
    labs(x = "Resolution",
         y = "Median of SC3 stability")
  
  clustree.analysis <- clustree.plot +
    ggtitle(patient) 
  return(clustree.analysis)
  
})
patchwork.clustree <- wrap_plots(list.clustree.plots) + plot_layout(guides = "collect")
ggsave(filename = "patchwork_clustree_patients.png",
       plot = patchwork.clustree,
       path = "./results/plots/reference/beyondcell/")

# DimPlot seurat clusters vs cell.types
list.dim.clusters <- lapply(patients.id, function(patient){
  object <- list.seurat.clusters[[patient]]
  dim.clusters <- DimPlot(object, group.by = "SCT_snn_res.0.2") +
    ggtitle(patient)
  dim.celltypes <- DimPlot(object, group.by = "celltype_major", cols = cell.colors) +
    ggtitle(paste(patient,"cell type major"))
  patch <- dim.clusters | dim.celltypes
  return(patch)
})
patchwork.dim.clusters <- wrap_plots(list.dim.clusters, ncol = 2)
ggsave(filename = "patchwork_dim_clusters.png",
       plot = patchwork.dim.clusters,
       path = "./results/plots/reference/beyondcell/")

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/reference/beyondcell/"), recursive = TRUE)

ggsave(filename = "CID3921.barplot.celltypes.png",
       plot = barplot.celltypes,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))
ggsave(filename = "CID3921.dimplot.normalize.png",
       plot = dimplot.3921.norm,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))
ggsave(filename = "CID3921.cellcycle.without.png",
       plot = patch.cell.cycle.without,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))
ggsave(filename = "CID3921.cellcycle.with.png",
       plot = patch.cell.cycle.with,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))
ggsave(filename = "CID3921.dimplot.phase.png",
       plot = dimplot.3921.phase,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))
ggsave(filename = "CID3921.clustersc3.kparam.png",
       plot = CID3921.clustersc3.kparam,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))
ggsave(filename = "CID3921.clustree.analysis.png",
       plot = clustree.analysis,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))
ggsave(filename = "CID3921.dimplot.clusters_vs_celltype.png",
       plot = dim.clusters,
       path = paste0(out.dir,"/plots/reference/beyondcell/"))

# Save data
saveRDS(list.seurat.clusters, file = "./results/analysis/reference/list.seurat.reference.clusters.rds")
