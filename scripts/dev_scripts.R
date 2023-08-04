rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

## Functions

CellCyclebyIdent <- function(x, split, ident = "Phase") {
  seuratobjs <- SplitObject(x, split.by = split)
  title <- lapply(names(seuratobjs), FUN = function(x) ggtitle(x))
  p <- lapply(seq_along(seuratobjs), FUN = function(i) {
    DimPlot(seuratobjs[[i]], group.by = ident, split.by = "Phase") + title[[i]]
  })
  return(p)
}

## CODE ##
# Random seed
set.seed(1)

# Load RDS
seuratobj.deconvoluted <- readRDS(file = "./results/analysis/seuratobj.deconvoluted.rds")

# Perform SCT normalization without cell cycle regression.
# return.only.var.genes = FALSE: Output the resulting scaled matrix with all genes.
# variable.features.n = NULL; 
# variable.features.rv.th = 1.1: Use the genes with 
# a residual variance > 1.1
## Note: 20230706 - change parameters of SCT to default

seuratobj.deconvoluted <- SCTransform(seuratobj.deconvoluted, 
                                      assay = "Spatial", 
                                      return.only.var.genes = TRUE, 
                                      verbose = TRUE)
seuratobj.deconvoluted <- RunPCA(seuratobj.deconvoluted, 
                                 assay = "SCT", 
                                 npcs = 50, 
                                 features = VariableFeatures(seuratobj.deconvoluted))



### Prueba 

a <- SCTransform(seuratobj.deconvoluted,
                 assay = "Spatial",
                 return.only.var.genes = TRUE, 
                 verbose = TRUE,
                 min_cells = 100)

which(VariableFeatures(a) == "DCD")
which(VariableFeatures(seuratobj.deconvoluted) == "DCD")

length(which(VariableFeatures(a) %in% VariableFeatures(seuratobj.deconvoluted)))

dim.pre.slide.regress <- DimPlot(a) +
  ggtitle("DimPlot without slide regression")

# Regression by slide (orig.ident)
new.slideregress <- SCTransform(a, 
                                      assay = "Spatial", 
                                      vars.to.regress = "orig.ident",
                                      return.only.var.genes = TRUE, 
                                      verbose = TRUE,
                                min_cells = 100)

new.slideregress <- RunPCA(new.slideregress, 
                                 assay = "SCT", 
                                 npcs = 50, features = VariableFeatures(new.slideregress))

dim.with.slide.regress <- DimPlot(new.slideregress) +
  ggtitle("DimPlot slide regressed")

patch.slide.regress <- dim.pre.slide.regress | dim.with.slide.regress


# Cell cycle effect
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
seuratobj.phase <- CellCycleScoring(new.slideregress, 
                                    g2m.features = g2m.genes, 
                                    s.features = s.genes, 
                                    assay = "SCT")
# Run PCA without cell cycle regression
seuratobj.phase <- RunPCA(seuratobj.phase, 
                          assay = "SCT", 
                          npcs = 50, 
                          features = VariableFeatures(seuratobj.phase))

# PCA plot without cell cycle regression
Idents(seuratobj.phase) <- "Phase"
cell.cycle.without <- DimPlot(seuratobj.phase)
cell.cycle.phase.without <- DimPlot(seuratobj.phase, split.by = "Phase")
cell.cycle.phase.without <- cell.cycle.phase.without + 
  geom_vline(xintercept = 80, linetype = "dotted", size = 0.5)
cell.cycle.slide.without <- CellCyclebyIdent(seuratobj.phase, 
                                             split = "orig.ident")

layout <- "
AAABBB
AAACCC
"

patch.cell.cycle.without <- cell.cycle.without + cell.cycle.phase.without + wrap_plots(cell.cycle.slide.without) +
  plot_layout(design = layout) +
  plot_annotation(title = "No cell cycle regression")



head(seuratobj.phase@meta.data)

spots.by.phase <- SpatialDimPlot(seuratobj.phase)

# Distribution of Phase by cell type

df.celltype.phase <- seuratobj.phase@meta.data %>%
  select(c(B.cells:T.cells,Phase)) %>%
  rownames_to_column(var = "spotid") %>%
  pivot_longer(cols = c(B.cells:T.cells),
               names_to = "Cell.type",
               values_to = "Freq") %>%
  group_by(spotid) %>%
  filter(Freq == max(Freq)) 

barplot.celltypes.phases <- ggplot(df.celltype.phase, aes(x = Cell.type)) +
  geom_bar(aes(fill = Phase), position=position_dodge(), colour = "black")

# Cell cycle regression
seuratobj.phase <- SCTransform(seuratobj.phase, 
                               assay = "Spatial", 
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               return.only.var.genes = TRUE, 
                               verbose = TRUE,
                               min_cells = 100)
seuratobj.phase <- RunPCA(seuratobj.phase, 
                          assay = "SCT", 
                          npcs = 50, 
                          features = VariableFeatures(seuratobj.phase))

cell.cycle.with <- DimPlot(seuratobj.phase)
cell.cycle.phase.with <- DimPlot(seuratobj.phase, split.by = "Phase")
cell.cycle.slide.with <- CellCyclebyIdent(seuratobj.phase, 
                                          split = "orig.ident")

patch.cell.cycle.with <- cell.cycle.with + cell.cycle.phase.with + wrap_plots(cell.cycle.slide.with) +
  plot_layout(design = layout) +
  plot_annotation(title = "Cycle regression")
SpatialDimPlot(seuratobj.phase)

# RunPCA made in previous script
# Elbow plot
elbow <- ElbowPlot(seuratobj.phase, 
                   ndims = 50, 
                   reduction = "pca")
elbow
#ggsave(elbow, filename = paste0(out.dir, "/plots/elbow.pdf")) 

# Selection of kNN at different resolution levels

k.param <- c(10,20,30,40,50)
l <- lapply(X = k.param, FUN = function(x){
  res <- c(0.1,0.2,0.3,0.4,0.5)
  a <- FindNeighbors(seuratobj.phase, 
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
       y = "Median of SC3 stability")

# Recompute Clustering at selected K.param
# (in this case k.param = 50)
seuratobj.clusters <- FindNeighbors(seuratobj.phase, 
                                    reduction = "pca", 
                                    dims = 1:20, 
                                    k.param = 50)
res <- c(0.1,0.2,0.3,0.4,0.5)
seuratobj.clusters <- FindClusters(seuratobj.clusters, 
                                   resolution = res, 
                                   verbose = FALSE)

# ReRun UMAP
seuratobj.clusters <- RunUMAP(seuratobj.clusters, 
                              reduction = "pca", 
                              dims = 1:20, 
                              n.components = 2)

# Clustree plots
clustree.plot <- clustree(seuratobj.clusters, 
                          prefix = "SCT_snn_res.",
                          node_colour = "sc3_stability") 

clustree.graph <- clustree(seuratobj.clusters, 
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

clustree.analysis <- (clustree.plot | max.stability.plot)

# Boxplot by cluster, celltypes proportion

cell.types <- names(seuratobj.clusters@meta.data)[6:14]
df.celltypes.clusters <- seuratobj.clusters@meta.data[6:14]
df.celltypes.clusters$cluster <- seuratobj.clusters@meta.data$SCT_snn_res.0.2
df.celltypes.clusters$cluster <- as.factor(df.celltypes.clusters$cluster)

df.celltypes.clusters <- as.data.frame(df.celltypes.clusters)
df.celltypes.clusters.pivot <- pivot_longer(df.celltypes.clusters, 
                                            cols = -cluster, 
                                            names_to = "cell.type", 
                                            values_to = "cell.prop")

boxplot.celltypes.clusters <- ggplot(df.celltypes.clusters.pivot, aes(x = cell.type, y = cell.prop)) +
  geom_boxplot(aes(fill = cell.type)) + 
  facet_wrap(~cluster) +
  ggtitle(label = "Cell types proportions within each cluster") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) 

spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2")
dim.clusters <- DimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1")
boxplot.celltypes.clusters / (spatial.clusters | dim.clusters)

breastcancermarkers.spatial <- SpatialFeaturePlot(seuratobj.clusters, features = c("ESR1", "PGR", "ERBB2"), ncol = 4)
breastcancermarkers.spatial / (spatial.clusters | dim.clusters)

DimPlot(seuratobj.clusters, group.by = "orig.ident")

clustree.plot | dim.clusters

spatial.clusters <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.1")
spatial.clusters2 <- SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2")

spatial.clusters / spatial.clusters2
head(seuratobj.clusters@meta.data)
