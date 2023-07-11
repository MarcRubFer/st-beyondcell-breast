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
## Update 2023-07-10: add min_cells = 100 (min number of cells that express a gene used for SCTprocessing, default 5)

seuratobj.deconvoluted <- SCTransform(seuratobj.deconvoluted, 
                                  assay = "Spatial", 
                                  return.only.var.genes = TRUE, 
                                  verbose = TRUE,
                                  min_cells = 100)
seuratobj.deconvoluted <- RunPCA(seuratobj.deconvoluted, 
                                 assay = "SCT", 
                                 npcs = 50, 
                                 features = VariableFeatures(seuratobj.deconvoluted))

dim.pre.slide.regress <- DimPlot(seuratobj.deconvoluted) +
  ggtitle("DimPlot without slide regression")

# Regression by slide (orig.ident)
seuratobj.slideregress <- SCTransform(seuratobj.deconvoluted, 
                 assay = "Spatial", 
                 vars.to.regress = "orig.ident",
                 return.only.var.genes = TRUE, 
                 verbose = TRUE,
                 min_cells = 100)

seuratobj.slideregress <- RunPCA(seuratobj.slideregress, 
                                 assay = "SCT", 
                                 npcs = 50, features = VariableFeatures(seuratobj.slideregress))

dim.with.slide.regress <- DimPlot(seuratobj.slideregress) +
  ggtitle("DimPlot slide regressed")

patch.slide.regress <- dim.pre.slide.regress | dim.with.slide.regress


# Cell cycle effect
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
seuratobj.phase <- CellCycleScoring(seuratobj.slideregress, 
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

## Alternative regression
## (see alternative workflow here : https://satijalab.org/seurat/articles/cell_cycle_vignette.html#:~:text=This%20means%20that%20signals%20separating,regressed%20out%20of%20the%20data)
## C
seuratobj.phase.alt <- seuratobj.phase
seuratobj.phase.alt$CC.difference <- seuratobj.phase.alt$S.Score - seuratobj.phase.alt$G2M.Score
seuratobj.phase.alt <- SCTransform(seuratobj.phase.alt, 
                 assay = "Spatial", 
                 vars.to.regress = "CC.difference",
                 return.only.var.genes = TRUE, 
                 verbose = TRUE,
                 min_cells = 100)
seuratobj.phase.alt <- RunPCA(seuratobj.phase.alt, 
            assay = "SCT", 
            npcs = 50, 
            features = VariableFeatures(seuratobj.phase.alt))

cell.cycle.with.alt <- DimPlot(seuratobj.phase.alt)
cell.cycle.phase.with.alt <- DimPlot(seuratobj.phase.alt, split.by = "Phase")
cell.cycle.slide.with.alt <- CellCyclebyIdent(seuratobj.phase.alt, 
                                          split = "orig.ident")

patch.cell.cycle.with.alt <- cell.cycle.with.alt + cell.cycle.phase.with.alt + wrap_plots(cell.cycle.slide.with.alt) +
  plot_layout(design = layout) +
  plot_annotation(title = "Cycle regression alternative (difference between the G2M and S)")
SpatialDimPlot(seuratobj.phase.alt)

# Save plots and ggplots

dir.create(path = paste0(out.dir,"/plots/SCT_cellcycle"), recursive = TRUE)
ggsave(filename = "patch_slide_regress.png", 
       plot = patch.slide.regress, 
       path = paste0(out.dir,"/plots/SCT_cellcycle"))
ggsave(filename = "cell_cycle_withoutregress.png", 
       plot = patch.cell.cycle.without, 
       path = "./results/plots/SCT_cellcycle/")
ggsave(filename = "spots_by_phases.png", 
       plot = spots.by.phase, 
       path = "./results/plots/SCT_cellcycle/")
ggsave(filename = "barplot_celltypes_phases.png",
       plot = barplot.celltypes.phases,
       path = "./results/plots/SCT_cellcycle/")
ggsave(filename = "cell_cycle_regressed.png", 
       plot = patch.cell.cycle.with, 
       path = "./results/plots/SCT_cellcycle/")
ggsave(filename = "cell_cycle_regressed_alternative.png", 
       plot = patch.cell.cycle.with.alt, 
       path = "./results/plots/SCT_cellcycle/")

all.plots <- list(cell.cycle.without, cell.cycle.phase.without, cell.cycle.slide.without, patch.cell.cycle.without,
                  spots.by.phase,barplot.celltypes.phases,
                  cell.cycle.with, cell.cycle.phase.with, cell.cycle.slide.with, patch.cell.cycle.with,
                  cell.cycle.with.alt, cell.cycle.phase.with.alt, cell.cycle.slide.with.alt, patch.cell.cycle.with.alt)
save(all.plots, file = paste0("./results/ggplots/cell_cycle.RData"))

# Save Data
saveRDS(seuratobj.phase, file = paste0(out.dir,"/analysis/seuratobj.phase.rds"))
saveRDS(seuratobj.phase.alt, file = paste0(out.dir,"/analysis/seuratobj.phase_alternative.rds"))
