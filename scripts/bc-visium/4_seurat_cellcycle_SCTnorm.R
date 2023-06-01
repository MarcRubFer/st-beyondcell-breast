rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)
# Functions

CellCyclebyIdent <- function(x, split, ident = "Phase") {
  seuratobjs <- SplitObject(x, split.by = split)
  title <- lapply(names(seuratobjs), FUN = function(x) ggtitle(x))
  p <- lapply(seq_along(seuratobjs), FUN = function(i) {
    DimPlot(seuratobjs[[i]], group.by = ident, split.by = "Phase") + title[[i]]
  })
  return(p)
}

# Load RDS
seuratobj.deconvoluted <- readRDS(file = "./results/analysis/seuratobj.deconvoluted.rds")

# Perform SCT normalization without cell cycle regression.
# return.only.var.genes = FALSE: Output the resulting scaled matrix with all genes.
# variable.features.n = NULL; 
# variable.features.rv.th = 1.1: Use the genes with 
# a residual variance > 1.1

seuratobj.deconvoluted <- SCTransform(seuratobj.deconvoluted, 
                                  assay = "Spatial", 
                                  return.only.var.genes = FALSE, 
                                  variable.features.n = NULL, 
                                  variable.features.rv.th = 1.1, 
                                  verbose = TRUE)
# Cell cycle effect
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
seuratobj.phase <- CellCycleScoring(seuratobj.deconvoluted, 
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


patch.cell.cycle <- (cell.cycle.without | cell.cycle.phase.without) & 
  theme(legend.position = "bottom")
patch.cell.cycle <- patch.cell.cycle + plot_layout(guides = "collect")
patch.cell.cycle

ggsave(filename = "cell_cycle_withoutregress.png", 
       plot = patch.cell.cycle, 
       path = "./results/plots/")

# Save Data
saveRDS(seuratobj.phase, file = paste0(out.dir,"/analysis/seuratobj.phase.rds"))
all.plots <- list(cell.cycle.without, cell.cycle.phase.without, cell.cycle.phase.without)
save(all.plots, file = paste0("./results/ggplots/cell_cycle.RData"))
