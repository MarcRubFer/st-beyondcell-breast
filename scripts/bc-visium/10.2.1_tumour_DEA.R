rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")

# stablish seed
set.seed(1)

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Read SeuratObjects
seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")

# Analysis of Tumour TCs (TCs 3 to 6)
TCs.tumours <- levels(seuratobj.tcs@meta.data$TCs_res.0.3)[3:6]
cells.TCs.tumours <- seuratobj.tcs@meta.data %>%
  filter(TCs_res.0.3 %in% TCs.tumours) %>%
  rownames_to_column(var = "spots") %>%
  pull(spots)

# Subset seurat object for Tumour TCs cells
seurat.tumour <- subset(seuratobj.tcs, cells = cells.TCs.tumours)

# Restablish levels properly
seurat.tumour$TCs_res.0.3 <- as.character(seurat.tumour$TCs_res.0.3)
seurat.tumour$TCs_res.0.3 <- as.factor(seurat.tumour$TCs_res.0.3)


DefaultAssay(seurat.tumour) <- "Spatial"
Idents(seurat.tumour) <- "TCs_res.0.3"
tumour.dea.gsea <- lapply(TCs.tumours, FUN = function(i) {
  markers <- FindMarkers(seurat.tumour, 
                         ident.1 = "TC-3", 
                         min.pct = 0, 
                         logfc.threshold = 0, 
                         test.use = "wilcox")
  markers <- markers %>%
    rownames_to_column("gene") %>%
    mutate(TC = i)
  return(markers)
}) %>%
  bind_rows()

write.table(tumour.dea.gsea, file = "./results/tumour_TC_DEA.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
