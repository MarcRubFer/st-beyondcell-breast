rm(list = ls())

#devtools::install_github("igrabski/sc-SHC")
library("scSHC")

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("viridis")
library("RColorBrewer")


out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# stablish seed
set.seed(1)

# Read sSeurat object.
seuratobj.aligned <- readRDS(file = "./results/analysis/seuratobj.aligned.rds")
bc.recomputed <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")

# Export TCs from beyondcell object (for both resolutions)
bc.recomputed@meta.data
metadata.tcs <- bc.recomputed@meta.data %>%
  select(TCs_res.0.3)
tcs <- sort(unique(metadata.tcs$TCs_res.0.3))


# Create new seurat object with TC metadata
seuratobj.tcs <- AddMetaData(seuratobj.aligned, metadata = metadata.tcs)
head(seuratobj.tcs@meta.data)

# Log normalize counts for DEA
DefaultAssay(seuratobj.tcs) <- "Spatial"
seuratobj.tcs <- NormalizeData(seuratobj.tcs)
seuratobj.tcs <-FindVariableFeatures(seuratobj.tcs, 
                                     selection.method = "vst",
                                     nfeatures = 2000)
seuratobj.tcs <- ScaleData(seuratobj.tcs)

# Perform DEA for each TC (for GSEA all genes)
# Resolution 0.1
Idents(seuratobj.tcs) <- "TCs_res.0.3"
dea.gsea <- lapply(tcs, FUN = function(i) {
  markers <- FindMarkers(seuratobj.tcs, 
                         ident.1 = i, 
                         min.pct = 0, 
                         logfc.threshold = 0, 
                         test.use = "wilcox")
  markers <- markers %>%
    rownames_to_column("gene") %>%
    mutate(TC = i)
  return(markers)
}) %>%
  bind_rows()


# Perfom DEA for Overrepresentation analysis (ORA) with logfc 0.5 and 1
dea.ora.0.5 <- FindAllMarkers(seuratobj.tcs, 
                              logfc.threshold = 0.5, 
                              min.pct = 0.1)
dea.ora.0.5 <- dea.ora.0.5 %>%
  relocate(gene) %>%
  rename(TC = cluster)

genes.by.cluster.ora05 <- dea.ora.0.5 %>%
  group_by(TC) %>%
  summarise(n_genes = n())

dea.ora.1 <- FindAllMarkers(seuratobj.tcs, 
                            logfc.threshold = 1, 
                            min.pct = 0.1)
dea.ora.1 <- dea.ora.1 %>%
  relocate(gene) %>%
  rename(TC = cluster)

# Save table
dir.create(path = "./data/dea_tables/")
write.table(dea.gsea, file = "./data/dea_tables/dea.gsea.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(dea.ora.0.5, file = "./data/dea_tables/dea.ora_0.5.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(dea.ora.1, file = "./data/dea_tables/dea.ora_1.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

saveRDS(seuratobj.tcs, file = paste0("./results/analysis/seuratobj.therapeutic.clusters.rds"))
