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
seuratobj.distances <- readRDS("./results/analysis/seuratobj.distances.rds")
bc.ranked.95 <- readRDS("./results/analysis/beyondcell_ranked95.rds")

metadata.tcs <- bc.ranked.95@meta.data %>%
  select(bc_clusters_new_renamed)
tcs <- sort(unique(metadata.tcs$bc_clusters))

seuratobj.tcs <- AddMetaData(seuratobj.distances, metadata = metadata.tcs)
seuratobj.tcs <- subset(seuratobj.tcs, subset = bc_clusters_new_renamed %in% tcs)

# Perform DEA for each TC (for GSEA all genes)
Idents(seuratobj.tcs) <- "bc_clusters_new_renamed"
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

genes.by.cluster <- dea.gsea %>% 
  group_by(TC) %>%
  summarise(n_genes = n())

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

genes.by.cluster.ora1 <- dea.ora.1 %>%
  group_by(TC) %>%
  summarise(n_genes = n())

# Save table
dir.create(path = "./data/dea_tables/")
write.table(dea.gsea, file = "./data/dea_tables/dea.gsea.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(dea.ora.0.5, file = "./data/dea_tables/dea.ora_0.5.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(dea.ora.1, file = "./data/dea_tables/dea.ora_1.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

saveRDS(seuratobj.tcs, file = paste0("./results/analysis/seuratobj.therapeutic.clusters.rds"))
