rm(list = ls())

library("beyondcell")
library("Seurat")
library("clustree")
library("tidyverse")
library("tidygraph")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

cell.colors <- c("B-cells" = "cornflowerblue",
                 "CAFs" = "goldenrod2",
                 "Cancer Epithelial" = "red3",
                 "Endothelial" = "seagreen4",
                 "Myeloid" = "darkorchid",
                 "Normal Epithelial" = "hotpink",
                 "Plasmablasts" = "greenyellow",
                 "PVL" = "darkorange1",
                 "T-cells" = "steelblue4")

set.seed(1)

# Read seurat object from each patient
list.seurat.clusters <- readRDS(file = "./results/analysis/reference/list.seurat.reference.clusters.rds")
patients.id <- names(list.seurat.clusters)

# Generate geneset object with one of the ready to use signature collections.
gs.ssc <- GetCollection(SSc)

# Set Assay.
DefaultAssay(list.seurat.clusters) <- "SCT"
# Compute score for the SSc. This might take a few minutes depending on the size of your dataset.
bc <- bcScore(seurat.3921.clusters, gs.ssc, expr.thres = 0.1) 


## Systematic analysis
list.bc <- lapply(list.seurat.clusters, function(object) {
  DefaultAssay(object) <- "SCT"
  bc <- bcScore(object, gs.ssc, expr.thres = 0.1) 
  n.NA <- data.frame(nNAs = colSums(is.na(bc@normalized)),
                     row.names = colnames(bc@normalized))
  bc <- bcAddMetadata(bc, n.NA)
  return(bc)
})
list.bc[[1]]@meta.data


# Barplots distribution nNAs
list.barplot.nas <- lapply(patients.id, function(patient) {
  bc <- list.bc[[patient]]
  df <- bc@meta.data %>%
    select(celltype_major, nNAs) %>%
    mutate(dropped = case_when((nNAs/dim(bc@normalized)[1]) > 0.95 ~ "NaN higher 0.95",
                               TRUE ~ "NaN lower 0.95")) %>%
    count(celltype_major, dropped) 
  
  totals_by_celltype <- df %>%
    group_by(celltype_major) %>%
    summarise(total = sum(n))
  totals_by_celltype
  
  df_percentage <- df %>%
    left_join(totals_by_celltype, by = "celltype_major") %>%
    mutate(percentage = round((n / total) * 100, digits = 0)) %>%
    ggplot(aes(x=factor(celltype_major), y=n, group=dropped)) +
    geom_col(aes(fill=dropped)) +
    geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), size = 3) +
    #ylim(c(0,1500)) +
    ylab("Count") +
    xlab("Cell Type") +
    ggtitle(paste("Number of filtrated cells by NaN - ", patient))
  
  return(df_percentage)
  
})
wrap_plots(list.barplot.nas)

# Distribution of Na in selected cells
list.hist.na <- lapply(patients.id, function(patient){
  bc <- list.bc[[patient]]
  distrib.na <- bc@meta.data %>%
    select(celltype_major, nNAs) %>%
    mutate(dropped = case_when((nNAs/dim(bc@normalized)[1]) > 0.95 ~ "NaN higher 0.95",
                               TRUE ~ "NaN lower 0.95")) %>%
    filter(dropped == "NaN lower 0.95") %>%
    ggplot(aes(x=(nNAs/dim(bc@normalized)[1])*100)) +
    geom_histogram(aes(fill=factor(celltype_major)), ) +
    facet_wrap(vars(celltype_major)) +
    ggtitle(paste("Distribution of Na in selected cells - ", patient))
  
  return(distrib.na)
})
wrap_plots(list.hist.na)


a <- data.frame(nNAs = rowSums(is.na(list.bc[[1]]@normalized)),
                sigs = rowSums(!is.na(list.bc[[1]]@normalized)),
                row.names = rownames(list.bc[[1]]@normalized)) 

a <- a %>%
  mutate(prop.NAs = round((nNAs/(nNAs + sigs)*100), digits = 0),
         prop.sigs = round((sigs/(nNAs + sigs)*100), digits = 0)) %>%
  filter(prop.sigs >= 50) %>%
  arrange(desc(prop.sigs))
a

list.signatures <- lapply(patients.id, function(patient) {
  object <- list.bc[[patient]]
  data <- data.frame(nNAs = rowSums(is.na(object@normalized)),
                  sigs = rowSums(!is.na(object@normalized)))
  data <- data %>%
    mutate(prop.NAs = round((nNAs/(nNAs + sigs)*100), digits = 0),
           prop.sigs = round((sigs/(nNAs + sigs)*100), digits = 0),
           patient_id = patient) %>%
    filter(prop.sigs >= 50) %>%
    select(prop.sigs, patient_id) %>%
    rownames_to_column("sig_id") %>%
    arrange(desc(prop.sigs))
})

list.signatures[[1]]


sig_id_por_df <- lapply(list.signatures, function(df) df$sig_id)
sig_id_comunes <- Reduce(intersect, sig_id_por_df)

df <- list.signatures %>%
  bind_rows() %>%
  filter(sig_id %in% sig_id_comunes)

ggplot(data = df, aes(x=factor(sig_id), y=prop.sigs)) +
  geom_point(aes(col=patient_id)) +
  geom_label(aes(label=patient_id, fill=patient_id), position = position_dodge2(width = 0.45, padding = 0.5)) +
  ylim(c(50,100)) +
  xlab("Cells with BCScore (%)") +
  ylab("BC Signatures") +
  


length(is.na(list.bc[[1]]@normalized))

dim(list.bc[[1]]@normalized)
dim(list.bc[[2]]@normalized)
dim(list.bc[[3]]@normalized)

# Filter out spots with a high percentage of NAs
bc.filtered <- bcSubset(bc, nan.cells = 0.95)
dim(bc@meta.data)
dim(bc.filtered@meta.data)

# Replace NAs by 0s

bc.filtered@normalized[is.na(bc.filtered@normalized)] <- 0
bc.recomputed <- bcRecompute(bc.filtered, slot = "normalized")
dim(bc.recomputed@meta.data)

#Selection of k parameters
k.param <- c(10, 20, 30, 40, 50)
l <- lapply(X = k.param, FUN = function(x){
  res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
  a <- bcUMAP(bc.recomputed, pc =20, k.neighbors = x, res = res)
  clustree.graph <- clustree(a@meta.data, 
                             prefix = "bc_clusters_res.",
                             prop_filter = 0.1, 
                             node_colour = "sc3_stability",
                             return = "graph")
  
  max.stability <- clustree.graph %>%
    activate(nodes) %>%
    as.data.frame() %>%
    group_by(bc_clusters_res.) %>%
    summarise(median.stability = median(sc3_stability),
              n.clusters = length(unique(cluster))) %>%
    mutate(k.param = x)
  return(max.stability)
})

l <- l %>%
  bind_rows() %>%
  mutate(k.param = as.factor(k.param))

clustersc3.kparam <- ggplot(data=l, aes(x=bc_clusters_res., y=median.stability, group = k.param, color = k.param)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability") +
  ggtitle(label = "SC3 stability at different kparam and resolutions")

# Run the bcUMAP function again, specifying the k.params
res <- c(0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
bc.recomputed <- bcUMAP(bc.recomputed, pc = 20, k.neighbors = 30, res = res)
head(bc.recomputed@meta.data)
# Reanalysis of sc3stability for plotting
clustree.plot <- clustree(bc.recomputed@meta.data, 
                          prefix = "bc_clusters_res.",
                          node_colour = "sc3_stability") 

clustree.graph <- clustree(bc.recomputed@meta.data, 
                           prefix = "bc_clusters_res.",
                           prop_filter = 0.1, 
                           node_colour = "sc3_stability",
                           return = "graph")

max.stability <- clustree.graph %>%
  activate(nodes) %>%
  as.data.frame() %>%
  group_by(bc_clusters_res.) %>%
  summarise(median.stability = median(sc3_stability),
            n.clusters = length(unique(cluster)))

max.stability.plot <- ggplot(data=max.stability, aes(x=bc_clusters_res., y=median.stability, group = 1)) + 
  geom_line()+
  geom_point() +
  labs(x = "Resolution",
       y = "Median of SC3 stability")

clustree.plot <- clustree.plot + ggtitle(label = "Clustree with 30 k.param")

patch.clustree <- clustersc3.kparam | clustree.plot

# Visualize whether the cells are clustered based on the number of genes detected per each cell.
bc.nFeatureRNA <- bcClusters(bc.recomputed, 
                             UMAP = "beyondcell", 
                             idents = "nFeature_RNA", 
                             pt.size = 1.5, 
                             factor.col = FALSE) +
  ggtitle("nFeatureRNA")
bc.nCountRNA <- bcClusters(bc.recomputed, 
                           UMAP = "beyondcell", 
                           idents = "nCount_RNA", 
                           pt.size = 1.5, 
                           factor.col = FALSE) +
  ggtitle("nCountsRNA")
bc.phases <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "Phase", pt.size = 1.5)

bc.clusters <- bcClusters(bc.recomputed, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("CID3921 - TC - UMAP Beyondcell (res 0.2)")
bc.clusters.seurat <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "bc_clusters_res.0.2", pt.size = 1) +
  ggtitle("CID3921 - TC - UMAP Seurat (res 0.2)")
bc.seurat <- bcClusters(bc.recomputed, UMAP = "Seurat", idents = "celltype_major", pt.size = 1) +
  ggtitle("CID3921 - Seurat (res 0.2)")
bc.nFeatureRNA | bc.nCountRNA

head(bc.recomputed@meta.data)

layout <- "
##BB
AABB
AACC
##CC
"

patch.bc.clusters <- bc.clusters + bc.clusters.seurat + bc.seurat + 
  plot_layout(design = layout)
plot(patch.bc.clusters)

bc.ranked <- bcRanks(bc.recomputed, idents = "bc_clusters_res.0.2")
list.4squares <- bc4Squares(bc.ranked, idents = "bc_clusters_res.0.2")
class(list.4squares)
(list.4squares[[1]] + list.4squares[[2]]) /
  (list.4squares[[3]] + list.4squares[[4]]) /
  (bc.clusters + bc.clusters.seurat + bc.seurat)

head(bc.ranked@ranks$bc_clusters_res.0.2)

  
wrap_plots(list.4squares)
bcSignatures(bc.ranked, UMAP = "beyondcell", signatures = list(values = "sig-21047"), pt.size = 1.5)
bcHistogram(bc.ranked, signatures = "sig-21047")
bcHistogram(bc.ranked, signatures = "sig-21047", idents = "bc_clusters_res.0.3")


top.diff <- bc.ranked@ranks$bc_clusters_res.0.2 %>%
  select(starts_with("group.")) %>%
  rownames_to_column("signatures") %>%
  pivot_longer(cols = starts_with("group."), names_to = "tcluster", values_to = "type") %>%
  filter(!is.na(type),
         tcluster == "group.2")
head(top.diff, 20)

normalized <- bc.ranked@normalized
dim(normalized)
head(bc.ranked@meta.data)
