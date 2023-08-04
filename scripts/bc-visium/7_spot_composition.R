rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("ggsci")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

############
## CODE ##
# Random seed
set.seed(1)

# Load seuratobject
seuratobj.aligned <- readRDS(file = "./results/analysis/seuratobj.aligned.rds")
head(seuratobj.aligned@meta.data)

# Subset cell types proportions metadata
cell.type <- seuratobj.aligned@meta.data %>%
  select(c("B.cells","CAFs","Cancer.Epithelial","Endothelial","Myeloid","Normal.Epithelial","Plasmablasts","PVL","T.cells","Cell.Type"))
head(cell.type)

# Reanalysis of first categorization
cell.type <- cell.type %>%
  mutate(Cell.Type = case_when(Cell.Type == "Tumour" ~ "Pure_Tumour",
                                (Cell.Type == "CAFs" & Cancer.Epithelial > 0.4 ~ "Tumour_from_CAF"),
                                (Cell.Type == "Others" & Cancer.Epithelial > 0.4 ~ "Tumour_from_Others"),
                                Cell.Type == "CAFs" ~ "CAFs",
                                Cell.Type == "Others" ~ "Others"))

table <- data.frame(table(factor(cell.type$Cell.Type))) %>%
  arrange(desc(Freq))
table
ggplot(data = cell.type, aes(x=Cell.Type)) +
  geom_bar()

# Collapse proportions of B.cells and T.cells in Lymphoid
cell.type <- cell.type %>%
  mutate(Lymphoid = B.cells + T.cells,
         B.cells = NULL,
         T.cells = NULL) %>%
  relocate(Lymphoid, .before = Cell.Type)
head(cell.type)
# Calculate ranked cell types in each spot (first, second, ...) and its
# relation with the first cell type (prop.celltype.x / prop.celltype1)
filas <- 1:nrow(cell.type)
cell.lines <- names(cell.type)[1:8]
ranked.cell.types <- lapply(filas, function(i) {
  props <- as.double(cell.type[i, cell.lines])
  rank_index <- order(props, decreasing = TRUE)
  main_celltype <- lapply(rank_index, function(index){
    cell.lines[index]
  })
  prop_celltype <- lapply(rank_index, function(index){
    round(props[index]/props[rank_index[1]], digits = 2)
  })
  return(list(main_celltype_1 = main_celltype[[1]], 
              prop_celltype_1 = prop_celltype[[1]],
              main_celltype_2 = main_celltype[[2]], 
              prop_celltype_2 = prop_celltype[[2]],
              main_celltype_3 = main_celltype[[3]], 
              prop_celltype_3 = prop_celltype[[3]],
              main_celltype_4 = main_celltype[[4]], 
              prop_celltype_4 = prop_celltype[[4]],
              main_celltype_5 = main_celltype[[5]], 
              prop_celltype_5 = prop_celltype[[5]],
              main_celltype_6 = main_celltype[[6]], 
              prop_celltype_6 = prop_celltype[[6]],
              main_celltype_7 = main_celltype[[7]], 
              prop_celltype_7 = prop_celltype[[7]],
              main_celltype_8 = main_celltype[[8]], 
              prop_celltype_8 = prop_celltype[[8]]))
})

# Convert list of celltypes and relation in data.frame
ranked_df <- do.call(rbind, lapply(ranked.cell.types, function(x) data.frame(x)))

# Add to general data.frame
cell.type <- cbind(cell.type, ranked_df) 
head(cell.type)

# Subset non-pure tumour cells
df <- cell.type %>%
  filter(Cell.Type != "Pure_Tumour")
head(df)

# Plot distribution of relation prop.celltypex / prop.celltype1
prop2_prop1 <- ggplot(data = df, mapping = aes(x  = prop_celltype_2)) +
  geom_histogram() +
  ggtitle("Celltype 2 / Celltype 1") +
  theme(axis.title.x = element_blank())

prop3_prop1 <- ggplot(data = df, mapping = aes(x  = prop_celltype_3)) +
  geom_histogram(position = position_dodge()) +
  ggtitle("Celltype 3 / Celltype 1") +
  theme(axis.title.x = element_blank())

prop4_prop1 <- ggplot(data = df, mapping = aes(x  = prop_celltype_4)) +
  geom_histogram(position = position_dodge()) +
  ggtitle("Celltype 4 / Celltype 1") +
  theme(axis.title.x = element_blank()) 

patch.props <- prop2_prop1 / prop3_prop1 / prop4_prop1

###############################################################################
# Note: after analyse distributions we stablish thresholds for unique, doublet,
# triplet and mixed spots:
# Unique spots: these have a relation between mayor celltype and second cell type
#               less than 0.5.
# Doublet spots: these have a relation between mayor celltype and second cell type
#               higher than 0.5, and a relation between mayor and third less than 
#               0.5.
# Triplet spots: these have a relation between mayor celltype and thrid cell type
#                higher 0.5
# Mixed spots: the rests of spots.
###############################################################################

# Create spot composition by major class (or sum of classes) and categorical (unique
# doublet, triplet or mixed).

cell.type <- cell.type %>%
  mutate(spot.composition = case_when(Cell.Type == "Pure_Tumour" ~ Cell.Type,
                                      prop_celltype_2 < 0.5 ~ main_celltype_1,
                                      (prop_celltype_2 > 0.5 & prop_celltype_3 < 0.5) ~ paste(main_celltype_1,main_celltype_2, sep = "+"),
                                      prop_celltype_3 > 0.5 ~ paste(main_celltype_1,main_celltype_2,main_celltype_3, sep = "+"),
                                      TRUE ~ "Mixed"),
         spot.composition = case_when(spot.composition == "B.cells" ~ "Lymphoid",
                                      spot.composition == "T.cells" ~ "Lymphoid",
                                      TRUE ~ spot.composition),
         cat.spot.composition = case_when(Cell.Type == "Pure_Tumour" ~ "Unique",
                                          prop_celltype_2 < 0.5 ~ "Unique",
                                          (prop_celltype_2 > 0.5 & prop_celltype_3 < 0.5) ~ "Doublet",
                                          prop_celltype_3 > 0.5 ~ "Triplet",
                                          TRUE ~ "Mixed"))

# Plot by categorical the type of composition and the number of spots in each level

unique_graph <- cell.type %>%
  filter(cat.spot.composition == "Unique") %>%
  ggplot(aes(y=spot.composition)) +
  geom_bar() +
  ggtitle("Composition of unique spots")
unique_graph
doublet_graph <- cell.type %>%
  filter(cat.spot.composition == "Doublet") %>%
  ggplot(aes(y=spot.composition)) +
  geom_bar() +
  ggtitle("Composition of doublet spots")
doublet_graph
triplet_graph <- cell.type %>%
  filter(cat.spot.composition == "Triplet") %>%
  ggplot(aes(y=spot.composition)) +
  geom_bar() +
  ggtitle("Composition of triplet spots")
triplet_graph
mixed_graph <- cell.type %>%
  filter(cat.spot.composition == "Mixed") %>%
  ggplot(aes(y=spot.composition)) +
  geom_bar() +
  ggtitle("Number of spots")
mixed_graph

patch.spot.comp <- (unique_graph | doublet_graph) / (triplet_graph | mixed_graph)

# Note: based on the methodology to create levels, some of them have the same 
# items but in different order. So, we depurate to collapse these levels in the same
# levels.

# Create variable depurated:
#   1) Split each level by symbol "+"
#   2) Sort each element by alphabetic order (function map - purrr package)
#   3) Concate sorted elements with symbol "+" (function map_chr)
cell.type <- cell.type %>%
  mutate(spot.composition.dep = spot.composition %>%
           str_split("\\+") %>%
           map(sort) %>%
           map_chr(paste, collapse = "+"))

# Replot depurated data
doublet_graph2 <- cell.type %>%
  filter(cat.spot.composition == "Doublet") %>%
  ggplot(aes(y=spot.composition.dep)) +
  geom_bar() +
  ggtitle("Composition of doublet spots (dep)")
doublet_graph2
triplet_graph2 <- cell.type %>%
  filter(cat.spot.composition == "Triplet") %>%
  ggplot(aes(y=spot.composition.dep)) +
  geom_bar() +
  ggtitle("Composition of triplet spots (dep)")
triplet_graph2

patch.spot.dep <- (unique_graph | doublet_graph2) / (triplet_graph2 | mixed_graph)

# Note: After analysis of depurated levels we observe that many of them contains
# a few number of spots. So, we decided to filter them stablishing a threshold in 25 spots

# Calculate frequency of each level and get the levels with less than 25 spots
freq.spot.comp.dep <- table(cell.type$spot.composition.dep)
spot.less25 <- names(freq.spot.comp.dep[freq.spot.comp.dep < 20])

# Filter and recategorize these spots to mixed category
cell.type <- cell.type %>%
  mutate(spot.composition.filter = ifelse(test = spot.composition.dep %in% spot.less25, yes = "Mixed", no = spot.composition.dep),
         cat.spot.composition = ifelse(test = spot.composition.filter == "Mixed", yes = "Mixed", no = cat.spot.composition))
head(cell.type)

# Plot filtered results
unique_graph3 <- cell.type %>%
  filter(cat.spot.composition == "Unique") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of unique spots (filtered)")
unique_graph3
doublet_graph3 <- cell.type %>%
  filter(cat.spot.composition == "Doublet") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of doublet spots (filtered)")
doublet_graph3
triplet_graph3 <- cell.type %>%
  filter(cat.spot.composition == "Triplet") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of triplet spots (filtered)")
triplet_graph3
mixed_graph3 <- cell.type %>%
  filter(cat.spot.composition == "Mixed") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Number of mixed spots (filtered)")
mixed_graph3

patch.spot.filtered <- (unique_graph3 | doublet_graph3) / (triplet_graph3 | mixed_graph3)

# Create metadata from cell types 1,2,3 and their relationships, spot comp filtered and categorical spot composition
metadata <- cell.type %>%
  select(main_celltype_1,
         prop_celltype_1,
         main_celltype_2,
         prop_celltype_2,
         main_celltype_3,
         prop_celltype_3,
         spot.composition.filter,
         cat.spot.composition)

# Add metadata to new seuratobject
seuratobj.spotcomp <- AddMetaData(seuratobj.aligned, metadata = metadata)
head(seuratobj.spotcomp@meta.data)


spatial.distrib.spotcomp <- SpatialDimPlot(seuratobj.spotcomp, group.by = "spot.composition.filter", combine = T) 

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/spot_composition/"), recursive = TRUE)
ggsave(filename = "patch_proportions_vs_celltype1.png",
       plot = patch.props,
       path = paste0(out.dir,"/plots/spot_composition/"))
ggsave(filename = "patch_spot_composition_raw.png",
       plot = patch.spot.comp,
       path = paste0(out.dir,"/plots/spot_composition/"))
ggsave(filename = "patch_spot_composition_depurated.png",
       plot = patch.spot.dep,
       path = paste0(out.dir,"/plots/spot_composition/"))
ggsave(filename = "patch_spot_composition_filtered.png",
       plot = patch.spot.filtered,
       path = paste0(out.dir,"/plots/spot_composition/"))
ggsave(filename = "spatial_spot_composition.png",
       plot = spatial.distrib.spotcomp,
       path = paste0(out.dir,"/plots/spot_composition/"))

all.plots <- list(patch.props,patch.spot.comp,patch.spot.dep,patch.spot.filtered,
                  spatial.distrib.spotcomp)
save(all.plots, file = paste0("./results/ggplots/spot_composition.RData"))

# Save data
saveRDS(seuratobj.spotcomp, file = "./results/analysis/seuratobj.spotcomp.rds")


list.spatial.clusters <- SpatialDimPlot(seuratobj.spotcomp, group.by = "SCT_snn_res.0.2", combine = F) 
list.spatial.clusters <- lapply(list.spatial.clusters, function(x){
  x + scale_fill_ucscgb()
})
list.spatial.clusters[[1]]

list.distrib.spotcomp <- SpatialDimPlot(seuratobj.distances, group.by = "spot.composition.filter", combine = F)
list.distrib.spotcomp <- lapply(list.distrib.spotcomp, function(x){
  x + scale_fill_ucscgb()
})
list.distrib.spotcomp[[2]]


