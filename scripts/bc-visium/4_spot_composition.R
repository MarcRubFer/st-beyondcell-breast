rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("ggsci")
library("RColorBrewer")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

############
## CODE ##
# Random seed
set.seed(1)

# Load seuratobject
seuratobj.deconvoluted <- readRDS(file = "./results/analysis/seuratobj.deconvoluted.rds")
head(seuratobj.deconvoluted@meta.data)

# Subset cell types proportions metadata
cell.type <- seuratobj.deconvoluted@meta.data %>%
  select(c("B.cells","CAFs","Cancer.Epithelial","Endothelial","Myeloid","Normal.Epithelial","Plasmablasts","PVL","T.cells","Cell.Type"))
head(cell.type)

# Reanalysis of first categorization
#cell.type <- cell.type %>%
#  mutate(Cell.Type = case_when(Cell.Type == "Tumour" ~ "Pure_Tumour",
#                                (Cell.Type == "CAFs" & Cancer.Epithelial > 0.4 ~ "Tumour_from_CAF"),
#                                (Cell.Type == "Others" & Cancer.Epithelial > 0.4 ~ "Tumour_from_Others"),
#                                Cell.Type == "CAFs" ~ "CAFs",
#                                Cell.Type == "Others" ~ "Others"))
#
#table <- data.frame(table(factor(cell.type$Cell.Type))) %>%
#  arrange(desc(Freq))
#table
#ggplot(data = cell.type, aes(x=Cell.Type)) +
#  geom_bar()

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
  mutate(spot.composition = case_when(Cancer.Epithelial > 0.65 ~ "Pure_Tumour",
                                      prop_celltype_2 < 0.5 ~ main_celltype_1,
                                      (prop_celltype_2 > 0.5 & prop_celltype_3 < 0.5) ~ paste(main_celltype_1,main_celltype_2, sep = "+"),
                                      prop_celltype_3 > 0.5 ~ paste(main_celltype_1,main_celltype_2,main_celltype_3, sep = "+"),
                                      TRUE ~ "Mixed"),
         spot.composition = case_when(spot.composition == "B.cells" ~ "Lymphoid",
                                      spot.composition == "T.cells" ~ "Lymphoid",
                                      TRUE ~ spot.composition),
         cat.spot.composition = case_when(spot.composition == "Pure_Tumour" ~ "Unique",
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

# Note:
# 1) Due to threshold of "pure tumour" some spots are unique for "cancer epithelial"
# but not are pure tumour. For this spots, we decided to re-convert them into "Doublet"

cell.type <- cell.type %>%
  mutate(spot.composition.res = case_when(spot.composition == "Cancer.Epithelial" ~ paste0(main_celltype_1,"+",main_celltype_2),
                                      TRUE ~ spot.composition),
         cat.spot.composition.res = case_when(spot.composition == "Cancer.Epithelial" ~ "Doublet",
                                          TRUE ~cat.spot.composition))

# 2)based on the methodology to create levels, some of them have the same 
# items but in different order. So, we depurate to collapse these levels in the same
# levels.

# Create variable depurated:
#   1) Split each level by symbol "+"
#   2) Sort each element by alphabetic order (function map - purrr package)
#   3) Concate sorted elements with symbol "+" (function map_chr)
cell.type <- cell.type %>%
  mutate(spot.composition.dep = spot.composition.res %>%
           str_split("\\+") %>%
           map(sort) %>%
           map_chr(paste, collapse = "+"))

# Replot depurated data
unique_graph2 <- cell.type %>%
  filter(cat.spot.composition.res == "Unique") %>%
  ggplot(aes(y=spot.composition)) +
  geom_bar() +
  ggtitle("Composition of unique spots")
unique_graph2
doublet_graph2 <- cell.type %>%
  filter(cat.spot.composition.res == "Doublet") %>%
  ggplot(aes(y=spot.composition.dep)) +
  geom_bar() +
  ggtitle("Composition of doublet spots (dep)")
doublet_graph2
triplet_graph2 <- cell.type %>%
  filter(cat.spot.composition.res == "Triplet") %>%
  ggplot(aes(y=spot.composition.dep)) +
  geom_bar() +
  ggtitle("Composition of triplet spots (dep)")
triplet_graph2

patch.spot.dep <- (unique_graph2 | doublet_graph2) / (triplet_graph2 | mixed_graph)

# Note: After analysis of depurated levels we observe that many of them contains
# a few number of spots. So, we decided to filter them stablishing a threshold in 20 spots

# Calculate frequency of each level and get the levels with less than 20 spots
freq.spot.comp.dep <- table(cell.type$spot.composition.dep)
spot.less20 <- names(freq.spot.comp.dep[freq.spot.comp.dep < 20])

# Filter and recategorize these spots to doublet category
cell.type <- cell.type %>%
  mutate(spot.composition.filter = ifelse(test = spot.composition.dep %in% spot.less20, yes = paste0(main_celltype_1,"+",main_celltype_2), no = spot.composition.dep),
         cat.spot.composition.res = ifelse(test = spot.composition.dep %in% spot.less20, yes = "Doublet", no = cat.spot.composition.res))
cell.type <- cell.type %>%
  mutate(spot.composition.filter = spot.composition.filter %>%
           str_split("\\+") %>%
           map(sort) %>%
           map_chr(paste, collapse = "+"))
head(cell.type)

# Plot filtered results
unique_graph3 <- cell.type %>%
  filter(cat.spot.composition.res == "Unique") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of unique spots (filtered)")
unique_graph3
doublet_graph3 <- cell.type %>%
  filter(cat.spot.composition.res == "Doublet") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of doublet spots (filtered)")
doublet_graph3
triplet_graph3 <- cell.type %>%
  filter(cat.spot.composition.res == "Triplet") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of triplet spots (filtered)")
triplet_graph3
mixed_graph3 <- cell.type %>%
  filter(cat.spot.composition.res == "Mixed") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Number of mixed spots (filtered)")
mixed_graph3

patch.spot.filtered <- (unique_graph3 | doublet_graph3) / (triplet_graph3 | mixed_graph3)

head(cell.type)

# Collapse categories
freq.spot.comp.collapse <- table(cell.type$spot.composition.filter)
spot.less100 <- names(freq.spot.comp.collapse[freq.spot.comp.collapse < 100])
## Exclude lymphoid (due to biological interest)
cell.type <- cell.type %>%
  mutate(spot.composition.collapse = ifelse(test = spot.composition.filter %in% spot.less100, yes = "Other", no = spot.composition.filter),
         cat.spot.composition.collapse = ifelse(test = spot.composition.filter %in% spot.less100, yes = "Unique", no = cat.spot.composition.res))

unique_graph4 <- cell.type %>%
  filter(cat.spot.composition.collapse == "Unique") %>%
  ggplot(aes(y=spot.composition.collapse)) +
  geom_bar() +
  ggtitle("Composition of unique spots (collapsed)")
unique_graph4
doublet_graph4 <- cell.type %>%
  filter(cat.spot.composition.collapse == "Doublet") %>%
  ggplot(aes(y=spot.composition.collapse)) +
  geom_bar() +
  ggtitle("Composition of doublet spots (collapsed)")
doublet_graph4
triplet_graph4 <- cell.type %>%
  filter(cat.spot.composition.collapse == "Triplet") %>%
  ggplot(aes(y=spot.composition.collapse)) +
  geom_bar() +
  ggtitle("Composition of triplet spots (collapsed)")
triplet_graph4
mixed_graph3 <- cell.type %>%
  filter(cat.spot.composition.res == "Mixed") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Number of mixed spots (filtered)")
mixed_graph3

patch.spot.collapsed <- (unique_graph4 | doublet_graph4) / (triplet_graph4 | plot_spacer())
patch.spot.collapsed

collapsed.graph <- cell.type %>%
  ggplot(aes(y=spot.composition.collapse)) +
  geom_bar()

prueba <- cell.type %>%
  filter(spot.composition.collapse == "Other") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar()
# Create metadata from cell types 1,2,3 and their relationships, spot comp filtered and categorical spot composition
metadata <- cell.type %>%
  select(main_celltype_1,
         prop_celltype_1,
         main_celltype_2,
         prop_celltype_2,
         main_celltype_3,
         prop_celltype_3,
         spot.composition.collapse,
         cat.spot.composition.collapse)

# Add metadata to new seuratobject
seuratobj.spotcomp <- AddMetaData(seuratobj.deconvoluted, metadata = metadata)
head(seuratobj.spotcomp@meta.data)

colors <- toupper(c("#4fafe3",
                    "#4e3f41",
                    "#bca251",
                    "#7ec967",
                    "#9451b2",
                    "#c4534e"))
levels.spot.comp <- levels(as.factor(seuratobj.spotcomp@meta.data$spot.composition.collapse))
names(colors) <- levels.spot.comp

spatial.distrib.spotcomp <- SpatialDimPlot(seuratobj.spotcomp, group.by = "spot.composition.collapse", combine = F) 
spatial.distrib.spotcomp <- lapply(X = seq_along(spatial.distrib.spotcomp), FUN = function(x) {
  spatial.distrib.spotcomp[[x]] +
    scale_fill_manual(values = colors)
  
})
spatial.distrib.spotcomp <- wrap_plots(spatial.distrib.spotcomp)
spatial.distrib.spotcomp

# Subset Pure_tumour spots
Idents(seuratobj.spotcomp) <- "spot.composition.filter"
pure.tumour <- subset(seuratobj.spotcomp, idents = "Pure_Tumour")

pure.cell.prop <- pure.tumour@meta.data %>%
  select(B.cells:T.cells) %>%
  pivot_longer(cols = everything(), names_to = "cell.type", values_to = "prop.cell.type") %>%
  mutate(cell.type = as.factor(cell.type),
         prop.cell.type = round(prop.cell.type, digits = 2))

p <- SpatialDimPlot(pure.tumour, ncol = 2, combine = F) 
p <- lapply(X = seq_along(p), FUN = function(i) {
  p[[i]] +
    scale_fill_manual(values = colors)
})
p <- wrap_plots(p)
p <- p + plot_layout(guides = "collect") &
  theme(legend.position = "top")
p[[2]] <- p[[2]] +
  theme(legend.position = "none")

q <- pure.cell.prop %>%
  filter(prop.cell.type > 0.00) %>%
  ggplot(aes(x=cell.type, y=prop.cell.type, fill=cell.type)) +
  #geom_violin(scale = "count", trim = F) #+
  #geom_jitter(alpha =0.2) +
  geom_boxplot()

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
r <- SpatialFeaturePlot(pure.tumour, features = c("Cancer.Epithelial","CAFs", "Myeloid"), ncol = 2, combine = F) 
#r <- SpatialFeaturePlot(pure.tumour, features = c("CAFs"), ncol = 2, combine = F) 
r <- lapply(X=seq_along(r), FUN = function(i){
  r[[i]]  +
    scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
    ggtitle(c("Hello","world")) +
    theme(legend.title = element_blank(),
          legend.direction = "vertical") 
})
r <- wrap_plots(r, ncol = 2, nrow = 3)
r + plot_layout(guides = "collect")
patch <- (p / q) | r
patch

# Subset non-pure tumour spots
spot.composition <- levels(as.factor(seuratobj.spotcomp@meta.data$spot.composition.filter))
non.pure.spots <- spot.composition[spot.composition != "Pure_Tumour"]
non.pure.tumour <- subset(seuratobj.spotcomp, idents = non.pure.spots)
non.pure.cell.prop <- non.pure.tumour@meta.data %>%
  select(B.cells:T.cells) %>%
  pivot_longer(cols = everything(), names_to = "cell.type", values_to = "prop.cell.type") %>%
  mutate(cell.type = as.factor(cell.type),
         prop.cell.type = round(prop.cell.type, digits = 2))

p2 <- SpatialDimPlot(non.pure.tumour, ncol = 2)

p2 <- p2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p2[[2]] <- p2[[2]] +
  theme(legend.position = "none")
p2
q2 <- ggplot(non.pure.cell.prop, aes(x=cell.type, y=prop.cell.type, fill=cell.type)) +
  geom_jitter(alpha =0.2) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0.0,1.0)
r2 <- SpatialFeaturePlot(non.pure.tumour, features = c("Cancer.Epithelial","CAFs", "Myeloid","B.cells","T.cells","Endothelial"), ncol = 4) 

r2.1 <- SpatialFeaturePlot(non.pure.tumour, features = c("B.cells","T.cells","Endothelial"), ncol = 4) 
patch2 <- (p2 / q2) | (r2)
patch2

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


