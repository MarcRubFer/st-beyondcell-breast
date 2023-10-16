rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("ggsci")
library("RColorBrewer")
library("ggplotify")
library("grid")
library("cowplot")

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
  select(B.cells:T.cells)
head(cell.type)

# Collapse proportions of B.cells and T.cells in Lymphoid
cell.type <- cell.type %>%
  mutate(Lymphoid = B.cells + T.cells,
         B.cells = NULL,
         T.cells = NULL) %>%
  relocate(Lymphoid, .after = PVL)
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
#df <- cell.type %>%
#  filter(Cell.Type != "Pure_Tumour")
#head(df)
#
## Plot distribution of relation prop.celltypex / prop.celltype1
#prop2_prop1 <- ggplot(data = df, mapping = aes(x  = prop_celltype_2)) +
#  geom_histogram() +
#  ggtitle("Celltype 2 / Celltype 1") +
#  theme(axis.title.x = element_blank())
#
#prop3_prop1 <- ggplot(data = df, mapping = aes(x  = prop_celltype_3)) +
#  geom_histogram(position = position_dodge()) +
#  ggtitle("Celltype 3 / Celltype 1") +
#  theme(axis.title.x = element_blank())
#
#prop4_prop1 <- ggplot(data = df, mapping = aes(x  = prop_celltype_4)) +
#  geom_histogram(position = position_dodge()) +
#  ggtitle("Celltype 4 / Celltype 1") +
#  theme(axis.title.x = element_blank()) 
#
#patch.props <- prop2_prop1 / prop3_prop1 / prop4_prop1

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
  mutate(spot.composition = case_when(Cancer.Epithelial > 0.65 ~ toupper("Pure Tumour"),
                                      prop_celltype_2 < 0.5 ~ toupper(main_celltype_1),
                                      (prop_celltype_2 > 0.5 & prop_celltype_3 < 0.5) ~ toupper(paste(main_celltype_1,main_celltype_2, sep = "+")),
                                      prop_celltype_3 > 0.5 ~ toupper(paste(main_celltype_1,main_celltype_2,main_celltype_3, sep = "+")),
                                      TRUE ~ "Mixed"),
         spot.composition = gsub(pattern = "\\.", replacement = " ", spot.composition),
         #spot.composition = case_when(spot.composition == "B.cells" ~ toupper("Lymphoid"),
         #                             spot.composition == "T.cells" ~ toupper("Lymphoid"),
         #                             spot.composition == toupper("Cancer.Epithelial") ~ toupper("Cancer Epithelial"),
         #                             TRUE ~ spot.composition),
         cat.spot.composition = case_when(spot.composition == toupper("Pure Tumour") ~ "Unique",
                                          prop_celltype_2 < 0.5 ~ "Unique",
                                          (prop_celltype_2 > 0.5 & prop_celltype_3 < 0.5) ~ "Doublet",
                                          prop_celltype_3 > 0.5 ~ "Triplet",
                                          TRUE ~ "Mixed"))
head(cell.type)
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
# but not are pure tumour (<0.65). For these spots, we decided to re-convert them into "Doublet"

cell.type <- cell.type %>%
  mutate(spot.composition.res = case_when(spot.composition == toupper("Cancer Epithelial") ~ toupper(paste(main_celltype_1,main_celltype_2, sep = "+")),
                                      TRUE ~ spot.composition),
         spot.composition.res = gsub(pattern = "\\.", replacement = " ", spot.composition.res),
         cat.spot.composition.res = case_when(spot.composition == toupper("Cancer Epithelial") ~ "Doublet",
                                          TRUE ~cat.spot.composition))
head(cell.type, n=10)


# Replot depurated data
unique_graph2 <- cell.type %>%
  filter(cat.spot.composition.res == "Unique") %>%
  ggplot(aes(y=spot.composition.res)) +
  geom_bar() +
  ggtitle("Composition of unique spots")
unique_graph2
doublet_graph2 <- cell.type %>%
  filter(cat.spot.composition.res == "Doublet") %>%
  ggplot(aes(y=spot.composition.res)) +
  geom_bar() +
  ggtitle("Composition of doublet spots (dep)")
doublet_graph2
triplet_graph2 <- cell.type %>%
  filter(cat.spot.composition.res == "Triplet") %>%
  ggplot(aes(y=spot.composition.res)) +
  geom_bar() +
  ggtitle("Composition of triplet spots (dep)")
triplet_graph2

patch.spot.dep <- (unique_graph2 | doublet_graph2) / (triplet_graph2 | mixed_graph)

# Note: After analysis of depurated levels we observe that many of them contains
# a few number of spots. So, we decided to filter them stablishing a threshold in 20 spots

# Calculate frequency of each level and get the levels with less than 20 spots
freq.spot.comp.res <- table(cell.type$spot.composition.res)
spot.less20 <- names(freq.spot.comp.res[freq.spot.comp.res < 20])

# Filter and recategorize these spots to doublet category
cell.type <- cell.type %>%
  mutate(spot.composition.filter = ifelse(test = spot.composition.res %in% spot.less20, yes = toupper(paste0(main_celltype_1,"+",main_celltype_2)), no = spot.composition.res),
         spot.composition.filter = gsub(pattern = "\\.", replacement = " ", spot.composition.filter),
         cat.spot.composition.filter = ifelse(test = spot.composition.res %in% spot.less20, yes = "Doublet", no = cat.spot.composition.res))
cell.type <- cell.type %>%
  mutate(spot.composition.filter = spot.composition.filter %>%
           str_split("\\+") %>%
           map(sort) %>%
           map_chr(paste, collapse = "+"))
head(cell.type)

# Plot filtered results
unique_graph3 <- cell.type %>%
  filter(cat.spot.composition.filter == "Unique") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of unique spots (filtered)")
unique_graph3
doublet_graph3 <- cell.type %>%
  filter(cat.spot.composition.filter == "Doublet") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar() +
  ggtitle("Composition of doublet spots (filtered)")
doublet_graph3
triplet_graph3 <- cell.type %>%
  filter(cat.spot.composition.filter == "Triplet") %>%
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

# Collapse categories: first filter of less than 20 spots creates many categories
# So, we decided to collapse in a new categorie ("Other") categories with less than 100 spots
freq.spot.comp.collapse <- table(cell.type$spot.composition.filter)
spot.less100 <- names(freq.spot.comp.collapse[freq.spot.comp.collapse < 100])

## Exclude lymphoid (due to biological interest)
spot.less100.nolymphoid <- spot.less100[which(spot.less100 != toupper("Lymphoid"))] 
cell.type <- cell.type %>%
  mutate(spot.composition.collapse = ifelse(test = spot.composition.filter %in% spot.less100.nolymphoid, yes = "OTHER", no = spot.composition.filter),
         cat.spot.composition.collapse = ifelse(test = spot.composition.filter %in% spot.less100.nolymphoid, yes = "Unique", no = cat.spot.composition.res))

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

cell.type <- cell.type %>%
  mutate(spot.composition.collapse = gsub(pattern = "\\+", replacement = " + ", spot.composition.collapse))

collapsed.graph <- cell.type %>%
  ggplot(aes(y=spot.composition.collapse)) +
  geom_bar()


# Create metadata from cell types 1,2,3 and their relationships and spot comp filtered 
metadata <- cell.type %>%
  select(main_celltype_1,
         prop_celltype_1,
         main_celltype_2,
         prop_celltype_2,
         main_celltype_3,
         prop_celltype_3,
         spot.composition.collapse)

# Add metadata to new seuratobject
seuratobj.spotcomp <- AddMetaData(seuratobj.deconvoluted, metadata = metadata)
head(seuratobj.spotcomp@meta.data)

# Create a named vector for a color palette with categories
# Palette creator: https://medialab.github.io/iwanthue/
colors <- toupper(c("#4fafe3",
                    "#be7adc",
                    "#dec36f",
                    "#7ec967",
                    "#f1703a",
                    "#79696B",
                    "#c4534e"))
levels.spot.comp <- levels(as.factor(seuratobj.spotcomp@meta.data$spot.composition.collapse))
names(colors) <- levels.spot.comp

# Plot spatial distribution of categories
spatial.distrib.spotcomp <- SpatialDimPlot(seuratobj.spotcomp, group.by = "spot.composition.collapse", combine = F) 
spatial.distrib.spotcomp <- lapply(X = seq_along(spatial.distrib.spotcomp), FUN = function(x) {
  spatial.distrib.spotcomp[[x]] +
    scale_fill_manual(name = "Categories",values = colors) +
    theme(legend.key = element_blank())
  
})
spatial.distrib.spotcomp <- wrap_plots(spatial.distrib.spotcomp, ncol = 2) +
  plot_layout(guides = "collect")
#& theme(legend.position = "none")

# Barplot number of spots / categories
collapsed.graph <- cell.type %>%
  ggplot(aes(y=spot.composition.collapse)) +
  geom_bar(aes(fill=spot.composition.collapse)) +
  ylab(NULL) +
  scale_y_discrete(limits = rev, labels = NULL) +
  scale_fill_manual(values = colors)

# Barplot number of spots / categories only for OTHER group (supplementary figure)
other.group.barplot <- cell.type %>%
  filter(spot.composition.collapse == "OTHER") %>%
  ggplot(aes(y=spot.composition.filter)) +
  geom_bar(fill = colors["OTHER"])

(collapsed.graph / other.group.barplot) | spatial.distrib.spotcomp

# Barplot

prop.data <- seuratobj.spotcomp@meta.data %>%
  select(B.cells:T.cells,spot.composition.collapse) %>%
  pivot_longer(cols = B.cells:T.cells, names_to = "cell.type", values_to = "prop.cell.type") %>%
  mutate(cell.type = as.factor(cell.type),
         prop.cell.type = round(prop.cell.type, digits = 2))
boxp.props <- prop.data %>%
  mutate(cell.type = gsub(pattern = "([B,T])\\.", replacement = "\\1-", cell.type),
         cell.type = gsub(pattern = "\\.", replacement = " ", cell.type),
         cell.type = toupper(cell.type)) %>%
  ggplot(aes(x=cell.type, y=prop.cell.type, fill=cell.type)) +
  geom_boxplot() +
  xlab(NULL) +
  ylab("Cell type proportion") +
  ylim(0.0,1.0) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") + 
  facet_grid(~ spot.composition.collapse) +
  theme(strip.text = element_text(size = 6,
                                  face = "bold"))
gt<-ggplot_gtable(ggplot_build(boxp.props))
strips <- which(startsWith(gt$layout$name,'strip'))
for (s in seq_along(strips)) {
  gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- colors[s]
}

boxp.props <- as.ggplot(gt)
collapsed.graph <- collapsed.graph + 
  xlim(c(0,5000)) +
  theme_minimal() + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey"))
upper <- plot_grid(spatial.distrib.spotcomp,NULL,collapsed.graph, ncol = 3, labels = c("A","","B"), rel_widths = c(1,0.01,1))
upper

figure1 <- plot_grid(upper,NULL,boxp.props, ncol = 1, labels = c("","","C"),rel_heights = c(1,0.01,1))

ggsave(filename = "Figure1.png",
       plot = figure1,
       path = paste0(out.dir,"/plots/spot_composition/"))

# Analysis in detail of pure tumor and non-pure tumour
# Subset PURE TUMOUR spots
Idents(seuratobj.spotcomp) <- "spot.composition.collapse"
pure.tumour <- subset(seuratobj.spotcomp, idents = "PURE TUMOUR")

# Subset NON-PURE TUMOUR spots
spot.composition <- levels(as.factor(seuratobj.spotcomp@meta.data$spot.composition.collapse))
non.pure.spots <- spot.composition[spot.composition != "PURE TUMOUR"]
non.pure.tumour <- subset(seuratobj.spotcomp, idents = non.pure.spots)

# Spatial Dimplot of pure and non-pure tumour spots
pure.spatial.dimplot <- SpatialDimPlot(pure.tumour, ncol = 2, combine = F) 
pure.spatial.dimplot <- lapply(X = seq_along(pure.spatial.dimplot), FUN = function(i) {
  pure.spatial.dimplot[[i]] +
    scale_fill_manual(values = colors) +
    theme(legend.position = "none")
})

non.pure.spatial.dimplot <- SpatialDimPlot(non.pure.tumour, ncol = 2, combine = F, image.alpha = 0.7)
non.pure.spatial.dimplot <- lapply(X = seq_along(non.pure.spatial.dimplot), FUN = function(i) {
  non.pure.spatial.dimplot[[i]] +
    scale_fill_manual(values = colors) +
    theme(legend.position = "none")
})

# SpatialFeature plot of cell type proportions
## We excluce Normal, PVL and Plasmablasts for principal figure (Supplementary?)

features <- c("Cancer.Epithelial","CAFs", "Myeloid","B.cells","T.cells","Endothelial")

# Format features for tittles
selected.cell.types <- sapply(seq_along(features), function(x){
  selected.type <- gsub(pattern = "([B,T])\\.", replacement = "\\1-", features[x])
  selected.type <- gsub(pattern = "\\.", replacement = " ", selected.type)
  selected.type <- toupper(selected.type)
  return(selected.type)
})
selected.cell.types <- rep(selected.cell.types, times=2)

#sections <- rep(c("Section 1","Section 2"), length = length(selected.cell.types), each = 6)

# SpatialFeature by images
pure.s1 <- SpatialFeaturePlot(pure.tumour, features = features, combine = F, images = "Section1") 
pure.s2 <- SpatialFeaturePlot(pure.tumour, features = features, combine = F, images = "Section2") 

non.pure.s1 <- SpatialFeaturePlot(non.pure.tumour, features = features, combine = F, images = "Section1", image.alpha = 0.7) 
non.pure.s2 <- SpatialFeaturePlot(non.pure.tumour, features = features, combine = F, images = "Section2", image.alpha = 0.7) 

# Combine in a unique list SpatialDim and SpatialFeature of pure or non-pure
spatial.pure <- append(list(pure.spatial.dimplot[[1]]),pure.s1)
spatial.pure <- append(spatial.pure,list(pure.spatial.dimplot[[2]]))
spatial.pure <- append(spatial.pure,pure.s2)

spatial.non.pure <- append(list(non.pure.spatial.dimplot[[1]]),non.pure.s1)
spatial.non.pure <- append(spatial.non.pure, list(non.pure.spatial.dimplot[[2]]))
spatial.non.pure <- append(spatial.non.pure,non.pure.s2)

## Adjust scales (0 to 1) for all plots and tag each plot
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
spatial.pure.tagged <- lapply(X=seq_along(spatial.pure), FUN = function(i){
  if (i == 1) {
    spatial.pure[[i]]  +
      ggtitle(label = toupper("Pure Tumour spots"), subtitle = "Section 1") +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else if (i %in% c(2:7)) {
    spatial.pure[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = NULL, subtitle = selected.cell.types[i - 1]) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else if (i == 8){
    spatial.pure[[i]]  +
      ggtitle(label = NULL, subtitle = "Section 2") +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else {
    spatial.pure[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = NULL, subtitle = NULL) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
})
pure.patch <- wrap_plots(spatial.pure.tagged, ncol = 7) + plot_layout(guides = "collect")

spatial.non.pure.tagged <- lapply(X=seq_along(spatial.non.pure), FUN = function(i){
  if (i == 1) {
    spatial.non.pure[[i]]  +
      ggtitle(label = toupper("Non-Pure Tumour spots"), subtitle = "Section 1") +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else if (i %in% c(2:7)) {
    spatial.non.pure[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = NULL, subtitle = selected.cell.types[i - 1]) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else if (i == 8){
    spatial.non.pure[[i]]  +
      ggtitle(label = NULL, subtitle = "Section 2") +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else {
    spatial.non.pure[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = NULL, subtitle = NULL) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
})
nonpure.patch <- wrap_plots(spatial.non.pure.tagged, ncol = 7) + plot_layout(guides = "collect") & theme(plot.title  = element_text(size = 12))

#Extract legend for categories
legend.categories <- get_legend(
  (spatial.distrib.spotcomp & theme(legend.position = "bottom")) & guides(fill = guide_legend(ncol = 7))
)


figure2 <- plot_grid(pure.patch, nonpure.patch, legend.categories, ncol = 1, labels = c("A","B",""), rel_heights = c(1,1,0.05), label_size = 18, hjust = -12, vjust = 1)
figure2

ggsave(filename = "Figure2.png",
       plot = figure2,
       path = paste0(out.dir,"/plots/spot_composition/"))

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


