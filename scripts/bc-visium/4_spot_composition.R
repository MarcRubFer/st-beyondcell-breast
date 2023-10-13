rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("ggsci")
library("RColorBrewer")
library("ggplotify")

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
& theme(legend.position = "none")

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

plot_grid(upper,NULL,boxp.props, ncol = 1, labels = c("","","C"),rel_heights = c(1,0.01,1))
(spatial.distrib.spotcomp | collapsed.graph) / boxp.props
(wrap_elements(full = spatial.distrib.spotcomp) | collapsed.graph) / boxp.props

# Analysis in detail of pure tumor and non-pure tumour
# Subset PURE TUMOUR spots
Idents(seuratobj.spotcomp) <- "spot.composition.collapse"
pure.tumour <- subset(seuratobj.spotcomp, idents = "PURE TUMOUR")

# Data frame of cell type and proportion
pure.cell.prop <- pure.tumour@meta.data %>%
  select(B.cells:T.cells,spot.composition.collapse) %>%
  pivot_longer(cols = B.cells:T.cells, names_to = "cell.type", values_to = "prop.cell.type") %>%
  mutate(cell.type = as.factor(cell.type),
         prop.cell.type = round(prop.cell.type, digits = 2))

# Spatial plot of pure tumour spots
p <- SpatialDimPlot(pure.tumour, ncol = 2, combine = F) 
p <- lapply(X = seq_along(p), FUN = function(i) {
  p[[i]] +
    scale_fill_manual(values = colors)
})
p <- wrap_plots(p)
p <- p + plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.key.size = unit(.5, 'cm'))

p[[2]] <- p[[2]] +
  theme(legend.position = "none")
p
# Boxplot of cell type proportions in pure tumour spot
q <- pure.cell.prop %>%
  mutate(cell.type = gsub(pattern = "([B,T])\\.", replacement = "\\1-", cell.type),
         cell.type = gsub(pattern = "\\.", replacement = " ", cell.type),
         cell.type = toupper(cell.type)) %>%
  ggplot(aes(x=cell.type, y=prop.cell.type, fill=cell.type)) +
  geom_boxplot() +
  facet_grid(~spot.composition.collapse) +
  theme_minimal() +
  ylab(label = "Cell type proportion") +
  labs(fill = "Cell type") +
  theme(strip.background = element_rect(fill = colors["PURE TUMOUR"]),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1)
        ) 
q
# Spatial plot of three main cell type proportions
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
r <- SpatialFeaturePlot(pure.tumour, features = c("Cancer.Epithelial","CAFs", "Myeloid"), ncol = 2, combine = F) 

selected.cell.types <- toupper(rep(c("Cancer Epithelial","CAFs", "Myeloid"), each=2))
sections <- rep(c("Section 1","Section 2"), length = length(selected.cell.types))
## Adjust scales (0 to 1) for all plots and tag each plot
r <- lapply(X=seq_along(r), FUN = function(i){
  if (i <=2) {
    r[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = sections[i], subtitle = selected.cell.types[i]) +
      ylab(label = "Hello") +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }else{
    r[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = NULL, subtitle = selected.cell.types[i]) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
   
})
r <- wrap_plots(r, ncol = 6, nrow = 1)
r <- r + plot_layout(guides = "collect")
r

patch.pure.tumour <- (p | q) / r
patch.pure.tumour

((collapsed.graph | spatial.distrib.spotcomp)  / ((p[[1]] / p[[2]]) | q)) | r

# Subset non-pure tumour spots
spot.composition <- levels(as.factor(seuratobj.spotcomp@meta.data$spot.composition.collapse))
non.pure.spots <- spot.composition[spot.composition != "PURE TUMOUR"]
non.pure.tumour <- subset(seuratobj.spotcomp, idents = non.pure.spots)
non.pure.cell.prop <- non.pure.tumour@meta.data %>%
  select(B.cells:T.cells,spot.composition.collapse) %>%
  pivot_longer(cols = B.cells:T.cells, names_to = "cell.type", values_to = "prop.cell.type") %>%
  mutate(spot.composition.collapse = as.factor(spot.composition.collapse),
         cell.type = as.factor(cell.type),
         prop.cell.type = round(prop.cell.type, digits = 2))

p2 <- SpatialDimPlot(non.pure.tumour, ncol = 2, combine = F, image.alpha = 0.7)
p2 <- lapply(X = seq_along(p2), FUN = function(i) {
  p2[[i]] +
    scale_fill_manual(values = colors)
})
p2 <- wrap_plots(p2)
p2 <- p2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p2

q2 <- non.pure.cell.prop %>%
  mutate(cell.type = gsub(pattern = "([B,T])\\.", replacement = "\\1-", cell.type),
         cell.type = gsub(pattern = "\\.", replacement = " ", cell.type),
         cell.type = toupper(cell.type)) %>%
  ggplot(aes(x=cell.type, y=prop.cell.type, fill=cell.type)) +
  geom_boxplot() +
  xlab(NULL) +
  ylim(0.0,1.0) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  facet_grid(~ spot.composition.collapse) +
  theme(strip.text = element_text(size = 6,
                                  face = "bold"))
q2

gt<-ggplot_gtable(ggplot_build(q2))
strips <- which(startsWith(gt$layout$name,'strip'))
for (s in seq_along(strips)) {
  gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- colors[s]
}
q2 <- plot(gt)
q2 <- as.ggplot(gt)
q2

r2 <- SpatialFeaturePlot(non.pure.tumour, features = c("Cancer.Epithelial","CAFs", "Myeloid","B.cells","T.cells","Endothelial"), combine = F) 
selected.non.pure <- rep(c("Cancer Epithelial","CAFs", "Myeloid","B cells","T cells","Endothelial"),each = 2)
sections.non.pure <- rep(c("Section 1","Section 2"), length = length(selected.non.pure))
colors.text <- rep(c(scales::hue_pal()(9)), each = 2)
r2 <- lapply(X=seq_along(r2), FUN = function(i){
  if (i %% 2 != 0) {
    r2[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = toupper(selected.non.pure[i])) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }else{
    r2[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = NULL, subtitle = selected.non.pure[i]) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
  
})

r2 <- lapply(X=seq_along(r2), FUN = function(i){
  r2[[i]]  +
      scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
      ggtitle(label = toupper(selected.non.pure[i]), subtitle = sections.non.pure[i]) +
      theme(legend.title = element_blank(),
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5,
                                         size = 6),
            legend.position = "none")
  })
r2 <- wrap_plots(r2, ncol = 6, nrow = 2)
r2 <- r2 + plot_layout(guides = "collect")

patch2 <- ((p2 | q2) / (r2)) 

a <- (((wrap_elements(panel = textGrob("Section 1")) / wrap_elements(panel = (textGrob("Section 2")))) | ((r2.s1 / r2.s1) + plot_layout(guides = "collect"))) + plot_layout(guides = "collect", widths = c(1,10)))
a
patch2 <- (((p2 | q2) + plot_layout(guides = "collect")) / (a)) 
patch2

ggsave(filename = "a.png",
       plot = a,
       path = paste0(out.dir,"/plots/spot_composition/"))
r2.s1 <- SpatialFeaturePlot(non.pure.tumour, features = c("Cancer.Epithelial","CAFs", "Myeloid","B.cells","T.cells","Endothelial"), images = "Section1",combine = F) 
selected.non.pure <- c("Cancer Epithelial","CAFs", "Myeloid","B cells","T cells","Endothelial")
r2.s1 <- lapply(X=seq_along(r2.s1), FUN = function(i){
  r2.s1[[i]]  +
    scale_fill_gradientn(limits = c(0.0,1), colours=SpatialColors(n=100)) +
    ggtitle(label = toupper(selected.non.pure[i])) +
    theme(legend.title = element_blank(),
          legend.direction = "vertical",
          plot.title = element_text(hjust = 0.5))
})
r2.s1 <- wrap_plots(r2.s1, ncol = 6) + plot_layout(guides = "collect")
library(grid)
(((wrap_elements(panel = textGrob("here\nhere")) | r2.s1) + plot_layout(widths = c(1,20))) & theme(plot.background = element_rect(fill = "blue"))) / ((wrap_elements(panel = textGrob("Section 2")) | r2.s1) + plot_layout(widths = c(1,20)))

wrap_elements(panel = (textGrob("here") + r2.s1))

((wrap_elements(panel = textGrob("Hello")) + plot_annotation(theme = theme(plot.margin = margin(.2,.2,.2,.2, "cm"))) & theme(plot.background = element_rect(fill = "yellow")) | wrap_elements(full = r2.s1) & theme(plot.background = element_rect(fill = "blue"))) + plot_layout(widths = c(1,10))) /
  ((wrap_elements(panel = textGrob("Hello")) | wrap_elements(full = r2.s1) & theme(plot.background = element_rect(fill = "blue"))) + plot_layout(widths = c(1,10)))

r2.s1 / r2.s1

plot_grid(p2,q2,r2, ncol = 2)
plot_grid(plot_grid(p2,q2),r2, ncol = 1, labels = "AUTO")
library(cowplot)

plots <- cowplot::align_plots(p2, r2, align = 'v', axis = 'l')
up_row <- plot_grid(plots[[1]], q2, labels = c('B', 'C'), label_size = 12)

plot_grid(up_row,plots[[2]], labels = c('A', ''), label_size = 12, ncol = 1)


align <- cowplot::align_plots(p2,q2, align = "h", axis = "b")
up_row <- plot_grid(p2,q2)

dim <- get_dim(up_row)
bottom_row <- plot_grid(r2)
bottom_align <- set_dim(bottom_row, dim)

fig <- plot_grid(up_row, bottom_align, ncol = 1)
fig
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


