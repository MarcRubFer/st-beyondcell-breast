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

# Variables
colors <- toupper(c(#"#4fafe3",
  "MYELOID" ="#be7adc", #violet
  "CAFS" = "#dec36f", #ocre
  "ENDOTHELIAL" = "#549f42", #green
  "LYMPHOID" = "#f1703a", #orange
  "OTHERS" ="#79696B", #grey
  "TUMOUR" = "#c4534e")) #dark.red
  #"PURE TUMOUR" = "#eb1d34")) # light.red
names(colors)
colors <- colors[order(names(colors))]
colors

############
## CODE ##
# Random seed
set.seed(1)

# Load seuratobject
seuratobj.deconvoluted <- readRDS(file = "./results/analysis/seuratobj.deconvoluted.rds")
head(seuratobj.deconvoluted@meta.data)


# Subset cell types proportions metadata
cell.proportions <- seuratobj.deconvoluted@meta.data %>%
  select(B.cells:T.cells) %>%
  mutate(
    spot.composition = case_when(
      Cancer.Epithelial >= 0.65 ~ "TUMOUR",
      #Cancer.Epithelial > 0.8 ~ "PURE TUMOUR",
      #Cancer.Epithelial >= 0.65 & Cancer.Epithelial <= 0.8 ~ "TUMOUR",
      TRUE ~ toupper(names(select(., -Cancer.Epithelial))[max.col(select(., -Cancer.Epithelial), "first")])
      #TRUE ~ names(.)[max.col(., "first")]
    ),
    spot.collapse = case_when(
      spot.composition == "B.CELLS" | spot.composition =="T.CELLS" ~ "LYMPHOID",
      spot.composition == "NORMAL.EPITHELIAL" | spot.composition == "PVL" ~ "OTHERS",
      TRUE ~ spot.composition
    )
  )
head(cell.proportions)

# Add metadata about cell spot composition
seuratobj.spotcomp <- AddMetaData(seuratobj.deconvoluted, metadata = cell.proportions)
head(seuratobj.spotcomp)

# BarPlot of number of spots by category


num.spots <- ggplot(cell.proportions, aes(y = spot.collapse, fill = spot.collapse)) +
  geom_bar(color = "#32382c") +
  geom_text(stat = "count", aes(label = after_stat(count), hjust = -0.3), color = "black", size = 3) +
  xlab("Number of spots") +
  ylab(NULL) +
  scale_x_continuous(limits = c(0,5000)) +
  scale_y_discrete(limits = rev(levels(as.factor(cell.proportions$spot.collapse)))) +
  scale_fill_manual(guide = guide_legend(title = "Categories"),
                    values = colors) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank())

# Mostrar el gráfico
num.spots
ggsave(filename = "number_spots_by_category.png",
       plot = num.spots,
       path = "./results/plots/spot_composition/")

# Boxplots of cell type proportion by category
df.boxplot <- cell.proportions %>%
  select(-spot.composition) %>%
  pivot_longer(cols = B.cells:T.cells, names_to = "cell.type", values_to = "prop.cell.type") %>%
  mutate(cell.type = as.factor(cell.type),
         prop.cell.type = round(prop.cell.type, digits = 2))
head(df.boxplot)

boxp.props <- df.boxplot %>%
  mutate(cell.type = gsub(pattern = "([B,T])\\.", replacement = "\\1-", cell.type),
         cell.type = gsub(pattern = "\\.", replacement = " ", cell.type),
         cell.type = toupper(cell.type)) %>%
  ggplot(aes(x=cell.type, y=prop.cell.type)) +
  geom_boxplot(aes(fill = cell.type)) +
  guides(fill = guide_legend(title = "Cell types")) +
  ylab(label = "Cell type proportion") +
  facet_grid(~spot.collapse) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        )
## For colour of each title in facet_grid we have to use next code:
gt<-ggplot_gtable(ggplot_build(boxp.props))
strips <- which(startsWith(gt$layout$name,'strip'))
for (s in seq_along(strips)) {
  gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- colors[s]
}

boxp.props <- as.ggplot(gt)
boxp.props 





spatial.distrib.spotcomp <- SpatialDimPlot(seuratobj.spotcomp, group.by = "spot.collapse", cols = colors) + plot_layout(guides = "collect") &
  guides(fill = guide_legend(title = "Categories"))

bottom <- plot_grid(num.spots,NULL,spatial.distrib.spotcomp, ncol = 3, labels = c("B","","C"), rel_widths = c(1,0.01,1))
bottom
figure1 <- plot_grid(boxp.props,NULL,bottom, ncol = 1, labels = c("A","",""),rel_heights = c(1,0.01,1))
figure1

ggsave(filename = "Figure1.png",
       plot = figure1,
       path = "./results/plots/Figures_TFM/")
ggsave(filename = "Figure1.svg",
       plot = figure1,
       path = "./results/plots/Figures_TFM/")

# Save data
saveRDS(seuratobj.spotcomp, file = "./results/analysis/seuratobj.spotcomp.rds")
