#install.packages("remotes")
#library("remotes")
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)


seuratobj.clusters <- readRDS(file = "./results/analysis/seuratobj.clusters.rds")
head(seuratobj.clusters@meta.data)

df <- seuratobj.clusters@meta.data %>%
  select(spot.composition.collapse, SCT_snn_res.0.2)

head(df)

df <- df %>%
  make_long(spot.composition.collapse, SCT_snn_res.0.2)

pl <- ggplot(df, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl <- pl +geom_sankey(flow.alpha = 0.5
                      , node.color = "black"
                      ,show.legend = T)
pl <- pl +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = 1)
pl

df.alt <- seuratobj.clusters.alt@meta.data %>%
  select(spot.composition.collapse, SCT_snn_res.0.2)

head(df)

df.alt <- df.alt %>%
  make_long(spot.composition.collapse, SCT_snn_res.0.2)

pl.alt <- ggplot(df.alt, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl.alt <- pl.alt +geom_sankey(flow.alpha = 0.5
                      , node.color = "black"
                      ,show.legend = T)
pl.alt <- pl.alt +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = 1)
pl.alt

(pl | pl.alt) & theme(legend.position = "none")

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "black") +
  scale_fill_viridis_d() +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")


(SpatialDimPlot(seuratobj.clusters, group.by = "SCT_snn_res.0.2") | SpatialDimPlot(seuratobj.clusters.alt, group.by = "SCT_snn_res.0.2")) / ((pl | pl.alt) & theme(legend.position = "none"))
                
                