dim(drugs.matrix.TME2)

metadata <- bc.ranked.TME@meta.data %>%
  select(spot.collapse, TCs_res.0.3) %>%
  rownames_to_column("spot")
head(metadata)

enrichment <- as.data.frame(t(drugs.matrix.TME2)) %>%
  rownames_to_column("spot")
head(enrichment)

data.raw <- right_join(metadata,enrichment, by = "spot")
head(data.raw)

data <- data.raw %>%
  pivot_longer(cols = -c(spot,spot.collapse,TCs_res.0.3), names_to = "signature", values_to = "enrichment")
head(data)

sigs.UP <- data %>%
  filter(TCs_res.0.3 == "TC-2",
         enrichment > 0) %>%
  select(signature) %>%
  unique() %>%
  pull()

sigs.DOWN <- data %>%
  filter(TCs_res.0.3 == "TC-2",
         enrichment < 0) %>%
  select(signature) %>%
  unique() %>%
  pull()

violin.UP <- ggplot(subset(data, signature %in% sigs.UP), aes(x = TCs_res.0.3, y = enrichment)) 
violin.UP + geom_violin(aes(fill = spot.collapse))

violin.DOWN <- ggplot(subset(data, signature %in% sigs.DOWN), aes(x = TCs_res.0.3, y = enrichment)) 
violin.DOWN + geom_violin(aes(fill = spot.collapse))


TC1.Down.in.UP <- data %>%
  filter(signature %in% sigs.UP,
         TCs_res.0.3 == "TC-1",
         enrichment < 0) %>%
  select(spot) %>%
  unique() %>%
  pull()

TC1.UP.in.DOWN <- data %>%
  filter(signature %in% sigs.DOWN,
         TCs_res.0.3 == "TC-1",
         enrichment > 0) %>%
  select(spot) %>%
  unique() %>%
  pull()

categories.down <- ggplot(subset(data.raw, spot %in% TC1.Down.in.UP), aes(x = spot.collapse)) +
  geom_bar() +
  ylim(NA,200) +
  ggtitle(label = "Categories in spots down in sigs UP")
categories.up <- ggplot(subset(data.raw, spot %in% TC1.UP.in.DOWN), aes(x = spot.collapse)) +
  geom_bar() +
  ylim(NA,200) +
  ggtitle(label = "Categories in spots up in sigs DOWN")

categories.figures <- categories.down | categories.up
categories.figures

intersect.spots <- intersect(TC1.Down.in.UP,TC1.UP.in.DOWN)


seurat.TME <- readRDS("./results/analysis/seuratobj.TME-TCs.rds")
SpatialDimPlot(seurat.TME, cells.highlight = TC1.Down.in.UP)

coords.DOWN.in.UP <- seurat.TME@meta.data %>%
  rownames_to_column("spot") %>%
  filter(spot %in% TC1.Down.in.UP)
coords.DOWN.in.UP.plot <- ggplot(coords.DOWN.in.UP, aes(x = corr_y, y = corr_x)) +
  geom_point(aes(col = spot.collapse)) +
  scale_y_reverse() +
  ggtitle(label = "Spots DOWN-enrich in UP-sigs", subtitle = "At least in one sig")

coords.UP.in.DOWN <- seurat.TME@meta.data %>%
  rownames_to_column("spot") %>%
  filter(spot %in% TC1.UP.in.DOWN)
coords.UP.in.DOWN.plot <- ggplot(coords.UP.in.DOWN, aes(x = corr_y, y = corr_x)) +
  geom_point(aes(col = spot.collapse)) +
  scale_y_reverse() +
  ggtitle(label = "Spots UP-enrich in DOWN-sigs", subtitle = "At least in one sig")

coords.intersect <- seurat.TME@meta.data %>%
  rownames_to_column("spot") %>%
  filter(spot %in% intersect.spots)
coords.intersect.plot <- ggplot(coords.intersect, aes(x = corr_y, y = corr_x)) +
  geom_point(aes(col = spot.collapse)) +
  scale_y_reverse() +
  ggtitle(label = "Intersect spots UP_DOWN", subtitle = "At least in one sig")

panel <- (coords.DOWN.in.UP.plot / plot_spacer() / coords.UP.in.DOWN.plot) + plot_layout(heights = c(1,0.3,1))
figure.spot <- (panel | coords.intersect.plot | SpatialDimPlot(seurat.TME,group.by = "spot.collapse", ncol = 1)) + plot_layout(widths = c(2,4,1))

ggsave(filename = "Coords_spots_TC1.png",
       plot = figure.spot,
       path = "./results/plots/TC_TME_analysis/")

cell.type.prop <- coords.intersect %>%
  select(spot, B.cells:T.cells)
head(cell.type.prop)
cell.type.prop <- cell.type.prop %>%
  pivot_longer(cols = -spot, names_to = "cell.type", values_to = "proportion")
head(cell.type.prop)
ggplot(cell.type.prop, aes(x = cell.type, y = proportion)) +
  geom_boxplot()
