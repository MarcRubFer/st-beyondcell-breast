rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Create ggplot2 theme
blank_x_theme <- theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank(),
                       plot.title = element_text(size = 14, face = "plain"))

# Create patchwork theme
title_theme <- theme(plot.title = element_text(size = 16, face = "bold"))

# Draws custom VlnPlot
customVln <- function(x, features) {
  VlnPlot(x, features = features, pt.size = 0.1) + blank_x_theme
}

#Load RDS
seuratobj <- readRDS(file = "./results/analysis/seuratobj.raw.rds")

# Create metadata about %Mitochondrial and %Ribosomal
seuratobj[["percent.mt"]] <- 
  PercentageFeatureSet(object = seuratobj, pattern = "^MT-")
seuratobj[["percent.rb"]] <- 
  PercentageFeatureSet(object = seuratobj, pattern = "^RP[SL][[:digit:]]")
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(percent.mt = if_else(is.na(percent.mt), 0, percent.mt),
         percent.rb = if_else(is.na(percent.rb), 0, percent.rb))
head(seuratobj@meta.data)

# ViolinPlots before spot filtering
vln.mt.before <- customVln(seuratobj, features = "percent.mt") +
  ggtitle("Mit Pre-filtering")
vln.rb.before <- customVln(seuratobj, features = "percent.rb") +
  ggtitle("Rib Pre-filtering") 
vln.count.before <- customVln(seuratobj, features = "nCount_Spatial") +
  ggtitle("Counts Pre-filtering")
vln.feat.before <- customVln(seuratobj, features = "nFeature_Spatial") +
  ggtitle("Features Pre-filtering")

pre.violin <- (vln.mt.before | vln.rb.before) / (vln.count.before | vln.feat.before) 
pre.violin <- pre.violin & theme(legend.position = "bottom") 
pre.violin <- pre.violin + plot_layout(guides = "collect")
pre.violin

# SpatialFeaturePlots before spot filtering
sf.mt.before <- SpatialFeaturePlot(seuratobj, features = "percent.mt")
sf.rb.before <- SpatialFeaturePlot(seuratobj, features = "percent.rb")
sf.count.before <- SpatialFeaturePlot(seuratobj, features = "nCount_Spatial")
sf.feat.before <- SpatialFeaturePlot(seuratobj, features = "nFeature_Spatial")

pre.spatial <- (sf.mt.before | sf.rb.before) / (sf.count.before | sf.feat.before)

#vln.mt.before + inset_element(sf.mt.before, left = 0.6, bottom = 0.6, right = 1, top = 1)

# Filtering seuratobject: mitochondrial percentaje (<10-15%), nCount (always
# >0, depends on status 200-500) and nFeature(always >0, depends on status
# 100-200)
seuratobj.filtered <- subset(seuratobj, subset = percent.mt < 15 & nCount_Spatial > 500 & nFeature_Spatial > 200)
dim(seuratobj)
dim(seuratobj.filtered)

# VlnPlots after spot filtering
vln.mt.after <- customVln(seuratobj.filtered, features = "percent.mt") + 
  ggtitle("Mit Filtered")
vln.rb.after <- customVln(seuratobj.filtered, features = "percent.rb") +
  ggtitle("Rb Filtered")
vln.count.after <- customVln(seuratobj.filtered, features = "nCount_Spatial") +
  ggtitle("Count Filtered")
vln.feat.after <- customVln(seuratobj.filtered, features = "nFeature_Spatial") +
  ggtitle("Feat Filtered")

post.violin <- (vln.mt.after | vln.rb.after ) / (vln.count.after | vln.feat.after) & theme(legend.position = "bottom")
post.violin <- post.violin + plot_layout(guides = "collect")

# SpatialFeaturePlots after spot filtering
sf.mt.after <- SpatialFeaturePlot(seuratobj.filtered, features = "percent.mt")
sf.rb.after <- SpatialFeaturePlot(seuratobj.filtered, features = "percent.rb")
sf.count.after <- SpatialFeaturePlot(seuratobj.filtered, features = "nCount_Spatial")
sf.feat.after <- SpatialFeaturePlot(seuratobj.filtered, features = "nFeature_Spatial")

post.spatial <- (sf.mt.after | sf.rb.after ) / (sf.count.after | sf.feat.after)

# Save data
dir.create(path = paste0(out.dir,"/analysis"), recursive = TRUE)
saveRDS(seuratobj.filtered, file = "./results/analysis/seuratobj.filtered.rds")

dir.create(path = paste0(out.dir,"/ggplots"), recursive = TRUE)
all.plots <- list(pre.violin, pre.spatial, post.violin, post.spatial)
save(all.plots, file = paste0("./results/ggplots/filter_plots.RData"))

dir.create(path = paste0(out.dir,"/filtering"), recursive = TRUE)
ggsave(filename = paste0(out.dir,"/filtering/patch.prefiltering.violin.pdf"),
       plot = pre.violin)
ggsave(filename = paste0(out.dir,"/filtering/patch.prefiltering.violin.jpg"),
       plot = pre.violin)
ggsave(filename = paste0(out.dir,"/filtering/patch.prefiltering.spatial.pdf"),
       plot = pre.spatial)
ggsave(filename = paste0(out.dir,"/filtering/patch.postfiltering.violin.pdf"),
       plot = post.violin)
ggsave(filename = paste0(out.dir,"/filtering/patch.postfiltering.spatial.pdf"),
       plot = post.spatial)
