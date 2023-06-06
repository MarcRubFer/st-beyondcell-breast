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
seuratobj.refHER2 <- readRDS(file = "./results/analysis/reference/seuratobj.refHER2.rds")
head(seuratobj.refHER2@meta.data)

# Create metadata %Ribosomal (mitochondrial is already made)
seuratobj.refHER2[["Percent_ribo"]] <- 
  PercentageFeatureSet(object = seuratobj.refHER2, pattern = "^RP[SL][[:digit:]]")
seuratobj.refHER2@meta.data <- seuratobj.refHER2@meta.data %>%
  mutate(Percent_ribo = if_else(is.na(Percent_ribo), 0, Percent_ribo))
head(seuratobj.refHER2@meta.data)

# ViolinPlots before spot filtering
Idents(seuratobj.refHER2) <- "celltype_major"
vln.mt.before <- customVln(seuratobj.refHER2, features = "Percent_mito") +
  ggtitle("Mit Pre-filtering")
vln.rb.before <- customVln(seuratobj.refHER2, features = "Percent_ribo") +
  ggtitle("Rib Pre-filtering") 
vln.count.before <- customVln(seuratobj.refHER2, features = "nCount_RNA") +
  ggtitle("Counts Pre-filtering")
vln.feat.before <- customVln(seuratobj.refHER2, features = "nFeature_RNA") +
  ggtitle("Features Pre-filtering")

pre.violin <- (vln.mt.before | vln.rb.before) / (vln.count.before | vln.feat.before) 
pre.violin <- pre.violin & theme(legend.position = "bottom") 
pre.violin <- pre.violin + plot_layout(guides = "collect")
pre.violin <- pre.violin + 
  plot_annotation(
    title = "Prefiltered QC parameters",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
    )


# params:
#   qc:
#   mt_limits: [0, 15]
# rb_limits: [0, 40]
# count_limits: [250, 50000]
# feat_limits: [200, 7500]
# Filtering seuratobject: mitochondrial percentaje (<10-15%), nCount (always
# >0, depends on status 200-500) and nFeature(always >0, depends on status
# 100-200)
seuratobj.refHER2.filtered <- subset(seuratobj.refHER2, 
                                 subset = (Percent_mito < 15) &
                                   (Percent_ribo < 40) & 
                                   (nCount_RNA > 250 & nCount_RNA < 50000) &  
                                   (nFeature_RNA > 200 & nFeature_RNA <7500))
dim(seuratobj.refHER2)
dim(seuratobj.refHER2.filtered)

# VlnPlots after spot filtering. Y limits adjust to pre-filtered

Idents(seuratobj.refHER2.filtered) <- "celltype_major"
vln.mt.after <- customVln(seuratobj.refHER2.filtered, features = "Percent_mito") + 
  ylim(0,20) +
  ggtitle("Mit Filtered")
vln.rb.after <- customVln(seuratobj.refHER2.filtered, features = "Percent_ribo") +
  ylim(0,60) +
  ggtitle("Rb Filtered")
vln.count.after <- customVln(seuratobj.refHER2.filtered, features = "nCount_RNA") +
  ylim(0, 150000) +
  ggtitle("Count Filtered")
vln.feat.after <- customVln(seuratobj.refHER2.filtered, features = "nFeature_RNA") +
  ylim(0, 10000) +
  ggtitle("Feat Filtered")

post.violin <- (vln.mt.after | vln.rb.after ) / (vln.count.after | vln.feat.after) & theme(legend.position = "bottom")
post.violin <- post.violin + plot_layout(guides = "collect")
post.violin <- post.violin + 
  plot_annotation(
    title = "Postfiltered QC parameters",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )
post.violin

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/reference"), recursive = TRUE)
ggsave(filename = "Violin.plots_prefiltered.reference.png", 
       plot = pre.violin, 
       path = paste0(out.dir,"/plots/reference"))
ggsave(filename = "Violin.plots_postfiltered.reference.png", 
       plot = post.violin, 
       path = paste0(out.dir,"/plots/reference"))

dir.create(path = paste0(out.dir,"/ggplots/reference"), recursive = TRUE)
all.plots <- list(pre.violin, post.violin)
save(all.plots, file = paste0("./results/ggplots/reference/filter_plots.RData"))

# Save data
dir.create(path = paste0(out.dir,"/analysis/reference"), recursive = TRUE)
saveRDS(seuratobj.refHER2.filtered, 
        file = "./results/analysis/reference/seuratobj.refHER2.filtered.rds")
