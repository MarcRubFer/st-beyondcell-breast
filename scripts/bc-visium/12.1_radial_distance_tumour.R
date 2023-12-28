rm(list = ls())

library(Seurat)
library(beyondcell)
library(semla)
library(tibble)
library("tidyverse")

library(ggplot2)
library(patchwork)
library(scico)
library(tidyr)
library(dplyr)

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

set.seed(1)

# Set dual colors Tumour/TME
dual.colors <- c(
  "TUMOUR" = "#c4534e",
  "TME" = "#098cf0")


seurat.tumor <- readRDS("./results/analysis/seuratobj.Tumour-TCs.rds")

seurat.tumor@meta.data

SpatialDimPlot(seurat.tumor, group.by = "spot.dual", cols = dual.colors) / SpatialDimPlot(seurat.tumor, group.by = "TCs_res.0.3")

# Create a Semla object from Seurat object
semlaobj.tumour <- UpdateSeuratForSemla(seurat.tumor, 
                                 image_type = "tissue_lowres",
                                 verbose = T)
semlaobj.tumour <- LoadImages(semlaobj.tumour)
# Calculate Radial distances 
semlaobj.tumour <- RadialDistance(semlaobj.tumour, 
                           column_name = "spot.dual",
                           selected_groups = "TUMOUR")
# Plot distribution of Tumour TME
tumour.tme.spots <- MapLabels(semlaobj.tumour, 
                              column_name = "spot.dual", 
                              override_plot_dims = TRUE, 
                              image_use = "raw", 
                              drop_na = TRUE, 
                              pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  scale_fill_manual(values = dual.colors) &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

ggsave(filename = "tumour_tme_spots.png",
       plot = tumour.tme.spots,
       path = "./results/plots/TC_Tumour_analysis/",
       scale = 0.45)
ggsave(filename = "tumour_tme_spots.svg",
       plot = tumour.tme.spots,
       path = "./results/plots/TC_Tumour_analysis/",
       scale = 0.45)

semlaobj.tumour@meta.data
# Scale radial distances with square root
semlaobj.tumour$r_dist_TUMOUR_scaled <- sign(semlaobj.tumour$r_dist_TUMOUR)*sqrt(abs(semlaobj.tumour$r_dist_TUMOUR))
radial.dist <- MapFeatures(semlaobj.tumour, 
                           features = "r_dist_TUMOUR_scaled", 
                           center_zero = TRUE, 
                           pt_size = 2, 
                           colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev(),
                           override_plot_dims = TRUE,
                           image_use = "raw")
ggsave(filename = "radial_dist.png",
       plot = radial.dist,
       path = "./results/plots/TC_Tumour_analysis/",
       scale = 0.5)
ggsave(filename = "radial_dist.svg",
       plot = radial.dist,
       path = "./results/plots/TC_Tumour_analysis/",
       scale = 0.5)

# Load beyondcell object for TC-Tumour
bc.tumour <- readRDS("./results/analysis/bc_ranked_Tumour.rds")

# Transposed enriched matrix
enrich.matrix <- t(bc.tumour@normalized)
# 
radial.dist.tumour <- semlaobj.tumour@meta.data$r_dist_TUMOUR

# Test Pearson's correlation between each enrichment column and distance
results.cor.pearson <- apply(enrich.matrix, 2, function(col) cor.test(col, radial.dist.tumour, method = "pearson"))

results.df <- data.frame(
  corr = sapply(results.cor.pearson, function(res) res$estimate),
  p.value = sapply(results.cor.pearson, function(res) res$p.value)
)

collapsed.moas.Tumour <- read_tsv(file = "./data/tsv/top.differential.drugs.TCs.tumour.tsv")
collapsed.moas.Tumour <- as.data.frame(collapsed.moas.Tumour)
rownames(collapsed.moas.Tumour) <- collapsed.moas.Tumour$top.diff

subset.moas <- collapsed.moas.Tumour %>%
  select(top.diff, preferred.drug.names, collapsed.MoAs)

top.diff <- collapsed.moas.Tumour$top.diff

results.top.diff <- results.df[top.diff,]
results.top.diff <- results.top.diff %>%
  rownames_to_column("top.diff") %>%
  mutate(top.diff = gsub("\\.cor$", "", top.diff)) 

results.top.diff <- right_join(y=results.top.diff, x=subset.moas, by = "top.diff")

results.top.diff.filtered <- results.top.diff %>%
  filter(p.value < 0.05) %>%
  arrange(corr)

correlation.plot <- ggplot(results.top.diff.filtered, aes(x=corr, y=reorder(preferred.drug.names, corr))) +
  geom_bar(aes(fill = corr), stat = "identity") +
  scale_fill_gradient2(limits = c(-1,1)) +
  xlim(-1,1) +
  labs(fill = "Pearson's correlation") +
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(filename = "correlation_plot.png",
       plot = correlation.plot,
       path = "./results/plots/TC_Tumour_analysis/")
ggsave(filename = "correlation_plot.svg",
       plot = correlation.plot,
       path = "./results/plots/TC_Tumour_analysis/")

top1_simvastatin <- bcSignatures(bc.tumour, spatial = T, mfrow = c(1,2), signatures = list(values = "simvastatin_CTRP_K22134346"))
ggsave(filename = "simvastatin.svg",
       plot = top1_simvastatin,
       path = "./results/plots/TC_Tumour_analysis/")
top2_IEM1754 <- bcSignatures(bc.tumour, spatial = T, mfrow = c(1,2), signatures = list(values = "IEM1754_PRISM_K80060353"))
ggsave(filename = "IEM1754.svg",
       plot = top2_IEM1754,
       path = "./results/plots/TC_Tumour_analysis/")
top1_VTP27999 <- bcSignatures(bc.tumour, spatial = T, mfrow = c(1,2), signatures = list(values = "VTP-27999_PRISM_K81694556"))
ggsave(filename = "VTP-27999.svg",
       plot = top1_VTP27999,
       path = "./results/plots/TC_Tumour_analysis/")
top2_valnemulin <- bcSignatures(bc.tumour, spatial = T, mfrow = c(1,2), signatures = list(values = "valnemulin_PRISM_K33813875"))
ggsave(filename = "valnemulin.svg",
       plot = top2_valnemulin,
       path = "./results/plots/TC_Tumour_analysis/")

SpatialDimPlot(seurat.tumor, group.by = "TCs_res.0.3") / bcSignatures(bc.tumour, spatial = T, mfrow = c(1,2), signatures = list(values = "valnemulin_PRISM_K33813875"))

