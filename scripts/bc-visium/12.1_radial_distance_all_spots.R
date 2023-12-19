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

# Load seurat object with TCs annotation and beyondcell with allspots
seuratobj.tcs <- readRDS(file = "./results/analysis/seuratobj.therapeutic.clusters.rds")
bc.allspots <- readRDS(file = "./results/analysis/beyondcell_allspots_breastsignature.rds")

# Semla object
semlaobj <- UpdateSeuratForSemla(seuratobj.tcs, 
                                 image_type = "tissue_lowres",
                                 verbose = T)

# Radial distances
semlaobj <- RadialDistance(semlaobj, 
                           column_name = "spot.dual",
                           selected_groups = "TUMOUR")

dual.colors <- c(
  "TUMOUR" = "#c4534e",
  "TME" = "#098cf0")
MapLabels(semlaobj, 
          column_name = "spot.dual", 
          override_plot_dims = TRUE, 
          image_use = NULL, 
          drop_na = TRUE, 
          pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  scale_fill_manual(values = dual.colors) &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

semlaobj@meta.data

MapFeatures(semlaobj, 
            features = "r_dist_TUMOUR", 
            center_zero = TRUE, 
            pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
semlaobj$r_dist_TUMOUR_scaled <- sign(semlaobj$r_dist_TUMOUR)*sqrt(abs(semlaobj$r_dist_TUMOUR))
MapFeatures(semlaobj, 
            features = "r_dist_TUMOUR_scaled", 
            center_zero = TRUE, 
            pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev(),
            override_plot_dims = TRUE)


r.dist <- semlaobj@meta.data$r_dist_TUMOUR_scaled
zebularine <- bc.allspots@normalized["zebularine_CTRP_A01145011", ]

cor.test(x = zebularine, y = r.dist, method = "pearson")

# Transposed enriched matrix
enrich.matrix <- t(bc.allspots@normalized)

# 
radial.dist.tumour <- semlaobj@meta.data$r_dist_TUMOUR

# Realizar la prueba de correlación de Pearson entre cada columna de datos y el vector de distancia
results.cor.pearson <- apply(enrich.matrix, 2, function(col) cor.test(col, radial.dist.tumour, method = "pearson"))

results.df <- data.frame(
  corr = sapply(results.cor.pearson, function(res) res$estimate),
  p.value = sapply(results.cor.pearson, function(res) res$p.value)
)

collapsed.moas <- read_tsv(file = "./data/tsv/collapsed.moas.top.differential.drugs - top.differential.drugs.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$top.diff
names.moas <- levels(factor(collapsed.moas$collapsed.MoAs))
length.moas <- length(names.moas)
col.moas <- c("#db4470",
              "#5cc151",
              "#9b58cf",
              "#9bb932",
              "#5c6fda",
              "#cca939",
              "#dd6fd1",
              "#508f36",
              "#c83f97",
              "#63c385",
              "#8850a1",
              "#da8b2f",
              "#5b70b5",
              "#c45926",
              "#43c4c4",
              "#d5433c",
              "#61a2da",
              "#926a2e",
              "#c190d8",
              "#317341",
              "#e286a5",
              "#4a9f7c",
              "#9c4c77",
              "#abb061",
              "#ae5050",
              "#697329",
              "#df936d")
names(col.moas) <- names.moas

subset.moas <- collapsed.moas %>%
  select(top.diff, preferred.drug.names, collapsed.MoAs)

top.diff <- collapsed.moas$top.diff

results.top.diff <- results.df[top.diff,]
results.top.diff <- results.top.diff %>%
  rownames_to_column("top.diff") %>%
  mutate(top.diff = gsub("\\.cor$", "", top.diff)) 
  
results.top.diff <- right_join(y=results.top.diff, x=subset.moas, by = "top.diff")


results.top.diff.filtered <- results.top.diff %>%
  filter(p.value < 0.05) %>%
  arrange(corr)

ggplot(results.top.diff.filtered, aes(x=corr, y=reorder(preferred.drug.names, corr))) +
  geom_bar(aes(fill = corr), stat = "identity") +
  scale_fill_gradient2(limits = c(-1,1)) +
  xlim(-1,1) +
  labs(fill = "Pearson's correlation") +
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank())
