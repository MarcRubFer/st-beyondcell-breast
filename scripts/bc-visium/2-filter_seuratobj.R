rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")

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

VlnPlot(seuratobj, features = "nCount_Spatial", pt.size = 0.1)
names(seuratobj)
