rm(list = ls())

library("spacexr")
library("Seurat")
library("tidyverse")
library("patchwork")

### FUNCTIONS ###

# Returns Seurat's SpatialFeaturePlot default colourscale
colorscale <- function(n = 100, rev = TRUE) {
  colors <- RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  if (rev) {
    colors <- colors |> 
      rev()
  }
  colors <- colors |>
    grDevices::colorRampPalette()
  return(colors(n = n))
}

### CODE ###

# Load Seurat cluster object
seuratobj.clusters <- readRDS(file = "./results/analysis/seuratobj.clusters.rds")

# Load Spacexr Reference
reference.HER2 <- readRDS(file = "./results/analysis/spacexr.reference.rds")

# Get raw counts and coordinates
counts <- seuratobj.clusters@assays$Spatial@counts
coords <- GetTissueCoordinates(seuratobj.clusters)

colnames(counts)

puck <- SpatialRNA(coords, counts)

myRCTD <- create.RCTD(puck, reference.HER2, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# Cell type proportions for each spot
cell.prop <- myRCTD@results$weights

# Normalize cell proportions to sum 1
cell.prop.normalized <- as.data.frame(normalize_weights(cell.prop))
head(cell.prop.normalized)

# Create categorical cell type
cell.prop.normalized <- cell.prop.normalized %>%
  mutate(Cell.Type = case_when(`Cancer Epithelial` > 0.65 ~ "Tumour", CAFs > 0.20 ~ "CAFs", TRUE ~ "Others"))
head(cell.prop.normalized)

# Add Metadata to seurat object
seuratobj.deconvoluted <- AddMetaData(seuratobj.clusters, 
                                      metadata = cell.prop.normalized)
head(seuratobj.deconvoluted@meta.data)




cell.types <- colnames(cell.prop.normalized) |>
  gsub(pattern = "-", replacement = ".") |>
  gsub(pattern = " ", replacement = ".")

# Plot cell types
l <- SpatialFeaturePlot(seuratobj.deconvoluted, features = cell.types,
                   combine = FALSE)
cell.types.tittles <- gsub(pattern = "\\.", replacement = " ", cell.types)
l2 <- lapply(seq_along(l), FUN = function(x) {
  l[[x]] +
    ggtitle(label = cell.types.tittles[x]) +
    scale_fill_gradientn(colours = colorscale(n = 100),
                          limits = c(0, 1)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
})

cell.type.prop <- wrap_plots(l2, ncol = 3, nrow = 3, guides = "collect") + 
  plot_annotation(
    title = 'Cell type proportion',
    caption = 'Section 1') &
  theme(legend.position = "bottom")
cell.type.prop

ggsave(filename = "cell.type.proportions-Section1.png", 
       plot = cell.type.prop, 
       path = "./results/plots/")

 
tumour.cafs.distr <- SpatialDimPlot(seuratobj.deconvoluted, 
                                    group.by = "Cell.Type") &
  ggtitle(label = "Distribution of Tumour cells and CAFs in section") &
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(1.5, "cm")) 

ggsave(filename = "Tumour.CAFs.distribution.png",
       plot = tumour.cafs.distr,
       path = "./results/plots/")

# Save data
saveRDS(myRCTD, file = "./results/analysis/myRCTD.rds")
