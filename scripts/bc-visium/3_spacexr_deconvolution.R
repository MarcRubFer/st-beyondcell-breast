rm(list = ls())

library("spacexr")
library("Seurat")
library("tidyverse")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

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
# Remember: seurat object origin is merged (section 1 and 2)
seuratobj.filtered <- readRDS(file = "./results/analysis/seuratobj.filtered.rds")

# Load Spacexr Reference
reference.HER2 <- readRDS(file = "./results/analysis/reference/spacexr.reference.rds")

# Get raw counts and coordinates from each Section from merged seurat object
# lapply return a list with deconvolution

cell.prop.normalized <- lapply(c("Section1", "Section2"), function(section){
  sc.cells <- seuratobj.filtered@meta.data %>%
    filter(orig.ident == section) %>%
    rownames()
  counts <- seuratobj.filtered@assays$Spatial@counts[,sc.cells]
  coords <- GetTissueCoordinates(seuratobj.filtered, image = section)
  
  # Create SpatialRNA object
  puck <- SpatialRNA(coords, counts)
  
  # Create RCTD object
  myRCTD <- create.RCTD(puck, reference.HER2, max_cores = 1)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  
  # Cell type proportions for each spot
  cell.prop <- myRCTD@results$weights
  
  # Normalize cell proportions to sum 1
  cell.prop.normalized <- as.data.frame(normalize_weights(cell.prop))
  return(cell.prop.normalized)
})

# Merge dataframes on list to a unique dataframe
cell.prop.normalized <- cell.prop.normalized %>%
  bind_rows()

cell.types <- colnames(cell.prop.normalized) |>
  gsub(pattern = "-", replacement = ".") |>
  gsub(pattern = " ", replacement = ".")
cell.types.tittles <- gsub(pattern = "\\.", replacement = " ", cell.types)
cell.types.tittles <- rep(cell.types.tittles, each = 2)

# Create categorical cell types
cell.prop.normalized <- cell.prop.normalized %>%
  mutate(Cell.Type = case_when(`Cancer Epithelial` > 0.65 ~ "Tumour", CAFs > 0.20 ~ "CAFs", TRUE ~ "Others"))
head(cell.prop.normalized)

colnames(cell.prop.normalized) <- gsub(pattern = "-",
                                       replacement = ".",
                                       x = colnames(cell.prop.normalized))
colnames(cell.prop.normalized) <- gsub(pattern = " ",
                                       replacement = ".",
                                       x = colnames(cell.prop.normalized))
head(cell.prop.normalized)
# Add Metadata to seurat object
seuratobj.deconvoluted <- AddMetaData(seuratobj.filtered, 
                                      metadata = cell.prop.normalized)

head(seuratobj.deconvoluted@meta.data)

# Plot cell types
# Cell types proportion plots
l <- SpatialFeaturePlot(seuratobj.deconvoluted, 
                        features = cell.types,
                        combine = FALSE)

l2 <- lapply(seq_along(l), FUN = function(x) {
  l[[x]] +
    ggtitle(label = cell.types.tittles[x]) +
    scale_fill_gradientn(colours = colorscale(n = 100),
                          limits = c(0, 1)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
})

cell.type.prop <- wrap_plots(l2, ncol = 6, nrow = 3, guides = "collect") + 
  plot_annotation(
    title = 'Cell type proportion',
    caption = 'Section 1') &
  theme(legend.position = "bottom")
cell.type.prop


# Distribution of tumour cells, CAFs and others
tumour.cafs.distr <- SpatialDimPlot(seuratobj.deconvoluted, 
                                    group.by = "Cell.Type") &
  ggtitle(label = "Distribution of Tumour cells and CAFs in section") &
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(1.5, "cm")) 
tumour.cafs.distr

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/deconvolution"), recursive = TRUE)
ggsave(filename = "cell.type.proportions.png", 
       plot = cell.type.prop, 
       path = paste0(out.dir,"/plots/deconvolution"))
ggsave(filename = "Tumour.CAFs.distribution.png",
       plot = tumour.cafs.distr,
       path = paste0(out.dir,"/plots/deconvolution"))

all.plots <- list(cell.type.prop, tumour.cafs.distr)
save(all.plots, file = paste0("./results/ggplots/deconvolution_plots.RData"))

# Save data
#saveRDS(myRCTD, file = "./results/analysis/myRCTD.rds")
saveRDS(seuratobj.deconvoluted, 
        file = "./results/analysis/seuratobj.deconvoluted.rds")
