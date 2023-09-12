rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("ggsci")

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

# Draws custom SpatialFeaturePlots
customSpatialFeaturePlot <- function(x, features, limits = NULL, 
                                     legend.title = NULL, share.scale = TRUE, 
                                     return.list = FALSE, ...) {
  nslides <- length(x@images)
  nrows <- nslides/3
  if (nrows < 1) {
    nrows <- 1
  }
  titles <- rep(names(x@images), length(features))
  if (is.null(limits) & share.scale) {
    values <- x@meta.data %>%
      select(all_of(features))
    limits.range <- pretty(na.omit(values[, , drop = TRUE]))
    limits <- c(min(limits.range), max(limits.range))
  }
  if (is.null(limits)) {
    scale <- NULL
  } else {
    scale <- scale_fill_gradientn(colours = colorscale(n = 100, ...),
                                  limits = limits)
  }
  if (is.null(legend.title)) {
    legend <- NULL
  } else {
    legend <- labs(fill = legend.title)
    titles <- paste(titles, rep(features, each = nslides), sep = " - ")
  }
  lplots <- SpatialFeaturePlot(x, features = features, pt.size.factor = 1, 
                               combine = FALSE, image.alpha = 0)
  lplots <- suppressMessages(
    lapply(seq_along(lplots), FUN = function(i) {
      lplots[[i]] + 
        scale + 
        ggtitle(titles[[i]]) +
        legend
    }))
  if(!return.list) {
    lplots <- lapply(seq(from = 1, to = length(lplots), by = nslides), 
                     FUN = function(i) {
                       wrap_plots(lplots[i:(i+nslides-1)], nrow = nrows) + 
                         plot_layout(guides = "collect")
                     })
  }
  return(lplots)
}

## CODE ##
# Random seed
set.seed(1)

# Load seuratobject
seuratobj.aligned <- readRDS(file = "./results/analysis/seuratobj.aligned.rds")
head(seuratobj.aligned@meta.data)

# Number of slides
nslides <- length(seuratobj.aligned@images)

# Number of rows to split the SpatialDimPlots into and height of the final plot
nrows <- nslides/3
h <- 7*nrows

# Get coordinates
coordinates <- seuratobj.aligned@meta.data %>%
  select(corr_x, corr_y, corr_z)

# Euclidean distances
distances <- as.matrix(dist(coordinates, method = "euclidean"))

# Scale distances
distances.scaled <- distances/min(distances[distances > 0], na.rm = TRUE)


head(seuratobj.aligned@meta.data)
# Split multiple spots in each cell type
spot.comp.df <- seuratobj.aligned@meta.data %>%
  select(spot.composition.filter) %>%
  rownames_to_column("spot")
spot.comp.splited <- data.frame()
for (i in 1:nrow(spot.comp.df)) {
  spot.composition <- spot.comp.df$spot.composition.filter[i]
  
  if (grepl(pattern = "\\+", spot.composition)) {
    terms <- strsplit(spot.composition, "\\+")[[1]]
    
    for (term in terms) {
      new_row <- spot.comp.df[i, ]
      new_row$cell.composition <- term
      spot.comp.splited <- rbind(spot.comp.splited, new_row)
    }
  } else {
    spot.comp.df[i, "cell.composition"] <- spot.composition
    spot.comp.splited <- rbind(spot.comp.splited, spot.comp.df[i, ])
  }
}
spot.comp.splited
rownames(spot.comp.splited) <- NULL
head(spot.comp.splited)

# Spot vectors for each cell type
## Function for get it
spot.vector <- function(dataframe, cell.type){
  vector <- dataframe %>%
    filter(cell.composition == cell.type) %>%
    pull(spot)
  return(vector)
}

levels(factor(spot.comp.splited$cell.composition))

puretumour.spots <- spot.vector(spot.comp.splited, "Pure_Tumour")
cancerepith.spots <- spot.vector(spot.comp.splited, "Cancer.Epithelial")
endothelial.spots <- spot.vector(spot.comp.splited, "Endothelial")
lymphoid.spots <- spot.vector(spot.comp.splited, "Lymphoid")
mixed.spots <- spot.vector(spot.comp.splited, "Mixed")
myeloid.spots <- spot.vector(spot.comp.splited, "Myeloid")
cafs.spots <- spot.vector(spot.comp.splited, "CAFs")

# Get the min distance of each spot to the closer Tumour spot
## Function for get it
celltype.distance <- function(dataframe, celltype.spots){
  cell.dist <- sapply(dataframe$spot, FUN = function(x) {
    dist <- min(distances[x, celltype.spots])
    ifelse(is.na(dist), 0, dist)
  })
  return(cell.dist)
}

tumour.distance <- celltype.distance(spot.comp.splited, puretumour.spots)
cancerepith.distance <- celltype.distance(spot.comp.splited, cancerepith.spots)
endothelial.distance <- celltype.distance(spot.comp.splited, endothelial.spots)
lymphoid.distance <- celltype.distance(spot.comp.splited, lymphoid.spots)
mixed.distance <- celltype.distance(spot.comp.splited, mixed.spots)
myeloid.distance <- celltype.distance(spot.comp.splited, myeloid.spots)
cafs.distance <- celltype.distance(spot.comp.splited, cafs.spots)

spot.comp.dist <- spot.comp.splited %>%
  cbind(data.frame(mindist.puretumour = tumour.distance),
        data.frame(mindist.cancerepith = cancerepith.distance),
        data.frame(mindist.endothelial = endothelial.distance),
        data.frame(mindist.lymphoid = lymphoid.distance),
        data.frame(mindist.mixed = mixed.distance),
        data.frame(mindist.myeloid = myeloid.distance),
        data.frame(mindist.cafs = cafs.distance)) 
head(spot.comp.dist)

# Plot density distances
mindist.types <- names(spot.comp.dist)[which(startsWith(names(spot.comp.dist), "mindist"))]

plot.list <- lapply(mindist.types, function(dist) {
  plot <- spot.comp.dist %>%
    ggplot(mapping = aes(x = !!sym(dist))) +
    geom_density(aes(fill = cell.composition)) +
    facet_grid(rows = vars(cell.composition),
               scales = "free_y") +
    ggtitle(paste(dist,"by cell type"))
})
plot.list[[1]] +
  xlim(c(mean(400), 1000))

combined_plot <- plot.list[[1]]
for (i in 2:length(plot.list)) {
  combined_plot <- combined_plot + plot.list[[i]]
}

combined_plot <- combined_plot + plot_layout(ncol = 2)  
combined_plot

# Add distance data to metadata
head(spot.comp.dist)
rwn <- spot.comp.dist %>%
  distinct(spot, .keep_all = TRUE) %>%
  pull(spot)
metadata <- spot.comp.dist %>%
  distinct(spot, .keep_all = TRUE) %>%
  mutate(spot = NULL,
         spot.composition.filter = NULL,
         cell.composition = NULL)
rownames(metadata) <- rwn
head(metadata)

seuratobj.distances <- AddMetaData(seuratobj.aligned, metadata = metadata)
head(seuratobj.distances@meta.data)


# Spatial features plots
spatial.puretumour.distance <- lapply(
  customSpatialFeaturePlot(seuratobj.distances, features = "mindist.puretumour", return.list = TRUE),
  FUN = function(x) {
    x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
  }) |>
  wrap_plots(nrow = 1, guides = "collect")
spatial.puretumour.distance

spatial.cancerep.distance <- lapply(
  customSpatialFeaturePlot(seuratobj.distances, features = "mindist.cancerepith", return.list = TRUE),
  FUN = function(x) {
    x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
  }) |>
  wrap_plots(nrow = 1, guides = "collect")
spatial.cancerep.distance

spatial.endothelial.distance <- lapply(
  customSpatialFeaturePlot(seuratobj.distances, features = "mindist.endothelial", return.list = TRUE),
  FUN = function(x) {
    x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
  }) |>
  wrap_plots(nrow = 1, guides = "collect")
spatial.endothelial.distance

spatial.lymphoid.distance <- lapply(
  customSpatialFeaturePlot(seuratobj.distances, features = "mindist.lymphoid", return.list = TRUE),
  FUN = function(x) {
    x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
  }) |>
  wrap_plots(nrow = 1, guides = "collect")
spatial.lymphoid.distance

spatial.myeloid.distance <- lapply(
  customSpatialFeaturePlot(seuratobj.distances, features = "mindist.myeloid", return.list = TRUE),
  FUN = function(x) {
    x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
  }) |>
  wrap_plots(nrow = 1, guides = "collect")
spatial.myeloid.distance

plot.list[[6]] +
  xlim(c(1000,2500))

# Save data
saveRDS(seuratobj.distances, file = "./results/analysis/seuratobj.distances.rds")
