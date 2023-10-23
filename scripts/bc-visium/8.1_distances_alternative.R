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
seuratobj.aligned.alt <- readRDS(file = "./results/analysis/seuratobj.aligned.alternative.rds")
head(seuratobj.aligned.alt@meta.data)
# Number of slides
nslides <- length(seuratobj.aligned@images)
nslides.alt <- length(seuratobj.aligned.alt@images)

# Number of rows to split the SpatialDimPlots into and height of the final plot
nrows <- nslides/3
h <- 7*nrows

# Get coordinates
coordinates <- seuratobj.aligned@meta.data %>%
  select(corr_x, corr_y, corr_z)
coordinates.alt <- seuratobj.aligned.alt@meta.data %>%
  select(corr_x, corr_y, corr_z)

# Euclidean distances
distances <- as.matrix(dist(coordinates, method = "euclidean"))
distances.alt <- as.matrix(dist(coordinates.alt, method = "euclidean"))

# Scale distances
distances.scaled <- distances/min(distances[distances > 0], na.rm = TRUE)


head(seuratobj.aligned@meta.data)
# Split multiple spots in each cell type
spot.comp.df <- seuratobj.aligned@meta.data %>%
  select(spot.composition.collapse) %>%
  rownames_to_column("spot")
spot.comp.splited <- data.frame()
for (i in 1:nrow(spot.comp.df)) {
  spot.composition <- spot.comp.df$spot.composition.collapse[i]
  
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
spot.comp.splited <- spot.comp.splited %>%
  mutate(cell.composition = str_squish(cell.composition))
rownames(spot.comp.splited) <- NULL
head(spot.comp.splited)

spot.comp.df.alt <- seuratobj.aligned.alt@meta.data %>%
  select(spot.composition.collapse) %>%
  rownames_to_column("spot")
spot.comp.splited.alt <- data.frame()
for (i in 1:nrow(spot.comp.df.alt)) {
  spot.composition <- spot.comp.df.alt$spot.composition.collapse[i]
  
  if (grepl(pattern = "\\+", spot.composition)) {
    terms <- strsplit(spot.composition, "\\+")[[1]]
    
    for (term in terms) {
      new_row <- spot.comp.df.alt[i, ]
      new_row$cell.composition <- term
      spot.comp.splited.alt <- rbind(spot.comp.splited.alt, new_row)
    }
  } else {
    spot.comp.df.alt[i, "cell.composition"] <- spot.composition
    spot.comp.splited.alt <- rbind(spot.comp.splited.alt, spot.comp.df.alt[i, ])
  }
}
spot.comp.splited.alt <- spot.comp.splited.alt %>%
  mutate(cell.composition = str_squish(cell.composition))
rownames(spot.comp.splited.alt) <- NULL
head(spot.comp.splited.alt)


# Spot vectors for each cell type
## Function for get it
spot.vector <- function(dataframe, cell.type){
  vector <- dataframe %>%
    filter(cell.composition == cell.type) %>%
    pull(spot)
  return(vector)
}

cell.types <- levels(factor(spot.comp.splited$cell.composition))
cell.types.alt <- levels(factor(spot.comp.splited.alt$cell.composition))
#puretumour.spots <- spot.vector(spot.comp.splited, "PURE TUMOUR")
#cancerepith.spots <- spot.vector(spot.comp.splited, "CANCER EPITHELIAL")
#endothelial.spots <- spot.vector(spot.comp.splited, "ENDOTHELIAL")
#lymphoid.spots <- spot.vector(spot.comp.splited, "LYMPHOID")
#mixed.spots <- spot.vector(spot.comp.splited, "OTHERS")
#myeloid.spots <- spot.vector(spot.comp.splited, "MYELOID")
#cafs.spots <- spot.vector(spot.comp.splited, "CAFS")

list.spots <- lapply(X = cell.types, FUN = function(type){
  cell.types.spot <- spot.vector(spot.comp.splited, type)
  return(cell.types.spot)
})

list.spots.alt <- lapply(X = cell.types.alt, FUN = function(type){
  cell.types.spot <- spot.vector(spot.comp.splited.alt, type)
  return(cell.types.spot)
})

# Get the min distance of each spot to the closer Tumour spot
## Function for get it
celltype.distance <- function(dataframe, celltype.spots){
  cell.dist <- sapply(dataframe$spot, FUN = function(x) {
    dist <- min(distances[x, celltype.spots])
    ifelse(is.na(dist), 0, dist)
  })
  return(cell.dist)
}

#tumour.distance <- celltype.distance(spot.comp.splited, puretumour.spots)
#cancerepith.distance <- celltype.distance(spot.comp.splited, cancerepith.spots)
#endothelial.distance <- celltype.distance(spot.comp.splited, endothelial.spots)
#lymphoid.distance <- celltype.distance(spot.comp.splited, lymphoid.spots)
#mixed.distance <- celltype.distance(spot.comp.splited, mixed.spots)
#myeloid.distance <- celltype.distance(spot.comp.splited, myeloid.spots)
#cafs.distance <- celltype.distance(spot.comp.splited, cafs.spots)

list.distances <- lapply(X = list.spots, FUN = function(spots){
  distances <- celltype.distance(spot.comp.splited, spots)
  return(distances)
})
list.distances.alt <- lapply(X = list.spots.alt, FUN = function(spots){
  distances <- celltype.distance(spot.comp.splited.alt, spots)
  return(distances)
})
#spot.comp.dist <- spot.comp.splited %>%
#  cbind(data.frame(mindist.puretumour = tumour.distance),
#        data.frame(mindist.cancerepith = cancerepith.distance),
#        data.frame(mindist.endothelial = endothelial.distance),
#        data.frame(mindist.lymphoid = lymphoid.distance),
#        data.frame(mindist.mixed = mixed.distance),
#        data.frame(mindist.myeloid = myeloid.distance),
#        data.frame(mindist.cafs = cafs.distance)) 
#head(spot.comp.dist)

distances.df <-  as.data.frame(do.call(cbind, list.distances))
names(distances.df) <- paste0("mindist.",cell.types)
spot.comp.dist <- cbind(spot.comp.splited,distances.df)
rownames(spot.comp.dist) <- NULL
head(spot.comp.dist)

distances.df.alt <-  as.data.frame(do.call(cbind, list.distances.alt))
names(distances.df.alt) <- paste0("mindist.",cell.types.alt)
spot.comp.dist.alt <- cbind(spot.comp.splited.alt,distances.df.alt)
rownames(spot.comp.dist.alt) <- NULL
head(spot.comp.dist.alt)

# Plot density distances
mindist.types <- names(spot.comp.dist)[which(startsWith(names(spot.comp.dist), "mindist"))]
mindist.types.alt <- names(spot.comp.dist.alt)[which(startsWith(names(spot.comp.dist.alt), "mindist"))]

plot.list <- lapply(mindist.types, function(dist) {
  type <- gsub(pattern = "mindist.", replacement = "",dist)
  plot <- spot.comp.dist %>%
    ggplot(mapping = aes(x = !!sym(dist))) +        #Option !!sym() let to pass arguments of a function
    geom_density(aes(fill = cell.composition)) +
    #xlim(0,12000) +
    facet_grid(rows = vars(cell.composition),
               scales = "free_y") +
    ggtitle(paste("Distance to",type,"by cell type")) +
    theme(axis.title.x = element_blank())
})
wrap_plots(plot.list, ncol = 4) + plot_layout(guides = "collect")

plot.list.alt <- lapply(mindist.types.alt, function(dist) {
  type <- gsub(pattern = "mindist.", replacement = "",dist)
  plot <- spot.comp.dist.alt %>%
    ggplot(mapping = aes(x = !!sym(dist))) +        #Option !!sym() let to pass arguments of a function
    geom_density(aes(fill = cell.composition)) +
    #xlim(0,12000) +
    facet_grid(rows = vars(cell.composition),
               scales = "free_y") +
    ggtitle(paste("Distance to",type,"by cell type")) +
    theme(axis.title.x = element_blank())
})
wrap_plots(plot.list.alt, ncol = 4) + plot_layout(guides = "collect")


# Add distance data to metadata
head(spot.comp.dist)
meta.distances <- function(data.frame){
  rwn <- spot.comp.dist %>%
    distinct(spot, .keep_all = TRUE) %>%
    pull(spot)
  metadata <- spot.comp.dist %>%
    distinct(spot, .keep_all = TRUE) %>%
    mutate(spot = NULL,
           spot.composition.collapse = NULL,
           cell.composition = NULL)
  rownames(metadata) <- rwn
  return(metadata)
}
metadata <- meta.distances(spot.comp.dist)
head(metadata)
metadata.alt <- meta.distances(spot.comp.dist.alt)
head(metadata.alt)

seuratobj.distances <- AddMetaData(seuratobj.aligned, metadata = metadata)
head(seuratobj.distances@meta.data)
seuratobj.distances.alt <- AddMetaData(seuratobj.aligned.alt, metadata = metadata.alt)
head(seuratobj.distances.alt)

# Spatial features plots
features <- names(seuratobj.distances@meta.data[grep("mindist.", x= names(seuratobj.distances@meta.data))])
SpatialFeaturePlot(seuratobj.distances, features = features, ncol = 8)

spatial.distances.plots <- lapply(X = features, FUN = function(feature){
  customSpatialFeaturePlot(seuratobj.distances, features = feature, return.list = T)
})
spatial.distances.plots <- lapply(spatial.distances.plots, function(x){
  x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
})
spatial.distances.plots[[1]]
wrap_plots(spatial.distances.plots)
spatial.distances.plots <- lapply(
  customSpatialFeaturePlot(seuratobj.distances, features = "mindist.ENDOTHELIAL", return.list = TRUE),
  FUN = function(x) {
    x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
  }) |>
  wrap_plots(nrow = 1, guides = "collect")

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
