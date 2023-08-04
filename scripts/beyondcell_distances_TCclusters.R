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

# Load beyondcell object
bc.object <- readRDS("./results/analysis/beyondcellobject.rds")
head(bc.object@meta.data)

# Number of slides
nslides <- length(bc.object@images)

# Number of rows to split the SpatialDimPlots into and height of the final plot
nrows <- nslides/3
h <- 7*nrows

# Get coordinates
coordinates <- bc.object@meta.data %>%
  select(corr_x, corr_y, corr_z)

# Euclidean distances
distances <- as.matrix(dist(coordinates, method = "euclidean"))

# Scale distances
distances.scaled <- distances/min(distances[distances > 0], na.rm = TRUE)


head(bc.object@meta.data)

# Split multiple spots in each cell type
spot.comp.df <- bc.object@meta.data %>%
  select(spot.composition.filter, bc_clusters_res.0.3, mindist.puretumour:mindist.cafs) %>%
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
levels(factor(spot.comp.splited$cell.composition))

plot1 <- spot.comp.splited %>%
  filter(cell.composition == "CAFs") %>%
  #filter(bc_clusters_res.0.3 == 0) %>%
  ggplot(mapping = aes(x = mindist.puretumour)) +
  geom_density(aes(fill = bc_clusters_res.0.3)) +
  facet_grid(rows = vars(bc_clusters_res.0.3),
             #rows = vars(cell.composition),
             scales = "free_y") +
  theme(legend.position = "none") +
  ggtitle("A) CAFs") 
plot1

plot2 <- spot.comp.splited %>%
  #filter(mindist.puretumour != 0) %>%
  #filter(bc_clusters_res.0.3 == 0) %>%
  ggplot() +
  geom_density(aes(x = mindist.puretumour, fill = "Min dist puretumour"), fill="grey") +
  geom_density(mapping = aes(x = mindist.cafs, fill = "Min cafs"), fill = "red", alpha = 0.5) +
  facet_grid(rows = vars(bc_clusters_res.0.3),
             scales = "free_y") +
  ggtitle(paste(" Minimum distance to pure tumour by cell type and Beyondcell TC ")) 
plot2

(plot1 | plot2) /
  (plot1 | plot2)
list.cell.types <- levels(factor(spot.comp.splited$cell.composition))

list.plots <- lapply(list.cell.types, function(x) {
  plot <- spot.comp.splited %>%
    filter(cell.composition == x) %>%
    ggplot(mapping = aes(x = mindist.puretumour)) +
    geom_density(aes(fill = cell.composition)) +
    facet_grid(rows = vars(bc_clusters_res.0.3),
               scales = "free_y") +
    theme(legend.position = "none") +
    ggtitle(x)
  return(plot)
})
(list.plots [[1]] | list.plots[[2]] | list.plots[[3]]) /
  (list.plots [[4]] | list.plots[[5]] | list.plots[[6]]) 


new.df <- spot.comp.splited %>%
  mutate(spot = NULL,
         spot.composition.filter = NULL) %>%
  relocate(cell.composition) %>%
  pivot_longer(cols = c(starts_with('mindist')),
               names_to = "type.distances",
               values_to = "dists")

head(new.df)

plot <- new.df %>%
  #filter(dists != 0) %>%
  ggplot(aes(x=dists)) +
  geom_histogram(aes(fill=type.distances)) +
  facet_grid(cols = vars(bc_clusters_res.0.3),
             rows = vars(type.distances),
             scales = "free_y") +
  theme(legend.position = "none")
plot
