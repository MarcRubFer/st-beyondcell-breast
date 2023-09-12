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

# Get the Tumour/Non-tumour spots
metadata <- seuratobj.aligned@meta.data %>%
  select(spot.composition.filter)

tumour.spots <- metadata %>%
  filter(spot.composition.filter == "Pure_Tumour") %>%
  rownames_to_column("spots") %>%
  pull(spots)

non.tumour.spots <- metadata %>%
  filter(spot.composition.filter != "Pure_Tumour") %>%
  rownames_to_column("spots") %>%
  pull(spots)

# Get the min distance of each spot to the closer Tumour spot
tumour.distance <- sapply(rownames(metadata), FUN = function(x) {
  dist <- min(distances[x, tumour.spots])
  ifelse(is.na(dist), 0, dist)
})

# For the Tumour spots, get the min distance to non-Tumour spots
tme.distance <- sapply(rownames(metadata), FUN = function(x) {
  dist <- min(distances[x, non.tumour.spots])
  ifelse(is.na(dist), 0, dist)
})

# New metadata
metadata2 <- metadata %>%
  cbind(data.frame(mindist.tumour = tumour.distance)) %>%
  cbind(data.frame(mindist.tme = tme.distance)) %>%
  mutate(dist = case_when(spot.composition.filter == "Pure_Tumour" ~ mindist.tme,
                          TRUE ~ - mindist.tumour)) 
head(metadata2)

# Metadata for plots
metadata3 <- metadata2 %>%
  rownames_to_column("spot")
## Split variable spot.composition in each term
new_df <- data.frame()
for (i in 1:nrow(metadata3)) {
  spot.composition <- metadata3$spot.composition.filter[i]
  
  if (grepl(pattern = "\\+", spot.composition)) {
    terms <- strsplit(spot.composition, "\\+")[[1]]
    
    for (term in terms) {
      new_row <- metadata3[i, ]
      new_row$cell.composition <- term
      new_df <- rbind(new_df, new_row)
    }
  } else {
    metadata3[i, "cell.composition"] <- spot.composition
    new_df <- rbind(new_df, metadata3[i, ])
  }
}
new_df
rownames(new_df) <- NULL
head(new_df)


ggplot(data = new_df, mapping = aes(x = factor(cell.composition), y = mindist.tumour)) +
  geom_boxplot(aes(factor(cell.composition))) +
  geom_jitter(aes(col=factor(cell.composition)), alpha = 0.5) 

ggplot(data = new_df, mapping = aes(x = mindist.tumour)) +
  geom_histogram(aes(fill=factor(cell.composition)), position = "dodge") +
  scale_x_continuous(breaks = seq(0,1000,100)) +
  facet_wrap(~cell.composition)

plot.dist.to.tumour <- new_df %>%
  #filter(cell.composition != "Pure_Tumour") %>%
  ggplot(mapping = aes(x = mindist.tumour, 
                       #after_stat(count)
                       )) +
  geom_density(aes(fill = cell.composition)) +
  scale_x_continuous(breaks = seq(0,1000,100)) +
  facet_grid(rows = vars(cell.composition),
             scales = "free_y") +
  ggtitle("Distances to pure tumour by cell type")
plot.dist.to.tumour

stats.new_df <- new_df %>%
  group_by(cell.composition) %>%
  summarise(median = median(mindist.tumour),
            min = min(mindist.tumour),
            max = max(mindist.tumour))
stats.new_df

# Add metadata distances
meta.toadd <- metadata2 %>%
  select(mindist.tumour:dist)
seuratobj.distances <- AddMetaData(seuratobj.aligned, metadata = meta.toadd)
head(seuratobj.distances@meta.data)


# Plot tumour distances
surrounding.tumour <- lapply(
  customSpatialFeaturePlot(seuratobj.distances, features = "mindist.tumour", return.list = TRUE),
  FUN = function(x) {
    x + scale_fill_gradientn(colours = colorscale(n = 100, rev = T)) 
  }) |>
  wrap_plots(nrow = 1, guides = "collect")
surrounding.tumour

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/distances/"), recursive = TRUE)
ggsave(filename = "distances_to_puretumour_by_celltype.png",
       plot = plot.dist.to.tumour,
       path = paste0(out.dir,"/plots/distances/"))
ggsave(filename = "distances_to_puretumour_spatial.png",
       plot = surrounding.tumour,
       path = paste0(out.dir,"/plots/distances/"))

all.plots <- list(plot.dist.to.tumour, surrounding.tumour)
save(all.plots, file = paste0("./results/ggplots/distances.RData"))

# Save data
saveRDS(seuratobj.distances, file = "./results/analysis/seuratobj.distances.rds")

