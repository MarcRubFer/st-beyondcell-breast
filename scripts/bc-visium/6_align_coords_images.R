rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)


## FUNCTIONS ##
# Obtain mode (the value that occurs most often)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object 
seuratobj.clusters <- readRDS(file = "./results/analysis/seuratobj.clusters.rds")

# Read and plot non-aligned coordinates
## Raw coords
cln <- c("spots","row","col")
### Slide 1
raw.coords.slide1 <- read_tsv(file = "./data/tsv/raw-coords-slide1.tsv", 
                              skip = 1, 
                              col_names = cln) %>%
  mutate(z.coord = 1,
         row = round(row, digits = 0),
         col = round(col, digits = 0))
rownames(raw.coords.slide1) <- raw.coords.slide1$spots

### Slide 2
raw.coords.slide2 <- read_tsv(file = "./data/tsv/raw-coords-slide2.tsv", 
                              skip = 1, 
                              col_names = cln) %>%
  mutate(z.coord = 2,
         row = round(row, digits = 0),
         col = round(col, digits = 0))
rownames(raw.coords.slide2) <- raw.coords.slide2$spots

### Merge Slide 1-2
raw.coords.merged <- rbind(raw.coords.slide1, raw.coords.slide2) 

### Plot merged coords
map.raw.coords <- ggplot(raw.coords.merged, 
                         aes(x=row, y=col, 
                             color = factor(z.coord))) +
  geom_point(aes(shape =factor(z.coord)),
             alpha = 0.9) +
  scale_color_manual(values = c("1" = "red",
                                "2" = "black")) +
  scale_shape_manual(values = c(0, 4)) +
  ggtitle("Raw coordinates")

## Scaled coords
cln <- c("spots","imagerow","imagecol")
### Slide 1
scaled.coords.slide1 <- read_tsv(file = "./data/tsv/coords-slide1.tsv", 
                                 skip = 1, 
                                 col_names = cln) %>%
  mutate(z.coord = 1,
         imagerow = round(imagerow, digits = 0),
         imagecol = round(imagecol, digits = 0))

rownames(scaled.coords.slide1) <- scaled.coords.slide1$spots

### Slide 2
scaled.coords.slide2 <- read_tsv(file = "./data/tsv/coords-slide2.tsv", 
                                 skip = 1, 
                                 col_names = cln) %>%
  mutate(z.coord = 2,
         imagerow = round(imagerow, digits = 0),
         imagecol = round(imagecol, digits = 0))

rownames(scaled.coords.slide2) <- scaled.coords.slide2$spots

### Merge slide 1-2
scaled.coords.merged <- rbind(scaled.coords.slide1, scaled.coords.slide2)

### Plot merged coords
map.scaled.coords <- ggplot(scaled.coords.merged, aes(x=imagerow, y=imagecol, color = factor(z.coord))) +
  geom_point(aes(shape =factor(z.coord)),alpha = 0.9) +
  scale_color_manual(values = c("1" = "red",
                                "2" = "black")) +
  scale_shape_manual(values = c(0, 4)) +
  ggtitle("Scaled coordinates")

###############################################################################

# Note: To proceed to the next step, it is necessary to run the "paste.py" script beforehand -
# This script process raw/aligned coordinates and generates aligned coordinates

###############################################################################

#Read and plot aligned coordinates
## Aligned raw coords
### Slide 1
aligned.raw.coords.slide1 <- read_tsv(file = "./data/tsv/aligned_raw_coords_sc1.tsv") %>%
  rename(spots = ...1) %>%
  mutate(z.coord = 1,
         adj_x = round(adj_x, digits = 0),
         adj_y = round(adj_y, digits = 0))
rownames(aligned.raw.coords.slide1) <- aligned.raw.coords.slide1$spots

### Slide 2
aligned.raw.coords.slide2 <- read_tsv(file = "./data/tsv/aligned_raw_coords_sc2.tsv") %>%
  rename(spots = ...1) %>%
  mutate(z.coord = 2,
         adj_x = round(adj_x, digits = 0),
         adj_y = round(adj_y, digits = 0)) 
rownames(aligned.raw.coords.slide2) <- aligned.raw.coords.slide2$spots

### Merge slides 1-2
aligned.raw.merged <- rbind(aligned.raw.coords.slide1,aligned.raw.coords.slide2)

### Plot merged coords
map.aligned.raw.coords <- ggplot(aligned.raw.merged, 
                                 aes(x=adj_x, y=adj_y, 
                                     color = factor(z.coord))) +
  geom_point(aes(shape =factor(z.coord)),alpha = 0.9) +
  scale_color_manual(values = c("1" = "red",
                                "2" = "black")) +
  scale_shape_manual(values = c(0, 4)) +
  ggtitle("Aligned raw coordinates")

## Aligned scaled coords
### Slide 1
aligned.scaled.coords.slide1 <- read_tsv(file = "./data/tsv/aligned_coords_sc1.tsv") %>%
  rename(spots = ...1) %>%
  mutate(z.coord = 1,
         adj_x = round(adj_x, digits = 0),
         adj_y = round(adj_y, digits = 0))
rownames(aligned.scaled.coords.slide1) <- aligned.scaled.coords.slide1$spots

### Slide 2
aligned.scaled.coords.slide2 <- read_tsv(file = "./data/tsv/aligned_coords_sc2.tsv") %>%
  rename(spots = ...1) %>%
  mutate(z.coord = 2,
         adj_x = round(adj_x, digits = 0),
         adj_y = round(adj_y, digits = 0)) 
rownames(aligned.scaled.coords.slide2) <- aligned.scaled.coords.slide2$spots

### Merge slides 1-2
aligned.scaled.merged <- rbind(aligned.scaled.coords.slide1,aligned.scaled.coords.slide2)

### Plot merged coordinates
map.aligned.scaled.coords <- ggplot(aligned.scaled.merged, aes(x=adj_x, y=adj_y, color = factor(z.coord))) +
  geom_point(aes(shape =factor(z.coord)),alpha = 0.9) +
  scale_color_manual(values = c("1" = "red",
                                "2" = "black")) +
  scale_shape_manual(values = c(0, 4)) +
  ggtitle("Aligned scaled coordinates")

### Pachtwork four plots

patch.coords <- (map.raw.coords | map.aligned.raw.coords) / 
  (map.scaled.coords | map.aligned.scaled.coords) 

####################
# Note: Alignment is similar between raw and scaled but raw is more easy to explain
# For this, aligned raw coords are selected 
###################

## Calculate the step between coords in x-axis and y-axis

step.raw.x <- aligned.raw.merged %>%
  mutate(spots = NULL) %>%
  group_by(adj_y, z.coord) %>%
  arrange(z.coord,adj_x, adj_y) %>%
  mutate(diff = c((sort(adj_x)[-1]), NA),
         step = abs(diff-adj_x))

step.raw.y <- aligned.raw.merged %>%
  mutate(spots = NULL) %>%
  group_by(adj_x, z.coord) %>%
  arrange(z.coord,adj_x, adj_y) %>%
  mutate(diff = c((sort(adj_y)[-1]), NA),
         step = abs(diff-adj_y))

step.x <- getmode(step.raw.x$step)
step.y <- getmode(step.raw.y$step)

## Calculate coords in the same scale: distance between spots is 100um and
## between slides 10um. To establish the same units we need to correct the 
## values of these coordinates.
corrected.raw.coords <- aligned.raw.merged %>%
  mutate(corr_x = (adj_x + abs(min(adj_x)) + step.x) * 100,
         corr_y = (adj_y + abs(min(adj_y)) + step.y) * 100,
         corr_z = z.coord * 10)

rnm <- corrected.raw.coords$spots
corrected.raw.coords$spots <- NULL
rownames(corrected.raw.coords) <- rnm

## Plot corrected coords
map.corrected.coords <- ggplot(corrected.raw.coords, aes(x=corr_x, y=corr_y, color = factor(z.coord))) +
  geom_point(aes(shape =factor(z.coord)),alpha = 0.9) +
  scale_color_manual(values = c("1" = "red",
                                "2" = "black")) +
  scale_shape_manual(values = c(0, 4)) +
  ggtitle("Corrected raw coordinates")

# Add metadata (aligned.raw.coords and corrected.coords) to seurat object
seuratobj.aligned <- AddMetaData(seuratobj.clusters, metadata = corrected.raw.coords)
seuratobj.aligned@meta.data

# Save plots and ggplots
dir.create(path = paste0(out.dir,"/plots/alignment/"), recursive = TRUE)
ggsave(filename = "raw_coordinates.png",
       plot = map.raw.coords,
       path = paste0(out.dir,"/plots/alignment/"))
ggsave(filename = "scaled_coordinates.png",
       plot = map.scaled.coords,
       path = paste0(out.dir,"/plots/alignment/"))
ggsave(filename = "aligned_raw_coordinates.png",
       plot = map.aligned.raw.coords,
       path = paste0(out.dir,"/plots/alignment/"))
ggsave(filename = "aligned_scaled_coordinates.png",
       plot = map.aligned.scaled.coords,
       path = paste0(out.dir,"/plots/alignment/"))
ggsave(filename = "patchwork_coordinates.png",
       plot = patch.coords,
       path = paste0(out.dir,"/plots/alignment/"))
ggsave(filename = "corrected_coordinates.png",
       plot = map.corrected.coords,
       path = paste0(out.dir,"/plots/alignment/"))

dir.create(paste0(out.dir, "/ggplots/"), recursive = TRUE, showWarnings = FALSE)
all.plots <- list(map.raw.coords,map.aligned.raw.coords,map.aligned.raw.coords,
                  map.aligned.scaled.coords,patch.coords,map.corrected.coords)
save(all.plots, file = paste0("./results/ggplots/alignment.RData"))

# Save object
saveRDS(seuratobj.aligned, file = "./results/analysis/seuratobj.aligned.rds")
