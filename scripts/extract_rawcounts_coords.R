rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")

# Load seuratobject
seuratobj.clusters <- readRDS(file = "./results/analysis/seuratobj.clusters.rds")

# Extract raw coordinates from each section
raw.coords.slide1 <- seuratobj.clusters@images$Section1@coordinates[,c("row","col")]
raw.coords.slide2 <- seuratobj.clusters@images$Section2@coordinates[,c("row","col")]

# Extract scaled coordinates from each section
l <- lapply(X = c("Section1", "Section2"), function(section){
  coords <- GetTissueCoordinates(seuratobj.clusters, image = section)
  return(coords)
})

scaled.coords.slide1 <- l[[1]]
scaled.coords.slide2 <- l[[2]]

# Extract raw counts for each slide
raw.counts <- as.data.frame(GetAssayData(seuratobj.clusters, slot = "counts"))
raw.counts.sc1 <- raw.counts %>%
  select(starts_with("Sc1_"))
raw.counts.sc2 <- raw.counts %>%
  select(starts_with("Sc2_"))

# Save dataframes as .tsv
write.table(raw.coords.slide1, file = "./data/tsv/raw-coords-slide1.tsv", sep = "\t")
write.table(raw.coords.slide2, file = "./data/tsv/raw-coords-slide2.tsv", sep = "\t")
write.table(scaled.coords.slide1, file = "./data/tsv/coords-slide1.tsv", sep = "\t")
write.table(scaled.coords.slide2, file = "./data/tsv/coords-slide2.tsv", sep = "\t")
write.table(raw.counts.sc1, file = "./data/tsv/raw-counts-sc1.tsv", sep = "\t")
write.table(raw.counts.sc2, file = "./data/tsv/raw-counts-sc2.tsv", sep = "\t")

