# conda activate beyondcell
rm(list = ls())

library("Seurat")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Load data matrix
mtx <- Read10X(data.dir = "./data/section1/raw_feature_bc_matrix/", gene.column = 2)

# Create Seurat object with data matrix
seuratobj <- CreateSeuratObject(counts = mtx, 
                                assay = "Spatial", 
                                project = "Section1")

# Load image to seuratobj, by default spatial image is not loaded.
image <- Read10X_Image(image.dir = "./data/section1/spatial/",
                       image.name = "tissue_lowres_image.png")
image@key <- paste0("Section1", "_")
# Correct the image data to match the Seurat object
image <- image[Cells(seuratobj)]
DefaultAssay(image) <- DefaultAssay(seuratobj)
# Add image to seuratobject
seuratobj@images[["Section1"]] <- image

# Save seuratobj as R object
dir.create(path = paste0(out.dir,"/analysis"), recursive = TRUE)
saveRDS(seuratobj, file = "./results/analysis/seuratobj.raw.rds")
