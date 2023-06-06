# conda activate beyondcell
rm(list = ls())

library("Seurat")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Load data matrix from section 1 and section 2
mtx <- Read10X(data.dir = "./data/section1/raw_feature_bc_matrix/", gene.column = 2)
mtx.2 <- Read10X(data.dir = "./data/section2/raw_feature_bc_matrix/", gene.column = 2)

# Create Seurat objects with data matrixes
seuratobj.s1 <- CreateSeuratObject(counts = mtx, 
                                assay = "Spatial", 
                                project = "Section1")

seuratobj.s2 <- CreateSeuratObject(counts = mtx.2,
                                   assay = "Spatial",
                                   project = "Section2")

# Load image to seuratobj, by default spatial image is not loaded.
image.s1 <- Read10X_Image(image.dir = "./data/section1/spatial/",
                       image.name = "tissue_lowres_image.png")
image.s1@key <- paste0("Section1", "_")

image.s2 <- Read10X_Image(image.dir = "./data/section2/spatial/",
                          image.name = "tissue_lowres_image.png")
image.s2@key <- paste0("Section2", "_")

# Correct the image data to match the Seurat object
image.s1 <- image.s1[Cells(seuratobj.s1)]
DefaultAssay(image.s1) <- DefaultAssay(seuratobj.s1)

image.s2 <- image.s2[Cells(seuratobj.s2)]
DefaultAssay(image.s2) <- DefaultAssay(seuratobj.s2)

# Add image to seuratobject
seuratobj.s1@images[["Section1"]] <- image.s1

seuratobj.s2@images[["Section2"]] <- image.s2

# Merge Seurat objects
seuratobj.merged <- merge(seuratobj.s1, y = seuratobj.s2, add.cell.ids = c("Sc1", "Sc2"), project = "Breast")

# Save seuratobjs as R object
dir.create(path = paste0(out.dir,"/analysis/section1"), recursive = TRUE)
saveRDS(seuratobj.s1, file = paste0(out.dir,"/analysis/section1/seuratobj.s1.raw.rds"))

dir.create(path = paste0(out.dir,"/analysis/section2"), recursive = TRUE)
saveRDS(seuratobj.s2, file = paste0(out.dir,"/analysis/section2/seuratobj.s2.raw.rds"))

dir.create(path = paste0(out.dir,"/analysis/merged"), recursive = TRUE)
saveRDS(seuratobj.merged, file = paste0(out.dir,"/analysis/merged/seuratobj.merged.raw.rds"))
