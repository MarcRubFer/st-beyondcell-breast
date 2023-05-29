# conda activate beyondcell
rm(list = ls())

library("Seurat")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# Load data matrix
mtx.ref <- Read10X(data.dir = "./data/reference/BrCa_Atlas_Count_out/")

# Read metadata
meta.ref <- read.csv(file = "./data/reference/Whole_miniatlas_meta.csv", 
                     row.names = 1)
meta.ref <- meta.ref[-1, ] %>%
  mutate(nCount_RNA = as.numeric(nCount_RNA),
         nFeature_RNA = as.numeric(nFeature_RNA),
         Percent_mito = as.numeric(Percent_mito))

# Subset matrix and metadata
intersect.cells <- intersect(colnames(mtx.ref), rownames(meta.ref))
mtx.ref <- mtx.ref[, intersect.cells]
meta.ref <- meta.ref[intersect.cells, ]

# Create Seurat object with data matrix
seuratobj.reference <- CreateSeuratObject(counts = mtx.ref, 
                                          assay = "RNA", 
                                          meta.data = meta.ref,
                                          project = "Single-cell")



# Save seuratobj as R object
dir.create(path = paste0(out.dir,"/analysis"), recursive = TRUE)
saveRDS(seuratobj.reference, file = "./results/analysis/seuratobj.reference.rds")
