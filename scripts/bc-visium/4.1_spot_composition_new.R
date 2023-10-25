rm(list = ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("ggsci")
library("RColorBrewer")
library("ggplotify")
library("grid")
library("cowplot")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

############
## CODE ##
# Random seed
set.seed(1)

# Load seuratobject
seuratobj.deconvoluted <- readRDS(file = "./results/analysis/seuratobj.deconvoluted.rds")
head(seuratobj.deconvoluted@meta.data)


# Subset cell types proportions metadata
cell.proportions <- seuratobj.deconvoluted@meta.data %>%
  select(B.cells:T.cells)
head(cell.proportions)

# Annotate spots as "Tumour" if Cancer Epitheli al prop >= 0.65
cell.proportions <- cell.proportions %>%
  mutate(spot.composition = case_when(Cancer.Epithelial >= 0.65 ~ "TUMOUR",
                                      TRUE ~ "Others"))


ggplot(cell.proportions, aes(x=spot.composition)) +
  geom_bar()

library(dplyr)

df <- cell.proportions %>%
  mutate(
    spot.composition = case_when(
      Cancer.Epithelial >= 0.65 ~ "TUMOUR",
      #Cancer.Epithelial > 0.8 ~ "PURE.TUMOUR",
      #Cancer.Epithelial >= 0.65 & Cancer.Epithelial <= 0.8 ~ "TUMOUR",
      TRUE ~ toupper(names(select(., -Cancer.Epithelial))[max.col(select(., -Cancer.Epithelial), "first")])
      #TRUE ~ names(.)[max.col(., "first")]
      ),
    spot.composition = case_when(
      spot.composition == "B.CELLS" | spot.composition =="T.CELLS" ~ "LYMPHOID",
      spot.composition == "NORMAL.EPITHELIAL" | spot.composition == "PVL" ~ "OTHERS",
      TRUE ~ spot.composition
    )
  )

head(df)
ggplot(df, aes(x=spot.composition)) +
  geom_bar()

objeto.prueba <- seuratobj.deconvoluted
objeto.prueba <- AddMetaData(objeto.prueba, metadata = df)
head(objeto.prueba)
cells <- WhichCells(objeto.prueba, expression = spot.composition == "TUMOUR")
lymphoid <- WhichCells(objeto.prueba, expression = spot.composition == "LYMPHOID")
list.cells <- list(cells,lymphoid)
SpatialDimPlot(objeto.prueba, group.by = "spot.composition", cells.highlight = list.cells, cols.highlight = c("red","blue"))

table(df$spot.composition)

df.prueba <- cell.proportions %>%
  filter(spot.composition != "TUMOUR") %>%
  select(-c(Cancer.Epithelial,spot.composition)) 
df.prueba <- df.prueba %>%
  mutate(Max_Var = names(.)[max.col(., "first")]) %>%
  select(Max_Var)

head(df.prueba)
df.prueba

merged.df <- merge(x=cell.proportions, y=df.prueba, by = "row.names")
head(merged.df)
head(cell.proportions)


right_join(x= cell.proportions, y=df.prueba$Max_Var)








