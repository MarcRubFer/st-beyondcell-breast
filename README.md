# st-beyondcell-breast
Author: Marcos Rubio-Fernández

Final Project Master

## Introduction

This repository contains the code for a spatial transcriptomic analysis in two breast cancer dataset. After processing and characterization, Beyondcell analysis will be carry out to determine therapeutic clusters.

## Breast cancer datasets

Breast cancer datasets were obtained from 10X genomics examples website (registration needed):

1. https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0
2. https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-2-1-standard-1-1-0

## Packages and versions

It is recommended to use conda/mamba enviroments located in envs directory

- R: version 4.3.1
- RStudio: version 2022.07.1 Build 554
- Seurat: version 4.3.0.1
- spacexr: version 2.2.1
- clustree: version 0.5.1
- scSHC: version 0.1.0
- ggsankey: version 0.0.99999
- PASTE: version
- fgsea: version 1.26.0
- ComplexHeatmap: 2.16.0
- Beyondcell: version 2.2.0

Some packages requiere devtools/remote installation:

To install spacexr package:

```
mamba install conda-forge::r-devtools
mamba install conda-forge::r-rfast
mamba install conda-forge::r-rcppgsl
mamba install conda-forge::r-rcppziggurat
```

In R terminal

```
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```