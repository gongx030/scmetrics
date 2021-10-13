# scmetrics
scmetrics is a package that implements a common work flow for several single cell RNA sequence metrics

## Prerequisites
The following R packages are required for installation of scmetrics:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "AnnotationDbi"))
devtools::install_github('skinnider/flavin')
```

The current available methods on the development branch are:
  * ARI
  * NMI
  * ASW
  * Kbet
  * Graph Connectivity
  * AUROC Score
  * PPI Z-Score
  * Reactome Z-Score
  
