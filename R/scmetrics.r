#' scmetrics
#'
#' An R package for measuring the batch correction in single cell data
#'
#' @import methods
#' @import Matrix
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges rowData colData assays rbind
#' @importFrom stats kmeans
#' @importFrom Seurat as.Seurat as.SingleCellExperiment FindClusters FindNeighbors Idents
#' @importFrom BiocParallel bplapply
#' @importFrom aricode NMI ARI
#' @importFrom cluster silhouette
#' @importFrom igraph graph.adjacency components
#' @importFrom kBET kBET
#' @docType package
#' @name scco
#'
NULL
# > NULL
