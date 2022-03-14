#' measure
#'
#' The generic function of measure
#'
#' @param x a data object
#' @param method a metric method object
#' @param ... Other arguments
#'
setGeneric('measure', function(x, method, ...) standardGeneric('measure'))

#' cluster_cells
#'
#' The generic function of cluster_cells
#'
#' @param x a data object
#' @param method a clustering method object
#' @param metrics a metrics method object
#' @param ... Other arguments
#'
setGeneric('cluster_cells', function(x, method, metrics, ...) standardGeneric('cluster_cells'))

