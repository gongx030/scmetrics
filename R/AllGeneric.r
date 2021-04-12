#' metrics
#'
#' The generic function of validate
#'
#' @param x a data object
#' @param method a metric method object
#' @param ... Other arguments
#'
setGeneric('metrics', function(x, method, ...) standardGeneric('metrics'))

#' validate
#'
#' The generic function of validate
#'
#' @param x a data object
#' @param ... Other arguments
#'
setGeneric('validate', function(x, ...) standardGeneric('validate'))


#' set_params
#'
#' The generic function of set_params
#'
#' @param method Method name
#' @param ... Other arguments
#'
setGeneric('set_params', function(method, ...) standardGeneric('set_params'))

#' cluster_cells
#'
#' The generic function of cluster_cells
#'
#' @param x a data object
#' @param method a metric method object
#' @param ... Other arguments
#'
setGeneric('cluster_cells', function(x, method, ...) standardGeneric('cluster_cells'))

