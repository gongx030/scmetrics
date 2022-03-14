#' measure
#'
#' @param x a Seurat object
#' @param method a MetricsNMICell object
#' @param label cell labels
#' @param cluster cell cluster
#' @importFrom aricode NMI
#' @export
#' 
setMethod(
	'measure',
	signature(
		x = 'Seurat',
		method = 'MetricsNMICell'
	),
	function(
		x,
		method,
		label = 'cell_type',
		cluster = 'cluster'
	){

		stopifnot(!is.null(x@meta.data[[label]]))
		stopifnot(!is.null(x@meta.data[[cluster]]))
		NMI(x@meta.data[[label]], x@meta.data[[cluster]])
	}
)


#' measure
#'
#' @param x a Seurat object
#' @param method a MetricsARICell object
#' @param label cell labels
#' @param cluster cell cluster
#' @importFrom aricode ARI
#' @export
#' 
setMethod(
	'measure',
	signature(
		x = 'Seurat',
		method = 'MetricsARICell'
	),
	function(
		x,
		method,
		label = 'cell_type',
		cluster = 'cluster'
	){

		stopifnot(!is.null(x@meta.data[[label]]))
		stopifnot(!is.null(x@meta.data[[cluster]]))
		ARI(x@meta.data[[label]], x@meta.data[[cluster]])
	}
)

#' measure
#'
#' @param x a Seurat object
#' @param method a MetricsARIBatch object
#' @param label cell labels
#' @param cluster cell cluster
#' @importFrom aricode ARI
#' @export
#' 
setMethod(
	'measure',
	signature(
		x = 'Seurat',
		method = 'MetricsARIBatch'
	),
	function(
		x,
		method,
		label = 'batch',
		cluster = 'cluster'
	){

		stopifnot(!is.null(x@meta.data[[label]]))
		stopifnot(!is.null(x@meta.data[[cluster]]))
		1 - ARI(x@meta.data[[label]], x@meta.data[[cluster]])
	}
)

#' measure
#'
#' @param x a Seurat object
#' @param method a MetricsASWCell object
#' @param label cell labels
#' @param reduction reduction name
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @export
#' 
setMethod(
	'measure',
	signature(
		x = 'Seurat',
		method = 'MetricsASWCell'
	),
	function(
		x,
		method,
		label = 'cell_type',
		reduction = 'pca'
	){

		stopifnot(!is.null(x@meta.data[[label]]))
		stopifnot(!is.null(x@reductions[[reduction]]))

		si <- silhouette(as.integer(factor(x@meta.data[[label]])), dist(x@reductions[[reduction]]@cell.embeddings))
		(summary(si)$avg.width + 1) / 2
	}
)


#' measure
#'
#' @param x a Seurat object
#' @param method a MetricsASWBatch object
#' @param label cell labels
#' @param reduction reduction name
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @export
#' 
setMethod(
	'measure',
	signature(
		x = 'Seurat',
		method = 'MetricsASWBatch'
	),
	function(
		x,
		method,
		label = 'batch',
		reduction = 'pca'
	){

		stopifnot(!is.null(x@meta.data[[label]]))
		stopifnot(!is.null(x@reductions[[reduction]]))

		si <- silhouette(as.integer(factor(x@meta.data[[label]])), dist(x@reductions[[reduction]]@cell.embeddings))
		1 - (summary(si)$avg.width + 1) / 2
	}
)



#' measure 
#'
#' @param x a Seurat object
#' @param method a MetricsGraphConnectivityBatch object
#' @param label cell labels
#' @param graph graph name
#' @importFrom igraph graph.adjacency components
#' @export
#' 
setMethod(
	'measure',
	signature(
		x = 'Seurat',
		method = 'MetricsGraphConnectivityBatch'
	),
	function(
		x,
		method,
		label = 'batch',
		graph = 'snn'
	){

		stopifnot(!is.null(x@meta.data[[label]]))
		stopifnot(!is.null(x@graphs[[graph]]))

		groups <- as.numeric(factor(x@meta.data[[label]]))
		group_size <- table(groups)
		n_groups <- max(groups)

		lcc <- sapply(seq_len(n_groups), function(i){
			j <- groups == i
			com <- x@graphs[[graph]][j, j] %>%
				graph.adjacency() %>%
				components()
			max(com$csize)
		})

		1 - sum(lcc / group_size) / n_groups

	}
)
