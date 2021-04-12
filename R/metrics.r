#' metrics
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsNMI object
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsNMI'
	),
	function(
		x,
		method,
		label = 'cell_type',
		cluster = 'cluster'
	){

		stopifnot(!is.null(colData(x)[[label]]))
		stopifnot(!is.null(colData(x)[[cluster]]))
		NMI(colData(x)[[label]],  colData(x)[[cluster]])
	}
)

#' metrics
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsARI object
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsARI'
	),
	function(
		x,
		method,
		label = 'cell_type',
		cluster = 'cluster'
	){

		stopifnot(!is.null(colData(x)[[label]]))
		stopifnot(!is.null(colData(x)[[cluster]]))
		ARI(colData(x)[[label]],  colData(x)[[cluster]])
	}
)

#' metrics
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsARI object
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsASW'
	),
	function(
		x,
		method,
		label = 'cell_type',
		reduction = 'umap'
	){

		stopifnot(!is.null(colData(x)[[label]]))
		stopifnot(reduction %in% reducedDimNames(x))

		si <- silhouette(as.integer(factor(colData(x)[[label]])), dist(reducedDim(x, reduction)))
		summary(si)$avg.width
	}
)


#' metrics
#'
#' Graph connectivity 
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsARI object
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsGraphConnectivity'
	),
	function(
		x,
		method,
		label = 'cell_type',
		graph = 'snn'
	){

		stopifnot(!is.null(colData(x)[[label]]))
		stopifnot(is(metadata(x)[[graph]], 'dgCMatrix'))

		groups <- as.numeric(factor(colData(x)[[label]]))
		group_size <- table(groups)
		n_groups <- max(groups)

		lcc <- sapply(seq_len(n_groups), function(i){
			j <- groups == i
			com <- metadata(x)[[graph]][j, j] %>% 
				graph.adjacency() %>%
				components()
			max(com$csize)
		})

		sum(lcc / group_size) / n_groups

	}
)

#' metrics
#'
#' kBET
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsKBET object
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsKBET'
	),
	function(
		x,
		method,
		label = 'cell_type',
		reduction = 'umap'
	){
	}
)
