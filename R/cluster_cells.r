#' cluster_cells
#'
#' @param x a SingleCellExperiment object
#' @param method a Method parameter
#' @param metrics a Metrics parameter
#' @param label cell labels
#' @param reduction The reudction used to find neighboring cells.
#' @param ... Additional arguments
#' @export
#' @importFrom methods is
#' @importFrom Seurat FindNeighbors FindClusters Idents
#' @importFrom aricode NMI ARI
#' 
setMethod(
	'cluster_cells',
	signature(
		x = 'Seurat',
		method = 'ClusterLouvain',
		metrics = 'MetricsMethod'
	),
	function(
		x,
		method,
		metrics,
		label = 'cell_type',
		reduction = 'umap',
		...
	){

		stopifnot(is(x[[reduction]], 'DimReduc'))
		stopifnot(!is.null(x@meta.data[[label]]))

		k <- ncol(x[[reduction]])

		x <- x %>%
			FindNeighbors(reduction = reduction, k.param = method@n_neighbors, dims = 1:k, verbose = FALSE, graph.name = 'graph')

		res <- 10^seq(log10(method@resolution_min), log10(method@resolution_max), length.out = method@resolution_length)
		
		sprintf('cluster_cells | %s | FindClusters | # resolutions=%d', class(method), length(res)) %>% message()

		w <- NULL
		cls <- list()
		for (i in 1:length(res)){
			# Returns a Seurat object where the idents have been updated with
			# new cluster info; latest clustering results will be stored in
		  # object metadata under 'seurat_clusters'. Note that
			# 'seurat_clusters' will be overwritten everytime FindClusters is run
			x <- FindClusters(x, algorithm = method@algorithm, resolution = res[i], verbose = FALSE, graph.name = 'graph') 
			cls[[i]] <- Idents(x)
			w <- c(w, measure(x, metrics, label = label, cluster = 'seurat_clusters'))
		}

		h <- which.max(w)
		for (i in 1:length(w)){
			sprintf('cluster_cells | %s=%.3f | resolution=%.3e%s', metrics@name, w[i], res[i], ifelse(i == h, '*', '')) %>% message()
		}

		x@meta.data[[method@name]] <- cls[[h]]
		x
	}
)

#' cluster_cells
#'
#' @param x a SingleCellExperiment object
#' @param method a Method parameter
#' @param metrics a Metrics parameter
#' @param ... Additional arguments to cluster_cells('Seurat', 'ANY' 'ANY', ...)
#' @importFrom Seurat as.Seurat
#' @export
#' 
setMethod(
	'cluster_cells',
	signature(
		x = 'SingleCellExperiment',
		method = 'ANY',
		metrics = 'ANY'
	),
	function(
		x,
		method, 
		metrics,
		...
	){
		x <- as.Seurat(x)
		cluster_cells(x, method, metrics, ...)
	}
)
