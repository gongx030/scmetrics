#' cluster_cells
#'
#' @param x a SingleCellExperiment object
#' @param method a 
#' @param ... Arguments passed to validate(...)
#' @export
#' 
setMethod(
	'cluster_cells',
	signature(
		x = 'SingleCellExperiment',
		method = 'ClusterLouvain'
	),
	function(
		x,
		method,
		reduction = 'umap',
		label = 'cell',
		cluster = 'cluster',
		...
	){

		stopifnot(reduction %in% reducedDimNames(x))
		stopifnot(!is.null(colData(x)[[label]]))

		k <- ncol(reducedDim(x, reduction))

		if (!is.null(assays(x)$counts) && !is.null(assays(x)$logcounts)){
			y <- as.Seurat(x, counts = 'counts', data = 'logcounts')
		}else if (is.null(assays(x)$counts) && !is.null(assays(x)$logcounts)){
			y <- as.Seurat(x, counts = NULL, data = 'logcounts')
		}else if (!is.null(assays(x)$counts) && is.null(assays(x)$logcounts)){
			y <- as.Seurat(x, counts = 'counts', data = NULL)
		}else
			stop(sprintf('assays(x)$counts or assays(x)$logcounts cannot be NULL'))

		y <- y %>%
			FindNeighbors(reduction = reduction, k.param = method@n_neighbors, dims = 1:k, verbose = FALSE, graph.name = 'graph')

		metadata(x)$snn <- as(y@graphs[['graph']], 'dgCMatrix')

		res <- 10^seq(log10(params@resolution_min), log10(params@resolution_max), length.out = params@resolution_length)
		
		sprintf('cluster_cells | %s | FindClusters | # resolutions=%d', class(method), length(res)) %>% message()

		cls <- lapply(res, function(r) FindClusters(y, algorithm = 3, resolution = r, verbose = FALSE, graph.name = 'graph') %>% Idents())
		nmi <- sapply(1:length(cls), function(i) NMI(cls[[i]], colData(x)[[label]]))

		h <- which.max(nmi)
		for (i in 1:length(res)){
			sprintf('cluster_cells | nmi=%.3f | resolution=%.3e%s', nmi[i], res[i], ifelse(i == h, '*', '')) %>% message()
		}

		colData(x)[[cluster]] <- cls[[h]]
		x
	}
)
