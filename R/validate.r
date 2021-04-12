#' validate
#'
#' @param a SingleCellExperiment object
#'
setMethod(
	'validate',
	signature(
		x = 'SingleCellExperiment'
	),
	function(
		x,
		batch_field = 'batch',
		label_field = 'cell_type',
		reduction = 'umap',
		...
	){

		stopifnot(!is.null(colData(x)[[batch_field]]))
		stopifnot(!is.null(colData(x)[[label_field]]))
		stopifnot(reduction %in% reducedDimNames(x))
		stopifnot(!is.null(assays(x)$counts))
		stopifnot(class(assays(x)$counts) == 'dgCMatrix')
	}
)
