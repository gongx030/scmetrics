#' set_params
#'
#' @param method The method name
#' @param ... Other arguments
#' 
setMethod(
	'set_params',
	signature(
		method  = 'character'
	),
	function(
		method,
		...
	){

		if (method == 'NMI'){
			new('MetricsNMI', ...)
		}else
			stop(sprintf('unknown method: %s', method))
	}
)
