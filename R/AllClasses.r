setClassUnion('listOrNULL', members = c('list', 'NULL'))
setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))

#' MetricsMethod
#'
setClass(
	'MetricsMethod',
	slot = c(
		name = 'character'
	),
	prototype = list(
	)
)

#' MetricsNMI
#'
setClass(
	'MetricsNMICell', 
	prototype = list(
		name = 'NMI_cell'
	),
	contains = 'MetricsMethod',
)

#' MetricsARI
#'
setClass(
	'MetricsARICell', 
	prototype = list(
		name = 'ARI_cell'
	),
	contains = 'MetricsMethod'
)

#' MetricsARIBatch
#'
setClass(
	'MetricsARIBatch', 
	prototype = list(
		name = 'ARI_batch'
	),
	contains = 'MetricsMethod'
)

#' MetricsASWCell
#'
setClass(
	'MetricsASWCell', 
	prototype = list(
		name = 'ASW_cell'
	),
	contains = 'MetricsMethod'
)

#' MetricsASWBatch
#'
setClass(
	'MetricsASWBatch', 
	prototype = list(
		name = 'ASW_batch'
	),
	contains = 'MetricsMethod'
)

#' MetricsGraphConnectivityBatch
#'
setClass(
	'MetricsGraphConnectivityBatch', 
	prototype = list(
		name = 'GC_batch'
	),
	contains = 'MetricsMethod'
)

#' ClusterMethod
#'
setClass(
	'ClusterMethod',
	slot = c(
		name = 'character'
	)
)

#' ClusterLouvain
#'
setClass(
	'ClusterLouvain', 
	contains = 'ClusterMethod',
	slot = c(
		algorithm = 'integer',
		n_neighbors = 'integer',
		resolution_min = 'numeric',
		resolution_max = 'numeric',
		resolution_length = 'integer'
	),
	prototype = list(
		name = 'ClusterLouvain',
		algorithm = 3L,
		n_neighbors = 15L,
		resolution_min =  1e-2,
		resolution_max =  1e-0,
		resolution_length = 20L
	)
)

