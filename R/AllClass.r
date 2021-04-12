setClassUnion('listOrNULL', members = c('list', 'NULL'))
setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))

#' MetricsMethod
#'
setClass(
	'MetricsMethod'
)


#' MetricsNMI
#'
setClass(
	'MetricsNMI', 
	contains = 'MetricsMethod'
)

#' MetricsARI
#'
setClass(
	'MetricsARI', 
	contains = 'MetricsMethod'
)

#' MetricsASW
#'
setClass(
	'MetricsASW', 
	contains = 'MetricsMethod'
)

#' MetricsGraphConnectivity
#'
setClass(
	'MetricsGraphConnectivity', 
	contains = 'MetricsMethod'
)

#' MetricsKBET
#'
setClass(
	'MetricsKBET',
	slot = c(
		k0 = 'integer'
	),
	prototype = list(
		k0 = 50L
	),
	contains = 'MetricsMethod'
)

#' ClusterMethod
#'
setClass(
	'ClusterMethod'
)

#' ClusterLouvain
#'
setClass(
	'ClusterLouvain', 
	contains = 'ClusterMethod',
	slot = c(
		n_neighbors = 'integer',
		resolution_min = 'numeric',
		resolution_max = 'numeric',
		resolution_length = 'integer'
	),
	prototype = list(
		n_neighbors = 15L,
		resolution_min =  1e-2,
		resolution_max =  1e-0,
		resolution_length = 20L
	)
)
