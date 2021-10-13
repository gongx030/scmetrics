#' metrics
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsNMI object
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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

#' metrics
#'
#' AUROC Score 
#'
#' @import irlba
#' @import SingleCellExperiment
#' @import coop
#' @import flavin
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsARI object
#' @export
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsAUROC'
	),
	function(
		x,
		method,
		label = 'cell_type',
		graph = 'snn'
	){

		stopifnot(!is.null(colData(x)[[label]]))
		stopifnot(is(metadata(x)[[graph]], 'dgCMatrix'))

    ##get ontology
    species = read_gpa("ontology/goa_human.gpa.gz", accession = "SYMBOL", 
                 database = org.Hs.eg.db, ontology = ontology, propagate = T)

    ##creating svd matrix
    data = assays(x)$counts
    data = data + 1 #Pseudocount for log normalization
    data = log(data)
    data = as.matrix(data)
    test = irlba::irlba(data)
    test = test$u
    rownames(test)=rownames(data)
    test = coop::cosine(t(test))

    ##running AUROC
    go = get(method@species)
		
		nrow = nrow(coexpr)
		genes = rownames(coexpr)
		net = matrix(rank(coexpr, na.last = "keep", ties.method = "average"), 
					 nrow = nrow, ncol = nrow)
		rownames(net) = colnames(net) = genes
		net = net / max(net, na.rm = T)
		diag(net) = 1
		
		go_subset = go %>%
		  dplyr::filter(SYMBOL %in% genes) %>%
		  dplyr::select(SYMBOL, GO.ID) %>%
		  as.matrix()
		  
		if (nrow(go_subset) == 0)
			stop("no GO terms annotated to network genes. check right GO file")
			
		go_terms = unique(go_subset[, "GO.ID"])
		annotations = make_annotations(go_subset, genes, go_terms)
		gba = run_GBA(net, annotations, min = 0, max = 1e4)
		
		aurocs = gba[[1]][, "auc"]
		all_terms = colSums(annotations)
		n_terms = all_terms[names(aurocs)]
		dataset = names(dataset_list[1])
		coef = metric_list[n]
		result = data.frame(dataset = dataset, coefficient = coef, 
							term = names(aurocs), auroc = aurocs, 
							n_proteins = n_terms, 
							pct_proteins = n_terms / length(genes))
		
		filtered = dplyr::filter(result, term %in% slim$id & 
							 dplyr::between(n_proteins, 10, 1000)) %>%
		dplyr::select(dataset, coefficient, term, auroc)
		
		mean_auc = mean(filtered$auroc)
    return(mean_auc)
	}
)

#' metrics
#'
#' PPI Z-Score 
#'
#' @import irlba
#' @import SingleCellExperiment
#' @import coop
#' @import flavin
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsARI object
#' @export
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsAUROC'
	),
	function(
		x,
		method,
		label = 'cell_type',
		graph = 'snn'
	){

		stopifnot(!is.null(colData(x)[[label]]))
		stopifnot(is(metadata(x)[[graph]], 'dgCMatrix'))

    ##get ontology
    species = read_gpa("ontology/goa_human.gpa.gz", accession = "SYMBOL", 
                 database = org.Hs.eg.db, ontology = ontology, propagate = T)

    ##creating svd matrix
    data = assays(x)$counts
    data = data + 1 #Pseudocount for log normalization
    data = log(data)
    data = as.matrix(data)
    test = irlba::irlba(data)
    test = test$u
    rownames(test)=rownames(data)
    test = coop::cosine(t(test))

    ##running PPI Z-Score
    		# read protein-coding genes
		coding = read.delim("protein_coding_genes.txt.gz")

		#file = test
		coexpr = test ## coexpr
		# replace missing values with median (ZI kendall)
		coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
		# replace infinite values with the minimum (binomial)
		coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)


		species = "human"

		# create ranked coexpression data frame
		tri = upper.tri(coexpr)
		idxs = which(tri, arr.ind = T)
		int = data.frame(id1 = rownames(coexpr)[idxs[,1]], 
						 id2 = rownames(coexpr)[idxs[,2]], 
						 correlation = coexpr[tri])
		int %<>% arrange(-correlation)


		# read original networks
		net_names = c("STRING")
		original_files = paste0("", net_names, "/", species, ".rds")
		nets = map(original_files, ~ graph_from_data_frame(drop_na(readRDS(.)))) %>%
		  setNames(net_names)
		  
		bootstraps = 100
		results = data.frame()
		rewire_dir = "./STRING"
		output_dir = "/scratch.global/skiex003/Etv2_Cancer/fig2/output/"
		dataset = names(dataset_list[k])
		coef = metric_list[n]
		output = file.path(output_dir, paste0(dataset, "_", "second_half", "_", as.character(seed), ".txt"))
		for (net_name in net_names) {
		  message("analyzing network ", net_name, " ...")
		  net = nets[[net_name]]
		  
		  # calculate observed overlap 
		  message("  calculating observed overlap with ", net_name, " ...")
		  observed = numeric(0)
		  cutoffs = c(2e4, 5e4, 1e5)
		  for (cutoff in cutoffs) {
			subset = graph_from_data_frame(int[seq_len(cutoff), ])
			obs = length(E(intersection(subset, net, keep.all.vertices = F)))
			observed %<>% c(obs)
		  }

		  # calculate overlap in rewired networks     
		  random = list()
		  for (bootstrap in seq_len(bootstraps)) {
			if (bootstrap %% 10 == 0)
			  message("  analyzing bootstrap ", bootstrap, " of ", bootstraps, " ...")
			
			# read rewired network
			rewire_file = file.path(rewire_dir, paste0(
			  net_name, "-", species, "-", bootstrap, ".rds"))
			rewire = graph_from_data_frame(drop_na(readRDS(rewire_file)))
			
			# calculate rewired network overlap at different cutoffs
			rnd = numeric(0)
			for (cutoff in cutoffs) {
			  subset = graph_from_data_frame(int[seq_len(cutoff), ])
			  rnd %<>% c(length(
				E(intersection(subset, rewire, keep.all.vertices = F))))
			}
			# record
			random[[bootstrap]] = rnd
		  }
		  
		  # calculate outcomes
		  rnd_mean = map_dbl(seq_len(length(cutoffs)), ~ 
							   mean(map_dbl(random, .), na.rm = T))
		  rnd_median = map_dbl(seq_len(length(cutoffs)), ~
								 median(map_dbl(random, .), na.rm = T))
		  rnd_sd = map_dbl(seq_len(length(cutoffs)), ~
							 sd(map_dbl(random, .), na.rm = T))
		  enrs = observed / rnd_median
		  z_scores = (observed - rnd_mean) / rnd_sd
	  }
    return(z_scores)
  }
)

#' metrics
#'
#' Reactome Z-Score 
#'
#' @import irlba
#' @import SingleCellExperiment
#' @import coop
#' @import flavin
#'
#' @param x a SingleCellExperiment object
#' @param method a MetricsARI object
#' @export
#' 
setMethod(
	'metrics',
	signature(
		x = 'SingleCellExperiment',
		method = 'MetricsReactomeScore'
	),
	function(
		x,
		method,
		label = 'cell_type',
		graph = 'snn'
	){

		stopifnot(!is.null(colData(x)[[label]]))
		stopifnot(is(metadata(x)[[graph]], 'dgCMatrix'))

    ##get ontology
    species = read_gpa("ontology/goa_human.gpa.gz", accession = "SYMBOL", 
                 database = org.Hs.eg.db, ontology = ontology, propagate = T)

    ##creating svd matrix
    data = assays(x)$counts
    data = data + 1 #Pseudocount for log normalization
    data = log(data)
    data = as.matrix(data)
    test = irlba::irlba(data)
    test = test$u
    rownames(test)=rownames(data)
    test = coop::cosine(t(test))

    ##running PPI Z-Score
    # read protein-coding genes
		coding = read.delim("data/protein_coding_genes.txt.gz")

		#file = test
		coexpr = test ## coexpr
		# replace missing values with median (ZI kendall)
		coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
		# replace infinite values with the minimum (binomial)
		coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)


		species = "human"

		# create ranked coexpression data frame
		tri = upper.tri(coexpr)
		idxs = which(tri, arr.ind = T)
		int = data.frame(id1 = rownames(coexpr)[idxs[,1]], 
						 id2 = rownames(coexpr)[idxs[,2]], 
						 correlation = coexpr[tri])
		int %<>% arrange(-correlation)


		# read original networks
		net_names = c("reactome")
		original_files = paste0("", net_names, "/", species, ".rds")
		nets = map(original_files, ~ graph_from_data_frame(drop_na(readRDS(.)))) %>%
		  setNames(net_names)
		  
		bootstraps = 100
		results = data.frame()
		rewire_dir = "./reactome"
		output_dir = "/scratch.global/skiex003/Etv2_Cancer/fig2/output/"
		dataset = names(dataset_list[k])
		coef = metric_list[n]
		output = file.path(output_dir, paste0(dataset, "_", "second_half", "_", as.character(seed), ".txt"))
		for (net_name in net_names) {
		  message("analyzing network ", net_name, " ...")
		  net = nets[[net_name]]
		  
		  # calculate observed overlap 
		  message("  calculating observed overlap with ", net_name, " ...")
		  observed = numeric(0)
		  cutoffs = c(2e4, 5e4, 1e5)
		  for (cutoff in cutoffs) {
			subset = graph_from_data_frame(int[seq_len(cutoff), ])
			obs = length(E(intersection(subset, net, keep.all.vertices = F)))
			observed %<>% c(obs)
		  }

		  # calculate overlap in rewired networks     
		  random = list()
		  for (bootstrap in seq_len(bootstraps)) {
			if (bootstrap %% 10 == 0)
			  message("  analyzing bootstrap ", bootstrap, " of ", bootstraps, " ...")
			
			# read rewired network
			rewire_file = file.path(rewire_dir, paste0(
			  net_name, "-", species, "-", bootstrap, ".rds"))
			rewire = graph_from_data_frame(drop_na(readRDS(rewire_file)))
			
			# calculate rewired network overlap at different cutoffs
			rnd = numeric(0)
			for (cutoff in cutoffs) {
			  subset = graph_from_data_frame(int[seq_len(cutoff), ])
			  rnd %<>% c(length(
				E(intersection(subset, rewire, keep.all.vertices = F))))
			}
			# record
			random[[bootstrap]] = rnd
		  }
		  
		  # calculate outcomes
		  rnd_mean = map_dbl(seq_len(length(cutoffs)), ~ 
							   mean(map_dbl(random, .), na.rm = T))
		  rnd_median = map_dbl(seq_len(length(cutoffs)), ~
								 median(map_dbl(random, .), na.rm = T))
		  rnd_sd = map_dbl(seq_len(length(cutoffs)), ~
							 sd(map_dbl(random, .), na.rm = T))
		  enrs = observed / rnd_median
		  z_scores = (observed - rnd_mean) / rnd_sd
    }
    return(z_scores)
	}
)