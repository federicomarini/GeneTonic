# similarity_matrix <- km_macro

#' Compute fuzzy clusters of gene sets
#' 
#' Compute fuzzy clusters of different gene sets, aiming to identify grouped 
#' categories that can better represent the distinct biological themes in the
#' enrichment results
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be displayed.
#' @param similarity_matrix A similarity matrix between gene sets. Can be e.g. 
#' computed with [create_kappa_matrix()] or [create_jaccard_matrix()] or a similar
#' function, returning a symmetric matrix with numeric values (max = 1). If not 
#' provided, this will be computed on the fly with [create_kappa_matrix()]
#' @param similarity_threshold A numeric value for the similarity matrix, used to 
#' determine the initial seeds as in the implementation of DAVID. Higher values
#' will lead to more genesets being initially unclustered, leading to a  functional 
#' classification result with fewer groups and fewer geneset members. Defaults to 0.35,
#' recommended to not go below 0.3 (see DAVID help pages)
#' @param fuzzy_seeding_initial_neighbors Integer value, corresponding to the minimum
#' geneset number in a seeding group. Lower values will lead to the inclusion of more
#' genesets in the functional groups, and may generate a lot of small size groups. 
#' Defaults to 3
#' @param fuzzy_multilinkage_rule Numeric value, comprised between 0 and 1. This 
#' parameter will determine how the seeding groups merge with each other, by specifying
#' the percentage of shared genesets required to merge the two subsets into one
#' group. Higher values will give sharper separation between the groups of genesets.
#' Defaults to 0.5 (50%)
#' 
#' @references
#' See https://david.ncifcrf.gov/helps/functional_classification.html#clustering
#' for details on the original implementation
#' 
#' @return A data frame, shaped in a similar way as the originally provided 
#' `res_enrich` object, containing two extra columns: `gs_fuzzycluster`, to specify 
#' the identifier of the fuzzy cluster of genesets, and `gs_cluster_status`, which
#' can specify whether the geneset is the "Representative" for that cluster or 
#' a simple "Member".
#' Notably, the number of rows in the returned object can be higher than the
#' original number of rows in `res_enrich`.
#' 
#' @export
#'
#' @examples
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#' # taking a smaller subset
#' res_enrich_subset <- res_enrich[1:100, ]
#' 
#' fuzzy_subset <- gs_fuzzyclustering(
#'   res_enrich = res_enrich_subset,
#'   n_gs = nrow(res_enrich),
#'   gs_ids = NULL,
#'   similarity_matrix = NULL,
#'   similarity_threshold = 0.35,
#'   fuzzy_seeding_initial_neighbors = 3,
#'   fuzzy_multilinkage_rule = 0.5)
#' 
#' # show all genesets members of the first cluster
#' fuzzy_subset[fuzzy_subset$gs_fuzzycluster == "1", ]
#' 
#' # list only the representative clusters
#' head(fuzzy_subset[fuzzy_subset$gs_cluster_status == "Representative", ], 10)
gs_fuzzyclustering <- function(res_enrich,
                               n_gs = nrow(res_enrich),
                               gs_ids = NULL,
                               similarity_matrix = NULL,
                               similarity_threshold = 0.35,
                               fuzzy_seeding_initial_neighbors = 3,
                               fuzzy_multilinkage_rule = 0.5) {
  
  stopifnot(is.numeric(similarity_threshold))
  stopifnot(is.numeric(fuzzy_seeding_initial_neighbors))
  stopifnot(is.numeric(fuzzy_multilinkage_rule))
  stopifnot(similarity_threshold >= 0 & similarity_threshold <= 1)
  stopifnot(fuzzy_seeding_initial_neighbors >= 2)
  stopifnot(fuzzy_multilinkage_rule > 0 & fuzzy_multilinkage_rule <= 1)
  
  if (similarity_threshold < 0.3)
    warning(
      "You selected a small value for the similarity threshold. ",
      "The resulting clusters might be driven by noisy similarity values."
    )
  
  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )
  # TODO: subset here res_enrich AND matrix if provided?
  
  # if simmat not provided, compute it on the fly
  if (is.null(similarity_matrix)) {
    similarity_matrix <- create_kappa_matrix(res_enrich,
                                             n_gs = n_gs,
                                             gs_ids = gs_ids)
    # rownames(similarity_matrix) <- colnames(similarity_matrix) <- res_enrich$gs_id
  }
  
  stopifnot(isSymmetric.matrix(similarity_matrix))
  # stopifnot(all(colnames(similarity_matrix) %in% res_enrich$gs_id))
  stopifnot(is.numeric(similarity_matrix))
  stopifnot(any(similarity_matrix <= 1))
  
  # reworking on the example depicted in https://david.ncifcrf.gov/helps/functional_classification.html#clustering
  # but using gene sets similarity instead of pure gene similarity
  
  # Fuzzy seeding
  fuzzy_seeds <- list()
  seed_id <- 1
  
  for (i_gs in seq_len(nrow(similarity_matrix))) {
    cur_gs <- rownames(similarity_matrix)[i_gs]
    cur_gs_simmat <- similarity_matrix[i_gs, ]
    gs_over_threshold <- cur_gs_simmat >= similarity_threshold
    if (sum(gs_over_threshold) > fuzzy_seeding_initial_neighbors) {
      similar_genesets <- names(cur_gs_simmat)[gs_over_threshold]
      all_seed_gs <- unique(c(cur_gs, similar_genesets)) # if similarity with itself is not reported
      # restricting on the initial seed members
      seed_kappa <- similarity_matrix[rownames(similarity_matrix) %in%
                                        all_seed_gs, colnames(similarity_matrix) %in% all_seed_gs]
      diag(seed_kappa) <- 0
      if (mean(seed_kappa >= similarity_threshold) >= fuzzy_multilinkage_rule) {
        # assign members and name the seed
        fuzzy_seeds[[seed_id]] <- similar_genesets
        names(fuzzy_seeds)[seed_id] <- cur_gs
        seed_id <- seed_id + 1
      }
    }
  }
  
  # Iteratively merging the above qualified fuzzy seeds
  fuzzy_seeds_unique <- unique(fuzzy_seeds)
  i_seed <- 1
  j_seed <- i_seed + 1
  while (i_seed < length(fuzzy_seeds_unique)) {
    shared_genesets <- intersect(fuzzy_seeds_unique[[i_seed]], fuzzy_seeds_unique[[j_seed]])
    all_genesets <- union(fuzzy_seeds_unique[[i_seed]], fuzzy_seeds_unique[[j_seed]])
    if (length(shared_genesets)/length(all_genesets) > fuzzy_multilinkage_rule & i_seed != j_seed) {
      # message("merging...")
      # is sharing majority and not the same seed, merge and drop the seed
      fuzzy_seeds_unique[[i_seed]] <- all_genesets
      fuzzy_seeds_unique[[j_seed]] <- NULL
      i_seed <- 1
      j_seed <- i_seed + 1
    } else if (j_seed < length(fuzzy_seeds_unique)) {
      # message("updating j")
      j_seed <- j_seed + 1
    } else {
      # message("updating i")
      i_seed <- i_seed + 1
      j_seed <- 1
    }
  }
  
  # Rescuing the singletons
  gs_singletons <-
    res_enrich[ !(res_enrich$gs_id %in% unlist(fuzzy_seeds_unique)), "gs_id"] # TODO: gs_id
  for (gs_singleton in gs_singletons) {
    fuzzy_seeds_unique[[gs_singleton]] <- gs_singleton
  }
  names(fuzzy_seeds_unique) <- seq_len(length(fuzzy_seeds_unique))
  # return(fuzzy_seeds_unique)
  
  # Making the list into matrix
  fuzzy_clusters_mat <- matrix(FALSE,
                               nrow = nrow(res_enrich),
                               ncol = length(fuzzy_seeds_unique),
                               dimnames = list(res_enrich[,"gs_id"],
                                               names(fuzzy_seeds_unique)))
  for (i_clu in names(fuzzy_seeds_unique)) {
    # message(i_clu)
    cluster_gs <- fuzzy_seeds_unique[[i_clu]]
    fuzzy_clusters_mat[cluster_gs, i_clu] <- TRUE
  }
  
  # return(fuzzy_clusters_mat)
  
  gs_list <- list()
  for (i_gs in rownames(fuzzy_clusters_mat)) {
    gs_list[[i_gs]] <- which(fuzzy_clusters_mat[i_gs, ])
  }
  
  res_enrich_norownames <- res_enrich
  rownames(res_enrich_norownames) <- NULL
  buildup_res_enrich <- c()
  for (i_gs in seq_len(nrow(res_enrich_norownames))) {
    cur_gs <- res_enrich_norownames[i_gs, ]
    cur_clusters <- gs_list[[cur_gs[, "gs_id"]]]
    for (i_fc in cur_clusters) {
      buildup_res_enrich <- rbind(
        buildup_res_enrich,
        data.frame(cur_gs,
                   gs_fuzzycluster = i_fc)
      )
    }
  }
  
  res_enrich <- buildup_res_enrich
  
  gs_mostsig <- res_enrich %>%
    group_by(.data$gs_fuzzycluster) %>%
    arrange(.data$gs_pvalue) %>%
    slice(1) %>%
    select(.data$gs_id, .data$gs_fuzzycluster)
  # handling this with tuple to account where the gene set is most representative
  gs_mostsig_tuple <- paste0(gs_mostsig$gs_id, "|", gs_mostsig$gs_fuzzycluster)
  
  res_enrich <- res_enrich[order(res_enrich$gs_fuzzycluster,
                                 res_enrich$gs_pvalue, decreasing = FALSE), ]
  
  re_tuple <- paste0(res_enrich$gs_id, "|", res_enrich$gs_fuzzycluster)
  gs_cluster_status <- re_tuple %in% gs_mostsig_tuple
  res_enrich$gs_cluster_status <- ifelse(gs_cluster_status, "Representative", "Member")
  
  # restoring row names making them unique
  rownames(res_enrich) <- make.unique(res_enrich$gs_id, sep = "_")
  
  return(res_enrich)
}



