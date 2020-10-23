# similarity_matrix <- km_macro

#' Title
#' 
#' TODO
#'
#' @param res_enrich TODO
#' @param n_gs TODO
#' @param gs_ids TODO
#' @param similarity_matrix TODO
#' @param similarity_threshold TODO
#' @param fuzzy_seeding_initial_neighbors TODO
#' @param fuzzy_multilinkage_rule TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_fuzzyclustering <- function(res_enrich,
                               n_gs = nrow(res_enrich),
                               gs_ids = NULL,
                               similarity_matrix = NULL,
                               similarity_threshold = 0.35,
                               fuzzy_seeding_initial_neighbors = 3,
                               fuzzy_multilinkage_rule = 0.5) {
  
  # checks: simmat is matrix
  # checks: simmat is symmetric
  # checks: simmat is containing the same set of ids
  # checks: simmat is numeric, and all are less or equal 1
  
  stopifnot(is.numeric(similarity_threshold))
  stopifnot(is.numeric(fuzzy_seeding_initial_neighbors))
  stopifnot(is.numeric(fuzzy_multilinkage_rule))
  stopifnot(similarity_threshold >= 0 & similarity_threshold <= 1)
  stopifnot(fuzzy_multilinkage_rule > 0 & fuzzy_multilinkage_rule <= 1)
  
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



