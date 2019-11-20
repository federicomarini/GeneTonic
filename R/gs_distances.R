#' Compute the kappa matrix for enrichment results
#'
#' Compute the kappa matrix for enrichment results, as a measure of overlap
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to see the
#' formatting requirements.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included (from the top ranked ones). Defaults to the number of rows of
#' `res_enrich`
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be included, additionally to
#' the ones specified via `n_gs`. Defaults to NULL.
#'
#' @return A matrix with the kappa scores between gene sets
#'
#' @seealso [gs_mds()]
#'
#' @export
#'
#' @examples
#' # TODO
create_kappa_matrix <- function(res_enrich,
                                n_gs = nrow(res_enrich),
                                gs_ids = NULL

                                ) {

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )

  enrich2list <- lapply(gs_to_use, function(gs) {
    go_genes <- res_enrich[gs, "gs_genes"]
    go_genes <- unlist(strsplit(go_genes, ","))
    return(go_genes)
  })
  names(enrich2list) <- res_enrich[gs_to_use, "gs_id"]

  # creating binary matrix with gs-gene occurrency
  all_genes <- unique(unlist(enrich2list))
  binary_mat <- matrix(0,
                       nrow = length(gs_to_use),
                       ncol = length(all_genes),
                       dimnames = list(names(enrich2list), all_genes))
  # fill in the occurrency per gene set
  for (i in seq_len(length(gs_to_use))) {
    current_genes <- enrich2list[[i]]
    binary_mat[i, colnames(binary_mat) %in% current_genes] <- 1
  }

  # creating Kappa Matrix
  kappa_matrix <- matrix(0, nrow = nrow(binary_mat), ncol = nrow(binary_mat),
                         dimnames = list(rownames(binary_mat), rownames(binary_mat)))
  diag(kappa_matrix) <- 1
  N <- length(gs_to_use)

  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      genes_i <- enrich2list[[i]]
      genes_j <- enrich2list[[j]]

      set_both <- length(intersect(genes_i, genes_j))
      set_none <- sum(!all_genes %in% genes_i & !all_genes %in% genes_j)
      set_only_i <- sum(all_genes %in% genes_i  & !all_genes %in% genes_j)
      set_only_j <- sum(!all_genes %in% genes_i & all_genes %in% genes_j)

      tot <- sum(set_none, set_only_i, set_only_j, set_both)

      observed <- (set_both + set_none) / tot
      chance <- (set_both + set_only_i) * (set_both + set_only_j) + (set_only_j + set_none) * (set_only_i + set_none)
      chance <- chance / tot^2
      kappa_matrix[j,i] <- kappa_matrix[i,j] <- (observed - chance) / (1 - chance)
    }
  }

  return(kappa_matrix)
}





#' Compute the overlap matrix for enrichment results
#'
#' Compute the overlap matrix for enrichment results, based on the Jaccard Index
#' between each pair of sets
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to see the
#' formatting requirements.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included (from the top ranked ones). Defaults to the number of rows of
#' `res_enrich`
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be included, additionally to
#' the ones specified via `n_gs`. Defaults to NULL.
#' @param return_sym Logical, whether to return the symmetrical matrix or just the
#' upper triangular - as needed by [enrichment_map()], for example.
#'
#' @return A matrix with the kappa scores between gene sets
#'
#' @seealso [gs_mds()], [enrichment_map()]
#'
#' @export
#'
#' @examples
#' # TODO
create_jaccard_matrix <- function(res_enrich,
                                  n_gs = nrow(res_enrich),
                                  gs_ids = NULL,
                                  return_sym = FALSE) {

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )

  enrich2list <- lapply(gs_to_use, function(gs) {
    go_genes <- res_enrich[gs, "gs_genes"]
    go_genes <- unlist(strsplit(go_genes, ","))
    return(go_genes)
  })
  names(enrich2list) <- res_enrich[gs_to_use, "gs_id"]

  n <- length(enrich2list)
  overlap_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- names(enrich2list)

  for (i in 1:n) {
    # no need to work on full mat, it is simmetric
    for (j in i:n) {
      overlap_matrix[i, j] <-
        overlap_jaccard_index(unlist(enrich2list[gs_to_use[i]]),
                              unlist(enrich2list[gs_to_use[j]]))
    }
  }

  if(return_sym) {
    # symmetryze :)
    overlap_matrix[lower.tri(overlap_matrix)] = t(overlap_matrix)[lower.tri(overlap_matrix)]
  }
  return(overlap_matrix)
}
