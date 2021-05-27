#' Simplify results from functional enrichment analysis
#'
#' Simplify results from functional enrichment analysis, removing genesets that
#' are redundant to enhance interpretation of the results
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param gs_overlap Numeric value, which defines the threshold for removing
#' terms that present an overlap greater than the specified value. Changing its
#' value can control the granularity of how redundant terms are removed from the
#' original `res_enrich` for the next steps, e.g. plotting this via [gs_volcano()]
#'
#' @return A `data.frame` with a subset of the original gene sets
#'
#' @seealso [gs_volcano()] and [ggs_graph()] can e.g. show an overview on the
#' simplified table of gene sets
#'
#' @export
#'
#' @examples
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#'
#' dim(res_enrich)
#' res_enrich_simplified <- gs_simplify(res_enrich)
#' dim(res_enrich_simplified)
#' # and then use this further for all other functions expecting a res_enrich
gs_simplify <- function(res_enrich,
                        gs_overlap = 0.75) {
  genelists <- lapply(seq_len(nrow(res_enrich)), function(gs) {
    cur_set <- res_enrich[["gs_genes"]][gs]
    gene_vec <- unlist(strsplit(cur_set, ","))
  })
  names(genelists) <- res_enrich[["gs_id"]]

  ol_mat <- matrix(0, nrow(res_enrich), nrow(res_enrich),
    dimnames = list(
      res_enrich[["gs_description"]],
      res_enrich[["gs_description"]]
    )
  )
  N <- nrow(res_enrich)
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      genes_i <- genelists[[i]]
      genes_j <- genelists[[j]]
      ol_mat[i, j] <- round(sum(genes_i %in% genes_j) / length(genes_i), digits = 3)
    }
  }

  ol_mat2 <- ol_mat
  ol_mat2[upper.tri(ol_mat2)] <- 0

  for (j in seq_len(N)) {
    idx <- which(ol_mat2[, j] > gs_overlap)
    sel_col <- setdiff(idx, j)
    ol_mat2[, sel_col] <- 0
  }

  res_enrich_mod <- res_enrich[colSums(ol_mat2) != 0, ]

  return(res_enrich_mod)
}
