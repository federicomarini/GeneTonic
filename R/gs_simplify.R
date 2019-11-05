#' Simplify results from functional enrichment analysis
#'
#' Simplify results from functional enrichment analysis, removing genesets that
#' are redundant to enhance interpretation of the results
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to see the
#' formatting requirements.
#' @param gs_overlap Numeric value, which defines the threshold for removing
#' terms that present an overlap greater than the specified value. Changing its
#' value can control the granularity of how redundant terms are removed from the
#' original `res_enrich` for the next steps, e.g. plotting this via [gs_volcano()]
#' @param genes_colname Character, specifying which column of the `res_enrich`
#' object contains the genes assigned to each gene set, detected as differentially
#' expressed. Defaults to `genes`.
#' @param genesetname_colname Character, specifies which column of the `res_enrich`
#' object contains a description of the gene set. Defaults to `Term`.
#' @param genesetid_colname Character, specifies which column of the `res_enrich`
#' object contains a unique identifier of the gene set. Defaults to `GO.ID`.
#' @param genes_separator Character, specifying which separator is used in the
#' column defined by `genes_colname` to split the character of features.
#'
#' @return A `data.frame` with a subset of the original gene sets
#'
#' @seealso [gs_volcano()] and [ggs_graph()] can e.g. show an overview on the
#' simplified table of gene sets
#'
#' @export
#'
#' @examples
#' # TODO
gs_simplify <- function(res_enrich,
                        gs_overlap = 0.75,
                        genes_colname = "genes",
                        genesetname_colname = "Term",
                        genesetid_colname = "GO.ID",
                        genes_separator = ","
                        ) {

  genelists <- lapply(seq_len(nrow(res_enrich)), function(gs) {
    cur_set <- res_enrich[[genes_colname]][gs]
    gene_vec <- unlist(strsplit(cur_set, genes_separator))
  })
  names(genelists) <- res_enrich[[genesetid_colname]]

  ol_mat <- matrix(0, nrow(res_enrich), nrow(res_enrich),
                   dimnames = list(res_enrich[[genesetname_colname]],
                                   res_enrich[[genesetname_colname]]))
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


