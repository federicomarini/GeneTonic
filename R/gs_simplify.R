

#' Title TODO
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param gs_overlap TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_simplify <- function(res_enrich,
                        gs_overlap = 0.75,
                        genes_colname = "genes",
                        genesetname_colname = "Term"
                        # ,
                        # genesetid_colname = "GO.ID",
                        # genes_separator = ","
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
      ol_mat[i, j] <- round(sum(genes_i %in% genes_j)/length(genes_i), digits = 3)
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


