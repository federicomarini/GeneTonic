


#' Title TODO
#'
#' @param res_enrich TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
create_kappa_matrix <- function(res_enrich,
                                genes_colname = "genes",
                                genesetname_colname = "Term",
                                genesetid_colname = "GO.ID"
                                ) {

  # initial checks
  ## more than one row
  ## genes column is present

  genelists <- lapply(seq_len(nrow(res_enrich)), function(gs) {
    cur_set <- res_enrich[[genes_colname]][gs]
    gene_vec <- unlist(strsplit(cur_set,","))
  })

  names(genelists) <- res_enrich[[genesetid_colname]]

  # creating binary matrix with gs-gene occurrency
  all_genes <- unique(unlist(genelists))
  binary_mat <- matrix(0,
                       nrow = nrow(res_enrich),
                       ncol = length(all_genes),
                       dimnames = list(res_enrich[[genesetid_colname]], all_genes))
  # fill in the occurrency per gene set
  for (i in seq_len(nrow(res_enrich))) {
    current_genes <- genelists[[i]]
    binary_mat[i, colnames(binary_mat) %in% current_genes] <- 1
  }

  # creating Kappa Matrix
  kappa_matrix <- matrix(0, nrow = nrow(binary_mat), ncol = nrow(binary_mat),
                      dimnames = list(rownames(binary_mat), rownames(binary_mat)))
  diag(kappa_matrix) <- 1
  N <- nrow(res_enrich)

  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      genes_i <- genelists[[i]]
      genes_j <- genelists[[j]]

      set_both <- length(intersect(genes_i, genes_j))
      set_none <- sum(!all_genes %in% genes_i & !all_genes %in% genes_j)
      set_only_i <- sum(all_genes %in% genes_i  & !all_genes %in% genes_j)
      set_only_j <- sum(!all_genes %in% genes_i & all_genes %in% genes_j)

      tot <- sum(set_none, set_only_i, set_only_j, set_both)

      observed <- (set_both + set_none)/tot
      chance <- (set_both + set_only_i)*(set_both + set_only_j) + (set_only_j + set_none)*(set_only_i + set_none)
      chance <- chance/tot^2
      kappa_matrix[j,i] <- kappa_matrix[i,j] <- (observed - chance)/(1 - chance)
    }
  }

  return(kappa_matrix)
}




