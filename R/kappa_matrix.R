


create_kappa_matrix <- function(res_enrich,
                                genes_colname = "genes",
                                genesetname_colname = "Term",
                                genesetid_colname = "GO.ID",
                                use_names = FALSE,
                                use_active_snw_genes = FALSE) {

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
  kappa_mat <- matrix(0, nrow = nrow(binary_mat), ncol = nrow(binary_mat),
                      dimnames = list(rownames(binary_mat), rownames(binary_mat)))
  diag(kappa_mat) <- 1

  for (i in 1:(nrow(binary_mat) - 1)) {
    for (j in (i+1):nrow(binary_mat)) {
      gene_vec_i <- binary_mat[i, ]
      gene_vec_j <- binary_mat[j, ]
      cross_tbl <- table(gene_vec_i, gene_vec_j)

      observed<- (cross_tbl[1, 1] + cross_tbl[2, 2]) / sum(cross_tbl)
      chance <- (sum(cross_tbl[1, ]) * sum(cross_tbl[, 1]) + sum(cross_tbl[2, ]) * sum(cross_tbl[, 2])) / (sum(cross_tbl)^2)
      kappa_mat[j, i] <- kappa_mat[i, j] <- (observed - chance) / (1 - chance)
    }
  }

  return(kappa_mat)
}




