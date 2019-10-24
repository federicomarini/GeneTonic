context("Testing create_kappa_matrix and related functionality")

test_that("Kappa matrix is created", {
  kmat <- create_kappa_matrix(topgoDE_macrophage_IFNg_vs_naive,
                              genes_colname = "genes",
                              genesetname_colname = "Term",
                              genesetid_colname = "GO.ID"
  )
  expect_true(all(diag(kmat) == 1))
})
