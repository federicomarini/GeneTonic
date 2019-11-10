context("Testing create_kappa_matrix and related functionality")

test_that("Kappa matrix is created", {
  kmat <- create_kappa_matrix(res_enrich_IFNg_vs_naive)
  expect_true(all(diag(kmat) == 1))
})
