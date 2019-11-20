context("Testing create_kappa_matrix and related functionality")

test_that("Kappa matrix is created", {
  kmat <- create_kappa_matrix(res_enrich_IFNg_vs_naive, n_gs = 30)
  expect_true(all(diag(kmat) == 1))

  kmat2 <- create_kappa_matrix(res_enrich_IFNg_vs_naive,
                               n_gs = 20,
                               gs_ids = res_enrich_IFNg_vs_naive$gs_id[21:30])
  expect_true(identical(kmat, kmat2))
})
