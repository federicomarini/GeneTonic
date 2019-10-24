context("Testing enhance_tables and related functionality")

test_that("Enhanced table is created", {
  p <- enhance_table(topgoDE_macrophage_IFNg_vs_naive,
                     res_macrophage_IFNg_vs_naive,
                     annotation_obj = anno_df,
                     n_gs = 50,
                     genes_colname = "genes",
                     genesetname_colname = "Term",
                     genesetid_colname = "GO.ID",
                     chars_limit = 60)
  expect_is(p, "gg")
  # TODOTODO, expect I get htmlwidgets...
})
