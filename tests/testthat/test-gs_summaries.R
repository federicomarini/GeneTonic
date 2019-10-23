context("Testing the gene set summaries and related functionality")

test_that("Plot is generated", {
  p <- gs_summary_heat(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                       res_de = res_macrophage_IFNg_vs_naive,
                       annotation_obj = anno_df,
                       n_gs = 30)
  expect_is(p, "gg")
})

