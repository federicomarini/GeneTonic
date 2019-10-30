context("Testing the gene set summaries and related functionality")

test_that("summary_heat plot is generated", {
  p <- gs_summary_heat(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                       res_de = res_macrophage_IFNg_vs_naive,
                       annotation_obj = anno_df,
                       n_gs = 30)
  expect_is(p, "gg")
})

test_that("summary plots are generated", {
  res_enrich_withscores <- get_aggrscores(topgoDE_macrophage_IFNg_vs_naive,
                                          res_macrophage_IFNg_vs_naive,
                                          annotation_obj = anno_df,
                                          aggrfun = mean)

  p1 <- gs_summary_overview(res_enrich_withscores)
  p2 <- gs_summary_overview_pair(res_enrich_withscores)
  p3 <- gs_horizon(res_enrich_withscores)
  expect_is(p1, "gg")
  expect_is(p2, "gg")
  expect_is(p3, "gg")
})
