context("Testing the gene set radar plot and co.")

test_that("radar plot is generated", {
  expect_warning(gs_radar(res_enrich_IFNg_vs_naive))
  res_enrich_withscores <- get_aggrscores(res_enrich_IFNg_vs_naive,
                                          res_macrophage_IFNg_vs_naive,
                                          annotation_obj = anno_df,
                                          aggrfun = mean)
  p <- gs_radar(res_enrich_withscores)
  expect_is(p, "plotly")
})
