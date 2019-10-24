context("Testing the GO volcano plot and related functionality")

test_that("Plot is generated", {
  expect_error(
    go_volcano(topgoDE_macrophage_IFNg_vs_naive,
               labels_to_use = "Term",
               pvals_to_use = "p.value_elim"))

  res_enrich_withscores <- get_aggrscores(topgoDE_macrophage_IFNg_vs_naive,
                                          res_macrophage_IFNg_vs_naive,
                                          annotation_obj = anno_df,
                                          aggrfun = mean)
  expect_is(
    gs_volcano(res_enrich_withscores,
               labels_to_use = "Term",
               pvals_to_use = "p.value_elim"),
    "gg")
})
