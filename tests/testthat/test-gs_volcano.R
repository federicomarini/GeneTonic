context("Testing the GO volcano plot and related functionality")

test_that("Plot is generated", {
  expect_error(
    gs_volcano(res_enrich_IFNg_vs_naive))

  res_enrich_withscores <- get_aggrscores(res_enrich_IFNg_vs_naive,
                                          res_macrophage_IFNg_vs_naive,
                                          annotation_obj = anno_df,
                                          aggrfun = mean)
  expect_is(gs_volcano(res_enrich_withscores), "gg")
})

test_that("mds plot with custom genesets", {
  mygenesets <- res_enrich_IFNg_vs_naive$gs_id[c(1, 10, 20)]
  res_enrich_withscores <- get_aggrscores(res_enrich_IFNg_vs_naive,
                                          res_macrophage_IFNg_vs_naive,
                                          annotation_obj = anno_df,
                                          aggrfun = mean)
  expect_is(
    gs_volcano(res_enrich_withscores,
               gs_ids = mygenesets),
    "gg")
})
