context("Testing the gene sets dendrogram")

test_that("Gene set dendrogram is created", {
  res_enrich_withscores <- get_aggrscores(res_enrich_IFNg_vs_naive,
                                          res_macrophage_IFNg_vs_naive,
                                          annotation_obj = anno_df,
                                          aggrfun = mean)
  my_dend <- gs_dendro(res_enrich_withscores, n_gs = 30)
  expect_is(my_dend, "dendrogram")
  my_dend2 <- gs_dendro(res_enrich_withscores, n_gs = 10,
                        gs_dist_type = "jaccard", color_branches_by = NULL)
  my_dend3 <- gs_dendro(res_enrich_withscores, n_gs = 20,
                        color_leaves_by = NULL, size_leaves_by = NULL)
  my_dend4 <- gs_dendro(res_enrich_withscores, n_gs = 20,
                        color_leaves_by = NULL, size_leaves_by = NULL,
                        create_plot = FALSE)
  expect_is(my_dend2, "dendrogram")
  expect_is(my_dend3, "dendrogram")
  expect_is(my_dend4, "dendrogram")

  expect_error(
    gs_dendro(res_enrich_withscores, n_gs = 20,
              color_leaves_by = "mean_score")
  )
  expect_error(
    gs_dendro(res_enrich_withscores, n_gs = 20,
              size_leaves_by = "pvalue")
  )
})
