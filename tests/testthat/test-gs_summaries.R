context("Testing the gene set summaries and related functionality")

test_that("summary_heat plot is generated", {
  p <- gs_summary_heat(res_enrich = res_enrich_IFNg_vs_naive,
                       res_de = res_macrophage_IFNg_vs_naive,
                       annotation_obj = anno_df,
                       n_gs = 30)
  expect_is(p, "gg")
})

test_that("summary plots are generated", {
  expect_error(gs_summary_overview(res_enrich_IFNg_vs_naive))
  expect_error(gs_summary_overview_pair(res_enrich_IFNg_vs_naive))
  res_enrich_withscores <- get_aggrscores(res_enrich_IFNg_vs_naive,
                                          res_macrophage_IFNg_vs_naive,
                                          annotation_obj = anno_df,
                                          aggrfun = mean)

  # generating a shuffled dataset
  res_enrich2 <- res_enrich_withscores[1:20, ]
  set.seed(42)
  shuffled_ones <- sample(seq_len(20)) # to generate permuted p-values
  res_enrich2$gs_pvalue <- res_enrich2$gs_pvalue[shuffled_ones]
  res_enrich2$z_score <- res_enrich2$z_score[shuffled_ones]
  res_enrich2$aggr_score <- res_enrich2$aggr_score[shuffled_ones]

  p1 <- gs_summary_overview(res_enrich_withscores)
  expect_is(p1, "gg")

  p2 <- gs_summary_overview_pair(res_enrich_withscores, res_enrich2)
  expect_is(p2, "gg")

  res_enrich2 <- res_enrich_withscores[1:42, ]
  res_enrich3 <- res_enrich_withscores[1:42, ]
  res_enrich4 <- res_enrich_withscores[1:42, ]

  set.seed(2*42)
  shuffled_ones_2 <- sample(seq_len(42)) # to generate permuted p-values
  res_enrich2$gs_pvalue <- res_enrich2$gs_pvalue[shuffled_ones_2]
  res_enrich2$z_score <- res_enrich2$z_score[shuffled_ones_2]
  res_enrich2$aggr_score <- res_enrich2$aggr_score[shuffled_ones_2]

  set.seed(3*42)
  shuffled_ones_3 <- sample(seq_len(42)) # to generate permuted p-values
  res_enrich3$gs_pvalue <- res_enrich3$gs_pvalue[shuffled_ones_3]
  res_enrich3$z_score <- res_enrich3$z_score[shuffled_ones_3]
  res_enrich3$aggr_score <- res_enrich3$aggr_score[shuffled_ones_3]

  set.seed(4*42)
  shuffled_ones_4 <- sample(seq_len(42)) # to generate permuted p-values
  res_enrich4$gs_pvalue <- res_enrich4$gs_pvalue[shuffled_ones_4]
  res_enrich4$z_score <- res_enrich4$z_score[shuffled_ones_4]
  res_enrich4$aggr_score <- res_enrich4$aggr_score[shuffled_ones_4]

  compa_list <- list(
    scenario2 = res_enrich2,
    scenario3 = res_enrich3,
    scenario4 = res_enrich4
  )

  p3a <- gs_horizon(res_enrich_withscores,
             compared_res_enrich_list = compa_list,
             n_gs = 50,
             sort_by = "clustered")
  p3b <- gs_horizon(res_enrich_withscores,
             compared_res_enrich_list = compa_list,
             n_gs = 20,
             sort_by = "first_set")
  expect_is(p3a, "gg")
  expect_is(p3b, "gg")
})
