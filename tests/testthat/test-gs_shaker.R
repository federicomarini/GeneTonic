context("Pre-shaking the enrichment results")

test_that("Converting from topGOtable results", {
  res_enrich_IFNg_vs_naive <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
  required_colnames <- c("gs_id", "gs_description", "gs_pvalue", "gs_genes", "gs_de_count", "gs_bg_count")
  expect_true(all(required_colnames %in% colnames(res_enrich_IFNg_vs_naive)))
})

test_that("Converting from clusterProfiler", {

  res_enrich_IFNg_vs_naive_cp <- shake_enrichResult(ego_IFNg_vs_naive)
  required_colnames <- c("gs_id", "gs_description", "gs_pvalue", "gs_genes", "gs_de_count", "gs_bg_count")
  expect_true(all(required_colnames %in% colnames(res_enrich_IFNg_vs_naive_cp)))

  expect_error(shake_enrichResult(topgoDE_macrophage_IFNg_vs_naive))
  expect_error(shake_enrichResult(as.data.frame(ego_IFNg_vs_naive)))
})
