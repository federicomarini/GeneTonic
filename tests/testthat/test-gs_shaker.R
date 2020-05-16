context("Pre-shaking the enrichment results")

test_that("Converting from topGOtable results", {
  res_enrich_IFNg_vs_naive <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
  required_colnames <- c("gs_id", "gs_description", "gs_pvalue", "gs_genes", "gs_de_count", "gs_bg_count")
  expect_true(all(required_colnames %in% colnames(res_enrich_IFNg_vs_naive)))

  topgo_not_all_columns <- topgoDE_macrophage_IFNg_vs_naive[,-1]
  expect_error(
    shake_topGOtableResult(topgo_not_all_columns)
  )
  expect_error(
    shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive, p_value_column = "p.value_weight")
  )
  topgo_nogenes <- topgoDE_macrophage_IFNg_vs_naive[,-9]
  expect_error(
    shake_topGOtableResult(topgo_nogenes)
  )
})

test_that("Converting from clusterProfiler", {

  res_enrich_IFNg_vs_naive_cp <- shake_enrichResult(ego_IFNg_vs_naive)
  required_colnames <- c("gs_id", "gs_description", "gs_pvalue", "gs_genes", "gs_de_count", "gs_bg_count")
  expect_true(all(required_colnames %in% colnames(res_enrich_IFNg_vs_naive_cp)))

  expect_error(shake_enrichResult(topgoDE_macrophage_IFNg_vs_naive))
  expect_error(shake_enrichResult(as.data.frame(ego_IFNg_vs_naive)))

  ego_mod <- ego_IFNg_vs_naive
  ego_mod@result$geneID <- NULL
  expect_error(shake_enrichResult(ego_mod))
})

test_that("Converting from the output of DAVID", {
  david_output <- system.file("extdata", "david_output_chart_BPonly_ifng_vs_naive.txt", package = "GeneTonic")
  res_enrich_IFNg_vs_naive_david <- shake_davidResult(david_output)
  required_colnames <- c("gs_id", "gs_description", "gs_pvalue", "gs_genes", "gs_de_count", "gs_bg_count")
  expect_true(all(required_colnames %in% colnames(res_enrich_IFNg_vs_naive_david)))
  
  expect_error(shake_davidResult("non_existing_file.txt"))
  expect_error(shake_davidResult(topgoDE_macrophage_IFNg_vs_naive))
  expect_error(shake_davidResult(as.data.frame(ego_IFNg_vs_naive)))
})
