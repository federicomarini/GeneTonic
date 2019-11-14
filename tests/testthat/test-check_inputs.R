context("Input main parameter checking works")

test_that("Early fails are triggered", {

  # providing a count matrix
  expect_error(GeneTonic(counts(dds_macrophage),
                         res_macrophage_IFNg_vs_naive,
                         res_enrich_IFNg_vs_naive,
                         annotation_obj = anno_df))

  # providing a simple data frame
  expect_error(GeneTonic(dds_macrophage,
                         deseqresult2df(res_macrophage_IFNg_vs_naive),
                         res_enrich_IFNg_vs_naive,
                         annotation_obj = anno_df))

  # providing data frame with missing key columns
  expect_error(GeneTonic(dds_macrophage,
                         res_macrophage_IFNg_vs_naive,
                         res_enrich_IFNg_vs_naive[, -1],
                         annotation_obj = anno_df))

  # providing data frame with missing key columns
  expect_error(GeneTonic(dds_macrophage,
                         res_macrophage_IFNg_vs_naive,
                         res_enrich_IFNg_vs_naive,
                         annotation_obj = anno_df[, -1]))
})
