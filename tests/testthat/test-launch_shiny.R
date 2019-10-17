test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

library("GeneTonic")

context("Check that the Shiny app is generated")

dds <- DESeq2::makeExampleDESeqDataSet(n = 100, m = 8)

test_that("Shiny app is generated", {
  # expect_is(GeneTonic(), "shiny.appobj")
  # expect_is(GeneTonic(dds), "shiny.appobj")
  # expect_is(GeneTonic(dds_macrophage, res_macrophage_IFNg_vs_naive), "shiny.appobj")
  expect_is(GeneTonic(dds_macrophage, res_macrophage_IFNg_vs_naive, topgoDE_macrophage_IFNg_vs_naive, annotation_obj = anno_df), "shiny.appobj")
})
