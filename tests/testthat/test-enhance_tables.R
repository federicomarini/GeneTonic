context("Testing enhance_tables and related functionality")

test_that("Enhanced table is created", {
  p <- enhance_table(res_enrich_IFNg_vs_naive,
                     res_macrophage_IFNg_vs_naive,
                     annotation_obj = anno_df,
                     n_gs = 50,
                     chars_limit = 60)
  expect_is(p, "gg")

  pl <- ggplotly(p)
  expect_is(pl, "plotly")
  expect_is(pl, "htmlwidget")

  p2 <- enhance_table(res_enrich_IFNg_vs_naive,
                     res_macrophage_IFNg_vs_naive,
                     annotation_obj = anno_df,
                     n_gs = 5,
                     chars_limit = 60,
                     plot_title = "My custom title")
  expect_is(p2, "gg")


})
