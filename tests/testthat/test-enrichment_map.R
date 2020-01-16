context("Testing enrichmentmap and related functionality")

test_that("Graph is generated", {
  g <- enrichment_map(res_enrich = res_enrich_IFNg_vs_naive,
                      res_de = res_macrophage_IFNg_vs_naive,
                      annotation_obj = anno_df,
                      n_gs = 50)
  expect_is(g, "igraph")

  pl <- ggplotly(g)
  expect_is(pl, "plotly")
  expect_is(pl, "htmlwidgets")
})
