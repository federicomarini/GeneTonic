context("Testing enrichmentmap and related functionality")

test_that("Graph is generated", {
  g <- enrichment_map(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                      res_de = res_macrophage_IFNg_vs_naive,
                      annotation_obj = anno_df,
                      n_gs = 50)
  expect_is(g, "igraph")
  # TODOTODO, expect I get htmlwidgets...
})
