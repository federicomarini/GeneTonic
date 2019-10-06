context("Testing enrich2graph and related functionality")

test_that("Graph is generated", {
  g <- enrich2graph(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                    res_de = res_macrophage_IFNg_vs_naive,
                    n_nodes = 20,
                    annotation_obj = anno_df)
  expect_is(g, "igraph")
  g2 <- enrich2graph(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                     res_de = res_macrophage_IFNg_vs_naive,
                     n_nodes = 30,
                     annotation_obj = anno_df,
                     prettify = FALSE)
  expect_is(g2, "igraph")
})
