context("Testing ggs_graph and related functionality")

test_that("Graph is generated", {
  g <- ggs_graph(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                 res_de = res_macrophage_IFNg_vs_naive,
                 annotation_obj = anno_df,
                 n_gs = 20)
  expect_is(g, "igraph")
  g2 <- ggs_graph(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                  res_de = res_macrophage_IFNg_vs_naive,
                  annotation_obj = anno_df,
                  n_gs = 30,
                  prettify = FALSE)
  expect_is(g2, "igraph")
})
