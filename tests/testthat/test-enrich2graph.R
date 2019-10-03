test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

context("Testing enrich2graph and related functionality")

test_that("Graph is generated", {
  g <- enrich2graph(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                    res_de = res_macrophage_IFNg_vs_naive,
                    n_nodes = 50,
                    annotation_obj = anno_df)
  expect_is(g, "igraph")
})
