context("Testing enrichmentmap and related functionality")

test_that("Graph is generated", {
  g <- enrichment_map(res_enrich = res_enrich_IFNg_vs_naive,
                      res_de = res_macrophage_IFNg_vs_naive,
                      annotation_obj = anno_df,
                      n_gs = 50)
  expect_is(g, "igraph")

  library(magrittr)
  vi <- visNetwork::visIgraph(g) %>%
    visOptions(highlightNearest = list(enabled = TRUE,
                                       degree = 1,
                                       hover = TRUE),
               nodesIdSelection = TRUE)

  expect_is(vi, "visNetwork")
  expect_is(vi, "htmlwidget")
})
