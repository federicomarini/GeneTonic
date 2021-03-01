context("Testing ggs_graph and related functionality")

test_that("Graph is generated", {
  g <- ggs_graph(res_enrich = res_enrich_IFNg_vs_naive,
                 res_de = res_macrophage_IFNg_vs_naive,
                 annotation_obj = anno_df,
                 n_gs = 20)
  expect_is(g, "igraph")
  g2 <- ggs_graph(res_enrich = res_enrich_IFNg_vs_naive,
                  res_de = res_macrophage_IFNg_vs_naive,
                  annotation_obj = anno_df,
                  n_gs = 30,
                  prettify = FALSE)
  expect_is(g2, "igraph")

  gtl_macrophage <- GeneTonic_list(dds = dds_macrophage,
                                   res_de = res_macrophage_IFNg_vs_naive,
                                   res_enrich = res_enrich_IFNg_vs_naive,
                                   annotation_obj = anno_df)
  g3 <- ggs_graph(gtl = gtl_macrophage,
                  n_gs = 20)
  expect_is(g3, "igraph")

  expect_true(identical_graphs(g, g3))

  alt_pal <- scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.4)
  g3 <- ggs_graph(res_enrich = res_enrich_IFNg_vs_naive,
                  res_de = res_macrophage_IFNg_vs_naive,
                  annotation_obj = anno_df,
                  n_gs = 10,
                  genes_graph_colpal = alt_pal)
  expect_is(g3, "igraph")

  expect_error(
    ggs_graph(res_enrich = res_enrich_IFNg_vs_naive,
              res_de = res_macrophage_IFNg_vs_naive,
              annotation_obj = anno_df,
              genes_graph_colpal = list("blue", "red"),
              n_gs = 20)
  )

  expect_error(
    ggs_graph(res_enrich = res_enrich_IFNg_vs_naive,
              res_de = res_macrophage_IFNg_vs_naive,
              annotation_obj = anno_df,
              genes_graph_colpal = c("blue", "whitesss", "red"),
              n_gs = 20)
  )
})
