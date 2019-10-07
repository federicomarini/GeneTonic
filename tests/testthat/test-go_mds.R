context("MDS plot of genesets works")

test_that("mds plot works", {
  p <- go_mds(res_enrich = topgoDE_macrophage_IFNg_vs_naive,
              res_de = res_macrophage_IFNg_vs_naive,
              annotation_obj = anno_df,
              genes_colname = "genes",
              genesetname_colname = "Term",
              genesetid_colname = "GO.ID",
              genes_separator = ",",
              similarity_measure = "kappa_matrix",
              mds_k = 2,
              mds_labels = 10,
              mds_colorby = "z_score")
  expect_is(p, "gg")
})
