context("Testing gene set scoring")

test_that("Scores are calculated and plotted", {
  gss_mat <- gs_scores(se = vst_macrophage,
                       res_de = res_macrophage_IFNg_vs_naive,
                       res_enrich = topgoDE_macrophage_IFNg_vs_naive,
                       annotation_obj = anno_df,
                       genes_colname = "genes",
                       genesetname_colname = "Term",
                       genesetid_colname = "GO.ID")
  expect_is(gss_mat, "matrix")
  p <- gs_ggheatmap(gss_mat)
  expect_is(p, "gg")
})
