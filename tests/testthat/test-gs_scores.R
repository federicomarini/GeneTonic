context("Testing gene set scoring")

test_that("Scores are calculated and plotted", {
  gss_mat <- gs_scores(vst_macrophage,
                       res_macrophage_IFNg_vs_naive,
                       topgoDE_macrophage_IFNg_vs_naive,
                       genes_colname = "genes",
                       genesetname_colname = "Term",
                       genesetid_colname = "GO.ID",
                       annotation_obj = anno_df)
  expect_is(gss_mat, "matrix")
  p <- gs_ggheatmap(gss_mat)
  expect_is(p, "gg")
})
