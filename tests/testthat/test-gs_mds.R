context("MDS plot of genesets works")

test_that("mds plot works", {
  p <- gs_mds(res_enrich = res_enrich_IFNg_vs_naive,
              res_de = res_macrophage_IFNg_vs_naive,
              annotation_obj = anno_df,
              similarity_measure = "kappa_matrix",
              mds_k = 2,
              mds_labels = 10,
              mds_colorby = "z_score")
  expect_is(p, "gg")
  myset <- res_enrich_IFNg_vs_naive$gs_id[c(1,10,20)]
})

test_that("mds plot with specified sets works", {
  mygenesets <- res_enrich_IFNg_vs_naive$gs_id[c(1,10,20)]
  p <- gs_mds(res_enrich = res_enrich_IFNg_vs_naive,
              res_de = res_macrophage_IFNg_vs_naive,
              annotation_obj = anno_df,
              gs_labels = mygenesets,
              mds_colorby = "z_score")
  expect_is(p, "gg")
  expect_warning(gs_mds(res_enrich = res_enrich_IFNg_vs_naive,
                        res_de = res_macrophage_IFNg_vs_naive,
                        annotation_obj = anno_df,
                        gs_labels = "a_random_name",
                        mds_colorby = "z_score"))
})


