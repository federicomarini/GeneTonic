context("Testing gene set heatmap and related functionality")

test_that("Geneset heatmap is created", {
  cur_gsid <- topgoDE_macrophage_IFNg_vs_naive$GO.ID[1]
  myvst <- vst(dds_macrophage)
  p <- gs_heatmap(myvst,
                  res_macrophage_IFNg_vs_naive,
                  topgoDE_macrophage_IFNg_vs_naive,
                  geneset_id = cur_gsid,
                  genes_colname = "genes",
                  genesetname_colname = "Term",
                  genesetid_colname = "GO.ID",
                  annotation_obj = anno_df,
                  FDR = 0.05,
                  de_only = FALSE,
                  cluster_rows = TRUE, # TODOTODO: options for the heatmap go on left side, as could be common to more!
                  cluster_cols = TRUE,
                  center_mean = TRUE,
                  scale_row = TRUE
                  # TODOTODO: use ellipsis for passing params to pheatmap?
                  # TODOTODO: option to just return the underlying data?s
                  # TODOTODO: options to subset to specific samples?
  )
  expect_is(p, "pheatmap")
  p2 <- gs_heatmap(myvst,
                  res_macrophage_IFNg_vs_naive,
                  topgoDE_macrophage_IFNg_vs_naive,
                  geneset_id = cur_gsid,
                  genes_colname = "genes",
                  genesetname_colname = "Term",
                  genesetid_colname = "GO.ID",
                  annotation_obj = anno_df,
                  FDR = 0.05,
                  de_only = TRUE,
                  cluster_rows = TRUE, # TODOTODO: options for the heatmap go on left side, as could be common to more!
                  cluster_cols = TRUE,
                  center_mean = TRUE,
                  scale_row = TRUE
                  # TODOTODO: use ellipsis for passing params to pheatmap?
                  # TODOTODO: option to just return the underlying data?s
                  # TODOTODO: options to subset to specific samples?
  )
  expect_is(p2, "pheatmap")
})
