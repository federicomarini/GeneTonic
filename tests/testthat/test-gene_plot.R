context("Testing function for plotting the gene expression levels")

test_that("Basic gene plot is generated", {
  p <- gene_plot(dds = dds_macrophage,
                 gene = "ENSG00000285982",
                 assay = "counts",
                 intgroup = "condition",
                 annotation_obj = anno_df,
                 transform = TRUE,
                 labels_repel = TRUE)
  expect_is(p, "gg")

  p2_noanno_normallabels_untransformed <- gene_plot(dds = dds_macrophage,
                 gene = "ENSG00000285982",
                 assay = "counts",
                 intgroup = "condition",
                 transform = FALSE,
                 labels_repel = FALSE)
  expect_is(p2_noanno_normallabels_untransformed, "gg")
})

test_that("Enforcing plot types", {
  p_jitter <- gene_plot(dds = dds_macrophage,
                 gene = "ENSG00000285982",
                 assay = "counts",
                 intgroup = "condition",
                 plot_type = "jitteronly")
  p_boxplot <- gene_plot(dds = dds_macrophage,
                 gene = "ENSG00000285982",
                 assay = "counts",
                 intgroup = "condition",
                 plot_type = "boxplot")
  p_violin <- gene_plot(dds = dds_macrophage,
                 gene = "ENSG00000285982",
                 assay = "counts",
                 intgroup = "condition",
                 plot_type = "violin")
  p_sina <- gene_plot(dds = dds_macrophage,
                 gene = "ENSG00000285982",
                 assay = "counts",
                 intgroup = "condition",
                 plot_type = "sina")
  expect_is(p_jitter, "gg")
  expect_is(p_boxplot, "gg")
  expect_is(p_violin, "gg")
  expect_is(p_sina, "gg")
})
