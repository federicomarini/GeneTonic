context("The happy hour is running!")

test_that("The happy hour starts", {
  my_html_report <- happy_hour(
    dds = dds_macrophage,
    res_de = res_macrophage_IFNg_vs_naive,
    res_enrich = res_enrich_IFNg_vs_naive,
    annotation_obj = anno_df,
    project_id = "testsuite",
    mygenesets = res_enrich_IFNg_vs_naive$gs_description[1:4],
    mygenes = deseqresult2df(res_de)$SYMBOL[1:5],
    force_overwrite = TRUE,
    open_after_creating = FALSE
  )
  expect_is(my_html_report, "character")
})

test_that("Errors and warnings are thrown", {
  expect_error(
    happy_hour(
      dds = dds_macrophage,
      res_de = res_macrophage_IFNg_vs_naive,
      res_enrich = res_enrich_IFNg_vs_naive,
      annotation_obj = anno_df,
      project_id = "testsuite",
      mygenesets = res_enrich_IFNg_vs_naive$gs_description[1:4],
      mygenes = deseqresult2df(res_de)$SYMBOL[1:5],
      output_file = "doc1.html",
      output_format = "pdf_document",
      force_overwrite = TRUE,
      open_after_creating = FALSE
    )
  )
  expect_error(
    happy_hour(
      dds = dds_macrophage,
      res_de = res_macrophage_IFNg_vs_naive,
      res_enrich = res_enrich_IFNg_vs_naive,
      annotation_obj = anno_df,
      project_id = "testsuite",
      mygenesets = res_enrich_IFNg_vs_naive$gs_description[1:4],
      mygenes = deseqresult2df(res_de)$SYMBOL[1:5],
      output_file = "doc1.ppt",
      force_overwrite = TRUE,
      open_after_creating = FALSE
    )
  )

  # file.create("decoy.html")
  # expect_warning(
  #   happy_hour(
  #     dds = dds_macrophage,
  #     res_de = res_macrophage_IFNg_vs_naive,
  #     res_enrich = res_enrich_IFNg_vs_naive,
  #     annotation_obj = anno_df,
  #     project_id = "testsuite",
  #     mygenesets = res_enrich_IFNg_vs_naive$gs_description[1:4],
  #     mygenes = deseqresult2df(res_de)$SYMBOL[1:5],
  #     output_file = "decoy.html",
  #     force_overwrite = TRUE,
  #     open_after_creating = FALSE
  #   )
  # )
  # expect_error(
  #   happy_hour(
  #     dds = dds_macrophage,
  #     res_de = res_macrophage_IFNg_vs_naive,
  #     res_enrich = res_enrich_IFNg_vs_naive,
  #     annotation_obj = anno_df,
  #     project_id = "testsuite",
  #     mygenesets = res_enrich_IFNg_vs_naive$gs_description[1:4],
  #     mygenes = deseqresult2df(res_de)$SYMBOL[1:5],
  #     output_file = "decoy.html",
  #     force_overwrite = FALSE,
  #     open_after_creating = FALSE
  #   )
  # )
  # file.remove("decoy.html")
})
