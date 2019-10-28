context("Testing extra functions/helpers for GeneTonic")

test_that("Overlap functions work", {
  set1 <- letters[1:10]
  set2 <- letters[1:15]
  set3 <- letters[5:20]
  ol_1_2 <- overlap_coefficient(set1,set2)
  ol_1_2_ji <- overlap_jaccard_index(set1,set2)

  expect_equal(ol_1_2, 1)
  expect_equal(ol_1_2_ji, 2/3)
})

test_that("map2color works", {
  mypal <- rev(scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.4))
  my_vals <- res_macrophage_IFNg_vs_naive$log2FoldChange[1:20]
  m2c <- map2color(x = my_vals, pal = mypal, limits = c(-4, 4))
  m2c_nolimits <- map2color(x = my_vals, pal = mypal)
  # plot(1:20, col = m2c, pch = 20, cex = 5)

  expect_length(m2c, 20)
  expect_length(m2c_nolimits, 20)
})

test_that("footer code is generated", {
  expect_is(footer(),"shiny.tag")
})

test_that("results to data frame conversion works", {
  res_df <- deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 1)
  res_df2 <- deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.05)
  res_df3 <- deseqresult2df(res_macrophage_IFNg_vs_naive)
  expect_is(res_df, "data.frame")
  expect_error(deseqresult2df(res_df))
})


test_that("Retrieving info on GO term", {
  out <- go_2_html("GO:0032729")
  expect_is(out, "character")
  expect_is(out, "html")
  expect_equal(go_2_html("GO:00"), HTML("GeneOntology term not found!"))
})

test_that("'Linking to AmiGO database", {
  out <- .link2amigo("GO:0032729")
  expect_is(out, "character")
})

test_that("'Linking to NCBI database", {
  out <- .link2ncbi("Actb")
  expect_is(out, "character")
})

test_that("'Linking to GeneCards database", {
  out <- .link2genecards("Gapdh")
  expect_is(out, "character")
})

