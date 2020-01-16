#' Radar (spider) plot for gene sets
#'
#' Radar (spider) plot for gene sets, either for one or more results from functional
#' enrichment analysis.
#'
#' @param res_enrich1 A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_enrich2 TODO
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed
#' @param p_value_column Character string, specifying the column of `res_enrich`
#' where the p-value to be represented is specified. Defaults to `gs_pvalue`
#' (it could have other values, in case more than one p-value - or an adjusted
#' p-value - have been specified).
#'
#' @return A `plotly` object
#' @export
#'
#' @examples
#'
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'                      keys = rownames(dds_macrophage),
#'                      column = "SYMBOL",
#'                      keytype = "ENSEMBL"),
#'   stringsAsFactors = FALSE,
#'   row.names = rownames(dds_macrophage)
#' )
#'
#' # res object
#' data(res_de_macrophage, package = "GeneTonic")
#' res_de <- res_macrophage_IFNg_vs_naive
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#' res_enrich <- get_aggrscores(res_enrich, res_de, anno_df)
#' gs_radar(res_enrich1 = res_enrich)
#' # or using the alias...
#' gs_spider(res_enrich1 = res_enrich)
gs_radar <- function(res_enrich1,
                     res_enrich2,
                     n_gs = 20,
                     p_value_column = "gs_pvalue") {
  # TODO: option to have more and compare them on the same plot!

  # TODO: option to have a list

  # TODO: check that the categories contained are the same
  # -- or option to merge w NAs or so

  # res_enrich has to contain the Z-score to be displayed

  if (!("z_score" %in% colnames(res_enrich1))) {
    warning("You need to add the z_score or the aggregated score")
  } # TODO: same for aggr_score

  # TODO: will remove
  set.seed(42)
  shuffled_ones <- sample(seq_len(n_gs))
  # TODO

  res_enrich1$logp10 <- -log10(res_enrich1[[p_value_column]])

  res_enrich1 <- res_enrich1[seq_len(n_gs), ]
  res_enrich2 <- res_enrich1[seq_len(n_gs), ]

  res_enrich2$logp10 <- -log10(res_enrich2[[p_value_column]])

  res_enrich2$logp10 <- res_enrich2$logp10[shuffled_ones]
  res_enrich2$z_score <- res_enrich2$z_score[shuffled_ones]

  log_smallest_p <- max(res_enrich1$logp10, res_enrich2$logp10)
  set_colors <- RColorBrewer::brewer.pal(n = 8, "Set1")

  p <- plot_ly(
    type = "scatterpolar",
    mode = "markers",
    fill = "toself"
  ) %>%
    add_trace(
      r = c(res_enrich1$logp10, res_enrich1$logp10[1]), # recycling the first element
      theta = c(res_enrich1[["gs_description"]], res_enrich1[["gs_description"]][1]),
      name = "scenario 1"
    ) %>%
    add_trace(
      r = c(res_enrich2$logp10, res_enrich2$logp10[1]),
      theta = c(res_enrich2[["gs_description"]], res_enrich2[["gs_description"]][1]),
      name = "scenario 2"
    ) %>%
    plotly::layout(
      polar = list(radialaxis = list(visible = TRUE,
                                     range = c(0, log_smallest_p))
                   )
      # ,
      # title = "Geneset Radar Chart", font = list(size = 10)
    )

  # TODO: do some kind of for loop/list approach?

  # TODO: some better tooltips
  # ideas: https://www.r-bloggers.com/radar-charts-in-r-using-plotly/

  return(p)
}

# TODO: do it like https://stackoverflow.com/questions/37670412/can-i-recreate-this-polar-coordinate-spider-chart-in-plotly/37778091#37778091   !!

#' @rdname gs_radar
#' @export
gs_spider <- gs_radar
