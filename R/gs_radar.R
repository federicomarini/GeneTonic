#' Radar (spider) plot for gene sets
#'
#' Radar (spider) plot for gene sets, either for one or more results from functional
#' enrichment analysis.
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_enrich2 Analogous to `res_enrich1`, another `data.frame` object,
#' storing the result of the functional enrichment analysis, but for a different
#' setting (e.g. another contrast).
#' Defaults to NULL (in this case, a single set of enrichment results is plotted).
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
#' gs_radar(res_enrich = res_enrich)
#' # or using the alias...
#' gs_spider(res_enrich = res_enrich)
#'
#' # with more than one set
#' res_enrich2 <- res_enrich[1:60, ]
#' set.seed(42)
#' shuffled_ones <- sample(seq_len(60)) # to generate permuted p-values
#' res_enrich2$gs_pvalue <- res_enrich2$gs_pvalue[shuffled_ones]
#' # ideally, I would also permute the z scores and aggregated scores
#' gs_radar(res_enrich = res_enrich,
#'          res_enrich2 = res_enrich2)
gs_radar <- function(res_enrich,
                     res_enrich2 = NULL,
                     n_gs = 20,
                     p_value_column = "gs_pvalue") {

  # res_enrich has to contain the Z-score to be displayed
  if (!("z_score" %in% colnames(res_enrich))) {
    warning("You need to add the z_score or the aggregated score")
  }

  if (!is.null(res_enrich2)) {
    if (!("z_score" %in% colnames(res_enrich))) {
      warning("You need to add the z_score or the aggregated score")
    }
  }

  # only one set
  if (is.null(res_enrich2)) {
    res_enrich$logp10 <- -log10(res_enrich[[p_value_column]])
    res_enrich <- res_enrich[seq_len(n_gs), ]
    log_smallest_p <- max(res_enrich$logp10)
    set_colors <- RColorBrewer::brewer.pal(n = 8, "Set1")

    p <- plot_ly(
      type = "scatterpolar",
      mode = "markers",
      fill = "toself"
    ) %>%
      add_trace(
        r = c(res_enrich$logp10, res_enrich$logp10[1]), # recycling the first element
        theta = c(res_enrich[["gs_description"]], res_enrich[["gs_description"]][1]),
        name = "scenario 1"
      ) %>%
      plotly::layout(
        polar = list(radialaxis = list(visible = TRUE,
                                       range = c(0, log_smallest_p))
        )
        # ,
        # title = "Geneset Radar Chart", font = list(size = 10)
      )
  } else {
    # if res_enrich2 is also provided
    gs_set1 <- res_enrich$gs_id
    gs_set2 <- res_enrich2$gs_id
    gs_common <- intersect(gs_set1, gs_set2)
    # restrict to the top common n_gs
    gs_common <- gs_common[seq_len(min(n_gs, length(gs_common)))]

    if (length(gs_common) == 0) {
      stop("No gene sets have been found in common to the two enrichment results")
    }

    common_re1 <- res_enrich[gs_common, ]
    common_re2 <- res_enrich2[gs_common, ]

    common_re1$logp10 <- -log10(common_re1[[p_value_column]])
    common_re2$logp10 <- -log10(common_re2[[p_value_column]])
    # if needed, I could access Z and aggregated scores

    log_smallest_p <- max(common_re1$logp10, common_re2$logp10)
    set_colors <- RColorBrewer::brewer.pal(n = 8, "Set1")

    p <- plot_ly(
      type = "scatterpolar",
      mode = "markers",
      fill = "toself"
    ) %>%
      add_trace(
        r = c(common_re1$logp10, common_re1$logp10[1]), # recycling the first element
        theta = c(common_re1[["gs_description"]], common_re1[["gs_description"]][1]),
        name = "scenario 1"
      ) %>%
      add_trace(
        r = c(common_re2$logp10, common_re2$logp10[1]),
        theta = c(common_re2[["gs_description"]], common_re2[["gs_description"]][1]),
        name = "scenario 2"
      ) %>%
      plotly::layout(
        polar = list(radialaxis = list(visible = TRUE,
                                       range = c(0, log_smallest_p))
        )
      )
  }

  return(p)
}


#' @rdname gs_radar
#' @export
gs_spider <- gs_radar
