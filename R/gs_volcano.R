#' Volcano plot for gene sets
#'
#' Volcano plot for gene sets, to summarize visually the functional enrichment
#' results
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' This object needs to be processed first by a function such as [get_aggrscores()]
#' to compute the term-wise `z_score` or `aggr_score`, which will be used for plotting
#' @param p_threshold Numeric, defines the threshold to be used for filtering the
#' gene sets to display. Defaults to 0.05
#' @param color_by Character specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults to `aggr_score`.
#' @param volcano_labels Integer, maximum number of labels for the gene sets to be
#' plotted as labels on the volcano scatter plot.
#' @param scale_circles TODO
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be labeled.
#' @param plot_title Character string, used as title for the plot. If left `NULL`,
#' it defaults to a general description of the plot and of the DE contrast
#'
#' @return A `ggplot` object
#'
#' @seealso [gs_simplify()] can be applied in advance to `res_enrich` to reduce
#' the redundancy of the displayed gene sets
#'
#' @export
#'
#' @examples
#' #  TODO
gs_volcano <- function(res_enrich,
                       p_threshold = 0.05,
                       color_by = "aggr_score",
                       volcano_labels = 10,
                       scale_circles = 1, # TODOTODO: see how to control point size
                       # TODO option to collapse similar terms?
                       gs_ids = NULL,
                       plot_title = NULL
) {
  # res_enrich has to contain the aggregated scores
  if (!all(c("z_score", "aggr_score") %in% colnames(res_enrich)))
    stop("You might need to compute the aggregated scores first")
  # TODO: or call in advance the get_aggr_scores function?

  volcano_labels <- min(volcano_labels, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(volcano_labels)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )


  volcano_df <- res_enrich
  volcano_df$logpval <- -log10(volcano_df[["gs_pvalue"]])
  volcano_df$gs_name <- volcano_df[["gs_description"]]
  volcano_df$`set members` <- volcano_df[["gs_de_count"]]

  volcano_df <- volcano_df[volcano_df[["gs_pvalue"]] <= p_threshold, ]
  max_z <- max(abs(range(volcano_df[[color_by]])))
  limit <- max_z * c(-1, 1)

  p <- ggplot(
    volcano_df,
    aes_string(x = "z_score", y = "logpval", size = "`set members`",  text = "gs_name")) +
    geom_point(aes_string(col = color_by), shape = 20, alpha = 1) +
    labs(x = "geneset Z score",
         y= "log10 p-value",
         size = "Gene set\nmembers",
         col = "Aggregated\nscore") +
    scale_x_continuous(limits = limit) +
    theme_bw() +
    scale_color_gradient2(limit = limit,
                          low = muted("deepskyblue"), high = muted("firebrick"), mid = "lightyellow")

  if (length(gs_to_use > 0)) {
    df_gs_labels <- volcano_df[volcano_df$gs_id %in% gs_to_use, ]

    p <- p + geom_label_repel(
      aes_string(label = "gs_name"), data = df_gs_labels, size = 4)
  }

  # handling the title
  if (is.null(plot_title)) {
    p <- p + ggtitle(paste0("Geneset volcano"))
  } else {
    p <- p + ggtitle(plot_title)
  }
  # beautify all I have in the plot...

  # something TODOTODO to prepare the labels for plotly, if I am using it
  ## one option: pass a text aes to geom_point, which is not used in gg but in plotly

  return(p)

}
