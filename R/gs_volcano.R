#' Volcano plot for gene sets
#'
#' Volcano plot for gene sets, to summarize visually the functional enrichment
#' results
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to see the
#' formatting requirements. This object needs to be processed first by a function
#' such as [get_aggrscores()] to compute the term-wise `z_score` or `aggr_score`,
#' which will be used for plotting
#' @param p_threshold Numeric, defines the threshold to be used for filtering the
#' gene sets to display. Defaults to 0.05
#' @param max_nr_labels Integer, maximum number of labels for the gene sets to be
#' plotted as labels on the volcano scatter plot.
#' @param scale_circles TODO
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
                       max_nr_labels = 10,
                       scale_circles = 1 # TODOTODO: see how to control point size
# TODO option to collapse similar terms?
                       ) {
  # res_enrich has to contain the aggregated scores
  if (!all(c("z_score", "aggr_score") %in% colnames(res_enrich)))
    stop("You might need to compute the aggregated scores first")
  # TODO: or call in advance the get_aggr_scores function?

  mydf <- res_enrich
  mydf$logpval <- -log10(mydf[["gs_pvalue"]])
  mydf$mylabels <- mydf[["gs_description"]]
  mydf$`set members` <- mydf[["gs_de_count"]]

  mydf <- mydf[mydf[["gs_pvalue"]] <= p_threshold, ]
  max_z <- max(abs(range(mydf$z_score)))
  limit <- max_z * c(-1, 1)

  p <- ggplot(
    mydf,
    aes_string(x = "z_score", y = "logpval", size = "`set members`",  text = "mylabels")) +
    # geom_point(aes(col = aggr_score),shape = 20, alpha = 1) +
    geom_point(aes_string(col = "aggr_score"), shape = 20, alpha = 1) +
    scale_x_continuous(limits = limit) +
    theme_bw() +
    scale_color_gradient2(limit = limit,
                          low = muted("deepskyblue"), high = muted("firebrick"), mid = "lightyellow")

  # TODO: condition this?
  p <- p + geom_label_repel(aes_string(label = "mylabels"), data = mydf[1:max_nr_labels, ], size = 4)

  # handling the title
  p <- p + ggtitle("TODOTODO some title related to the provided info")
  # and maybe a caption as well

  # beautify all I have in the plot...

  # something TODOTODO to prepare the labels for plotly, if I am using it
  ## one option: pass a text aes to geom_point, which is not used in gg but in plotly

  return(p)

}
