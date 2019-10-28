#' Title TODO
#'
#' @param res_enrich TODO
#' @param labels_to_use TODO
#' @param pvals_to_use TODO
#' @param p_threshold TODO
#' @param max_nr_labels TODO
#' @param scale_circles TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' #  TODO
gs_volcano <- function(res_enrich,
                       labels_to_use = "Term",
                       pvals_to_use = "p.value_elim",
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
  mydf$logpval <- -log10(mydf[[pvals_to_use]])
  mydf$mylabels <- mydf[[labels_to_use]]
  mydf$`set members` <- mydf$Significant

  mydf <- mydf[mydf[[pvals_to_use]] <= p_threshold, ]
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

  if (!is.null(labels_to_use)) {
    p <- p + geom_label_repel(aes_string(label = "mylabels"), data = mydf[1:max_nr_labels, ], size = 4)
  }

  # handling the title
  p <- p + ggtitle("TODOTODO some title related to the provided info")
  # and maybe a caption as well

  # beautify all I have in the plot...

  # something TODOTODO to prepare the labels for plotly, if I am using it
  ## one option: pass a text aes to geom_point, which is not used in gg but in plotly

  return(p)

}
