# something on the line of plotCounts, ggplotCounts, but with more pimpedity :D
## maybe even plotly-fied already, or pimped in gg so that it is readily plugged into ggplotly

#' Title TODO
#'
#' TODO
#'
#' @param dds TODO
#' @param gene TODO
#' @param assay TODO
#' @param intgroup TODO
#' @param annotation_obj TODO
#' @param transform TODO
#' @param labels_repel TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gene_plot <- function(dds,
                      gene,
                      assay = "counts",
                      intgroup = "condition",
                      annotation_obj = NULL,
                      transform = TRUE,
                      labels_repel = TRUE) {

  # TODO: enable the possibility of using other assays?

  # TODO: use a more general func to extract the values?
  df <- plotCounts(dds, gene, intgroup, returnData = TRUE)

  df$sampleID <- rownames(df)
  if (!is.null(annotation_obj)) {
    genesymbol <- annotation_obj$gene_name[match(gene, annotation_obj$gene_id)]
  } else {
    genesymbol <- ""
  }

  onlyfactors <- df[, match(intgroup, colnames(df))]
  df$plotby <- interaction(onlyfactors)

  min_by_groups <- min(table(df$plotby))
  # depending on this, use boxplots/nothing/violins/sina

  p <- ggplot(df, aes_string(x = "plotby", y = "count", col = "plotby")) +
    scale_x_discrete(name = "") +
    scale_color_discrete(name = "Experimental\ngroup") +
    theme_bw()

  # somewhat following the recommendations here
  # https://www.embopress.org/doi/full/10.15252/embj.201694659
  if (min_by_groups <= 3) {
    p <- p +
      geom_jitter(aes_string(x = "plotby", y = "count"),
                     position = position_jitter(width = 0.2, height = 0))
    # do nothing - or add a line for the median?
  } else if (min_by_groups > 3 & min_by_groups < 10) {
    p <- p +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(width = 0.2, height = 0))

  } else if (min_by_groups >= 11 & min_by_groups < 40) {
    p <- p +
      geom_violin() +
      geom_jitter(position = position_jitter(width = 0.2, height = 0)) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar", width = 0.3)
  } else {
    p <- p +
      ggforce::geom_sina() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar", width = 0.3)
  }

  # handling the labels
  if (labels_repel) {
    p <- p + ggrepel::geom_text_repel(aes_string(label = "sampleID"))
  }
  else {
    p <- p + geom_text(aes_string(label = "sampleID"), hjust = -0.1, vjust = 0.1)
  }

  # handling y axis transformation
  if (transform) {
    p <- p + scale_y_log10(name = "Normalized counts (log10 scale)")
  } else {
    p <- p + scale_y_continuous(name = "Normalized counts")
  }

  # handling the displayed names and ids
  if (!is.null(annotation_obj)) {
    p <- p + labs(title = paste0("Normalized counts for ", genesymbol, " - ", gene))
  } else {
    p <- p + labs(title = paste0("Normalized counts for ", gene))
  }

  return(p)
}
