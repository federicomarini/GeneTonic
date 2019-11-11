# something on the line of plotCounts, ggplotCounts, but with more pimpedity :D
## maybe even plotly-fied already, or pimped in gg so that it is readily plugged into ggplotly

#' Plot expression values for a gene
#'
#' Plot expression values (e.g. normalized counts) for a gene of interest, grouped
#' by experimental group(s) of interest
#'
#' @param dds A `DESeqDataSet` object, normally obtained after running your data
#' through the `DESeq2` framework.
#' @param gene Character, specifies the identifier of the feature (gene) to be
#' plotted
#' @param assay TODO
#' @param intgroup A character vector of names in `colData(dds)` to use for grouping.
#' Note: the vector components should be categorical variables.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param transform Logical value, corresponding whether to have log scale y-axis
#' or not. Defaults to TRUE.
#' @param labels_repel Logical value. Whether to use `ggrepel`'s functions to
#' place labels; defaults to TRUE
#' @param plot_type Character, one of "auto", "jitteronly", "boxplot", "violin",
#' or "sina". Defines the type of `geom_` to be used for plotting. Defaults to
#' `auto`, which in turn chooses one of the layers according to the number of
#' samples in the smallest group defined via `intgroup`
#'
#' @return A `ggplot` object
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
                      labels_repel = TRUE,
                      plot_type = "auto") {

  plot_type <- match.arg(plot_type,
                         c("auto", "jitteronly", "boxplot", "violin", "sina"))
  # TODO: enable the possibility of using other assays?

  # TODO: use a more general func to extract the values?
  df <- plotCounts(dds, gene, intgroup, returnData = TRUE)

  df$sample_id <- rownames(df)
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
  if (plot_type == "jitteronly" || (plot_type == "auto" & min_by_groups <= 3)) {
    p <- p +
      geom_jitter(aes_string(x = "plotby", y = "count"),
                     position = position_jitter(width = 0.2, height = 0))
    # do nothing - or add a line for the median?
  } else if (plot_type == "boxplot" || (plot_type == "auto" & (min_by_groups > 3 & min_by_groups < 10))) {
    p <- p +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(width = 0.2, height = 0))

  } else if (plot_type == "violin" || (plot_type == "auto" & (min_by_groups >= 11 & min_by_groups < 40))) {
    p <- p +
      geom_violin() +
      geom_jitter(position = position_jitter(width = 0.2, height = 0)) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar", width = 0.3)
  } else if (plot_type == "sina" || (plot_type == "auto" & (min_by_groups >= 40))) {
    p <- p +
      ggforce::geom_sina() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar", width = 0.3)
  }

  # handling the labels
  if (labels_repel) {
    p <- p + ggrepel::geom_text_repel(aes_string(label = "sample_id"))
  }
  else {
    p <- p + geom_text(aes_string(label = "sample_id"), hjust = -0.1, vjust = 0.1)
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
