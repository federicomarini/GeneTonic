#' Volcano plot for gene sets
#'
#' Volcano plot for gene sets, to summarize visually the functional enrichment
#' results
#'
#' It is also possible to reduce the redundancy of the input `res_enrich` object,
#' if it is passed in advance to the [gs_simplify()] function.
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
#' @param scale_circles A numeric value, to define the scaling factor for the
#' circle sizes. Defaults to 1.
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
#'
#' gs_volcano(res_enrich)
#'
gs_volcano <- function(res_enrich,
                       p_threshold = 0.05,
                       color_by = "aggr_score",
                       volcano_labels = 10,
                       scale_circles = 1,
                       gs_ids = NULL,
                       plot_title = NULL
) {
  # res_enrich has to contain the aggregated scores
  if (!all(c("z_score", "aggr_score") %in% colnames(res_enrich)))
    stop("You might need to compute the aggregated scores first")

  if (!color_by %in% colnames(res_enrich))
    stop("Your res_enrich object does not contain the ",
         color_by,
         " column.\n",
         "Compute this first or select another column to use for the color.")

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
  volcano_df$`set members` <- volcano_df[["gs_de_count"]] * scale_circles

  volcano_df <- volcano_df[volcano_df[["gs_pvalue"]] <= p_threshold, ]
  max_x <- max(abs(range(volcano_df[["z_score"]])))
  max_z <- max(abs(range(volcano_df[[color_by]])))
  limit_x <- max_x * c(-1, 1)
  limit_z <- max_z * c(-1, 1)

  p <- ggplot(
    volcano_df,
    aes_string(x = "z_score", y = "logpval", size = "`set members`",  text = "gs_name")) +
    geom_point(aes_string(col = color_by), shape = 20, alpha = 1) +
    labs(x = "geneset Z score",
         y = "-log10 p-value",
         size = "Gene set\nmembers",
         col = "Aggregated\nscore") +
    scale_x_continuous(limits = limit_x) +
    theme_bw() +
    scale_color_gradient2(limit = limit_z,
                          low = muted("deepskyblue"), high = muted("firebrick"), mid = "lightyellow")

  if (length(gs_to_use > 0)) {
    df_gs_labels <- volcano_df[volcano_df$gs_id %in% gs_to_use, ]

    p <- p + geom_label_repel(
      aes_string(label = "gs_name"), 
      data = df_gs_labels, 
      size = 4,
      min.segment.length = 0
    )
  }

  # handling the title
  if (is.null(plot_title)) {
    p <- p + ggtitle(paste0("Geneset volcano"))
  } else {
    p <- p + ggtitle(plot_title)
  }

  return(p)
}
