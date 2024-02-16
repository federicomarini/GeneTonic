# something on the line of plotCounts, ggplotCounts, but with more pimpedity :D
## maybe even plotly-fied already, or pimped in gg so that it is readily plugged into ggplotly

#' Plot expression values for a gene
#'
#' Plot expression values (e.g. normalized counts) for a gene of interest, grouped
#' by experimental group(s) of interest
#'
#' @details The result of this function can be fed directly to [plotly::ggplotly()]
#' for interactive visualization, instead of the static `ggplot` viz.
#'
#' @param dds A `DESeqDataSet` object, normally obtained after running your data
#' through the `DESeq2` framework.
#' @param gene Character, specifies the identifier of the feature (gene) to be
#' plotted
#' @param intgroup A character vector of names in `colData(dds)` to use for grouping.
#' Note: the vector components should be categorical variables.
#' @param assay Character, specifies with assay of the `dds` object to use for
#' reading out the expression values. Defaults to "counts".
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param normalized Logical value, whether the expression values should be
#' normalized by their size factor. Defaults to TRUE, applies when `assay` is
#' "counts"
#' @param transform Logical value, corresponding whether to have log scale y-axis
#' or not. Defaults to TRUE.
#' @param labels_display Logical value. Whether to display the labels of samples,
#' defaults to TRUE.
#' @param labels_repel Logical value. Whether to use `ggrepel`'s functions to
#' place labels; defaults to TRUE
#' @param plot_type Character, one of "auto", "jitteronly", "boxplot", "violin",
#' or "sina". Defines the type of `geom_` to be used for plotting. Defaults to
#' `auto`, which in turn chooses one of the layers according to the number of
#' samples in the smallest group defined via `intgroup`
#' @param return_data Logical, whether the function should just return the
#' data.frame of expression values and covariates for custom plotting. Defaults
#' to FALSE.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#'
#' @return A `ggplot` object
#' @export
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'     keys = rownames(dds_macrophage),
#'     column = "SYMBOL",
#'     keytype = "ENSEMBL"
#'   ),
#'   stringsAsFactors = FALSE,
#'   row.names = rownames(dds_macrophage)
#' )
#'
#' gene_plot(dds_macrophage,
#'   gene = "ENSG00000125347",
#'   intgroup = "condition",
#'   annotation_obj = anno_df
#' )
gene_plot <- function(dds,
                      gene,
                      intgroup = "condition",
                      assay = "counts",
                      annotation_obj = NULL,
                      normalized = TRUE,
                      transform = TRUE,
                      labels_display = TRUE,
                      labels_repel = TRUE,
                      plot_type = "auto",
                      return_data = FALSE,
                      gtl = NULL) {
  plot_type <- match.arg(
    plot_type,
    c("auto", "jitteronly", "boxplot", "violin", "sina")
  )

  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }
  
  if (!intgroup %in% colnames(colData(dds))) {
    stop("`intgroup` not found in the colData slot of the dds object",
         "\nPlease specify one of the following: \n",
         paste0(colnames(colData(dds)), collapse = ", ")
         )
  }

  df <- get_expression_values(
    dds = dds,
    gene = gene,
    intgroup = intgroup,
    assay = assay,
    normalized = normalized
  )

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

  if (return_data) {
    return(df)
  }

  p <- ggplot(df, aes(x = .data$plotby, 
                      y = .data$exp_value, 
                      col = .data$plotby)) +
    scale_x_discrete(name = "") +
    scale_color_discrete(name = "Experimental\ngroup") +
    theme_bw()

  # for connected handling of jittered points AND labels
  jit_pos <- position_jitter(width = 0.2, height = 0, seed = 42)

  # somewhat following the recommendations here
  # https://www.embopress.org/doi/full/10.15252/embj.201694659
  if (plot_type == "jitteronly" || (plot_type == "auto" & min_by_groups <= 3)) {
    p <- p +
      geom_point(aes(x = .data$plotby, 
                     y = .data$exp_value),
        position = jit_pos
      )
    # do nothing - or add a line for the median?
  } else if (plot_type == "boxplot" || (plot_type == "auto" & (min_by_groups > 3 & min_by_groups < 10))) {
    p <- p +
      geom_boxplot(outlier.shape = NA) +
      geom_point(aes(x = .data$plotby, 
                     y = .data$exp_value), position = jit_pos)
  } else if (plot_type == "violin" || (plot_type == "auto" & (min_by_groups >= 11 & min_by_groups < 40))) {
    p <- p +
      geom_violin() +
      geom_point(aes(x = .data$plotby, 
                     y = .data$exp_value), position = jit_pos) +
      stat_summary(
        fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", width = 0.3
      )
  } else if (plot_type == "sina" || (plot_type == "auto" & (min_by_groups >= 40))) {
    p <- p +
      ggforce::geom_sina() +
      stat_summary(
        fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", width = 0.3
      )
  }

  # handling the labels
  if (labels_display) {
    if (labels_repel) {
      p <- p + ggrepel::geom_text_repel(aes(label = .data$sample_id),
        min.segment.length = 0,
        position = jit_pos
      )
    }
    else {
      p <- p + geom_text(aes(label = .data$sample_id),
        hjust = -0.1, vjust = 0.1,
        position = jit_pos
      )
    }
  }

  y_label <- if (assay == "counts" & normalized) {
    "Normalized counts"
  } else if (assay == "counts" & !normalized) {
    "Counts"
  } else if (assay == "abundance") {
    "TPM - Transcripts Per Million"
  } else {
    assay
  }

  # handling y axis transformation
  if (transform) {
    p <- p + scale_y_log10(name = paste0(y_label, " (log10 scale)"))
  } else {
    p <- p + scale_y_continuous(name = y_label)
  }

  # handling the displayed names and ids
  if (!is.null(annotation_obj)) {
    p <- p + labs(title = paste0(genesymbol, " - ", gene))
  } else {
    p <- p + labs(title = paste0(gene))
  }

  return(p)
}

#' Get expression values
#'
#' Extract expression values, with the possibility to select other assay slots
#'
#' @param dds A `DESeqDataSet` object, normally obtained after running your data
#' through the `DESeq2` framework.
#' @param gene Character, specifies the identifier of the feature (gene) to be
#' extracted
#' @param intgroup A character vector of names in `colData(dds)` to use for grouping.
#' @param assay Character, specifies with assay of the `dds` object to use for
#' reading out the expression values. Defaults to "counts".
#' @param normalized Logical value, whether the expression values should be
#' normalized by their size factor. Defaults to TRUE, applies when `assay` is
#' "counts"
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#'
#' @return A tidy data.frame with the expression values and covariates for further
#' processing
#' @export
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' df_exp <- get_expression_values(dds_macrophage,
#'   gene = "ENSG00000125347",
#'   intgroup = "condition"
#' )
#' head(df_exp)
get_expression_values <- function(dds,
                                  gene,
                                  intgroup,
                                  assay = "counts",
                                  normalized = TRUE,
                                  gtl = NULL) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
  }
  
  if (!(assay %in% names(assays(dds)))) {
    stop(
      "Please specify a name of one of the existing assays: \n",
      paste(names(assays(dds)), collapse = ", ")
    )
  }

  # checking the normalization factors are in
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }

  if (assay == "counts") {
    exp_vec <- counts(dds, normalized = normalized)[gene, ]
  } else {
    exp_vec <- assays(dds)[[assay]][gene, ]
  }

  exp_df <- data.frame(
    exp_value = exp_vec,
    colData(dds)[intgroup]
  )

  return(exp_df)
}
