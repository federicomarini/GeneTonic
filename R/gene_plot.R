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
#' Note: the vector components should be categorical variables. Defaults
#' to NULL, which which would then select the first column of the `colData` slot.
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
                      intgroup = NULL,
                      assay = "counts",
                      annotation_obj = NULL,
                      normalized = TRUE,
                      transform = TRUE,
                      labels_display = TRUE,
                      labels_repel = TRUE,
                      plot_type = "auto",
                      return_data = FALSE,
                      gtl = NULL) {
  .Deprecated(old = "gene_plot", new = "mosdef::gene_plot", 
              msg = paste0(
                "Please use `mosdef::gene_plot()` in replacement of the `gene_plot()` function, ",
                "originally located in the GeneTonic package. \nCheck the manual page for ",
                "`?mosdef::gene_plot()` to see the details on how to use it"))
  
  p <- mosdef::gene_plot(de_container = dds,
                         gene = gene, 
                         intgroup = intgroup,
                         assay = assay,
                         annotation_obj = annotation_obj,
                         normalized = normalized,
                         transform = transform,
                         labels_display = labels_display,
                         labels_repel = labels_repel,
                         plot_type = plot_type,
                         return_data = return_data)
  
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
  
  .Deprecated(old = "get_expression_values", new = "mosdef::get_expr_values", 
              msg = paste0(
                "Please use `mosdef::get_expr_values()` in replacement of the `get_expression_values()` function, ",
                "originally located in the GeneTonic package. \nCheck the manual page for ",
                "`?mosdef::get_expr_values()` to see the details on how to use it"))
  
  expr_values <- mosdef::get_expr_values(
    de_container = dds,
    gene = gene, 
    intgroup = intgroup, 
    assay = assay, 
    normalized = normalized
  )
  
  return(expr_values)
}
