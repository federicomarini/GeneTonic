#' Plot a volcano plot of the gene signature on the data
#'
#' Plot a volcano plot for the selected gene signature on the provided data, with the possibility to compactly display also DE only genes
#'
#' @param res_de A `DESeqResults` object.
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param geneset_id Character specifying the gene set identifier to be plotted
#' @param genelist A vector of character strings, specifying the identifiers
#' contained in the row names of the `se` input object.
#' @param FDR Numeric value, specifying the significance level for thresholding
#' adjusted p-values. Defaults to 0.05.
#' @param color Character string to specify color of filtered points in the plot.
#'  Default color is red.
#' @param volcano_labels Integer, maximum number of labels for the gene sets to be
#' plotted as labels on the volcano scatter plot.
#' @param plot_title Character string, to specify the title of the plot,
#' displayed over the heatmap. If left to `NULL` as by default, it tries to use
#' the information on the geneset identifier provided
#'
#' @return A plot returned by the [ggplot()] function
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
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' vst_macrophage <- vst(dds_macrophage)
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
#' ggs_volcano(vst_macrophage,
#'            res_de,
#'            res_enrich,
#'            anno_df,
#'            geneset_id = res_enrich$gs_id[1]
#'            )
#'
ggs_volcano <- function(res_de,
                       res_enrich,
                       annotation_obj = NULL,
                       gtl = NULL,
                       geneset_id = NULL,
                       genelist = NULL,
                       FDR = 0.05,
                       color = "red",
                       volcano_labels = 10,
                       plot_title = NULL
                       ) {

  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  # Retrieve information about genes in geneset gs_id
  if (!is.null(geneset_id)) {
    if (geneset_id %in% res_enrich[["gs_id"]]) {
      thisset_name <- res_enrich[geneset_id, "gs_description"]
      thisset_members <- unlist(strsplit(res_enrich[geneset_id, "gs_genes"], ","))
      thisset_members_ids <- annotation_obj$gene_id[match(thisset_members, annotation_obj$gene_name)]
    }
  }else {
    # overridable via a list
    if (!all(genelist %in% rownames(se))) {
      not_there <- genelist[!(genelist %in% rownames(se))]
      warning("Some of the provided gene ids were not found in the SummarizedExperiment",
              "\nNot found: ", not_there)
    }
    thisset_members_ids <- intersect(rownames(se), genelist)
    thisset_name <- "Custom list"
  }

  thisset_members_ids
  thisset_members
  
  # retrieve adjusted pvalues and log2Fold Change for genes
  padj <- res_de[thisset_members_ids, "padj"]
  padjlog <- sapply(padj, function(x) -log10(x))
  log2FoldChange <- res_de[thisset_members_ids, "log2FoldChange"]
  filter_info <- sapply(padj, function(x) x <= FDR)
  
  # Prepare data to plot
  thisset_data <- data.frame(thisset_members_ids, padj, log2FoldChange, padjlog, filter_info)

  
  if (is.null(plot_title)) {
    title <- paste0("Volcano Plot - ", thisset_name, " - ", geneset_id)
  } else {
    title <- plot_title
  }
  
  colnames(thisset_data) <- c("genes", "pvalue", "log2FoldChange", "logTransformedpvalue", "significant")
  
  
  # Prepare plotting
  volcano_df <- thisset_data
  volcano_df$genes_name <- thisset_members
  max_x <- max(abs(range(thisset_data["log2FoldChange"])))
  limit_x <- max_x * c(-1, 1)


  # Plot data
   p <- ggplot(volcano_df, aes(x = log2FoldChange, y = logTransformedpvalue, label = genes_name)) +
    geom_point(aes(color = significant)) +
    labs(x = "log2Fold Change",
         y = "-log10 p-value", color = "p-value <= FDR") + 
    scale_x_continuous(limits = limit_x) +
    scale_color_manual(labels = c("significant", "not significant"), 
                       breaks = c("TRUE", "FALSE"), values = c(color, "grey25")) +
    theme(legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 10)) +
    geom_text_repel()

    
  
  # handling the title
  p <- p + ggtitle(title)
  return(p)
}
