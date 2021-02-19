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
#'  Default color is #1a81c2.
#' @param volcano_labels Integer, maximum number of labels for the gene sets to be
#' plotted as labels on the volcano scatter plot. Default value is 25.
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
#' ggs_volcano(res_de,
#'            res_enrich,
#'            anno_df,
#'            geneset_id = res_enrich$gs_id[1]
#'            )
#'            
#' alternatively
#'
#' # genelist <- sapply(anno_df$gene_name, function(x) str_contains(x, "CD"))
#' 
#' # ggs_volcano(res_de,
#' #             res_enrich,
#' #             anno_df,
#' #             genelist = genelist
#' #             )
#' 
ggs_volcano <- function(res_de,
                       res_enrich,
                       annotation_obj = NULL,
                       gtl = NULL,
                       geneset_id = NULL,
                       genelist = NULL,
                       FDR = 0.05,
                       color = "red",
                       volcano_labels = 25,
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
    # overwritable via a list
    if (!all(genelist %in% rownames(res_de))) {
      not_there <- genelist[!(genelist %in% rownames(res_de))]
      warning("Some of the provided gene ids were not found in the SummarizedExperiment",
              "\nNot found: ", not_there)
    }
    thisset_members_ids <- intersect(rownames(res_de), genelist)
    thisset_name <- "Custom list"
  }

  thisset_members_ids
  
  # Prepare the data
  complete_genes_ids <- rownames(res_de)
  complete_genes <- annotation_obj$gene_name[match(complete_genes_ids, annotation_obj$gene_id)]
  
  padj_complete <- res_de[complete_genes_ids, "padj"]
  filter_info_complete <- sapply(padj_complete, function(x) x <= FDR)
  padj_complete <- sapply(padj_complete, function(x) -log10(x))
  
  log2FoldChange_complete <- res_de[complete_genes_ids, "log2FoldChange"]
  
  gene_set_belong <- complete_genes_ids %in% thisset_members_ids
  filter_info_complete <- filter_info_complete & gene_set_belong

  
  thisset_complete_data <- data.frame(complete_genes_ids, padj_complete, log2FoldChange_complete, filter_info_complete, gene_set_belong)
  colnames(thisset_complete_data) <- c("genes", "logTransformedpvalue", "log2FoldChange", "significant", "belonging")
  
  
  # Prepare plotting
  volcano_df_complete <- thisset_complete_data
  volcano_df_complete$genes_name <- complete_genes
  max_x <- max(abs(range(thisset_complete_data["log2FoldChange"])))
  limit_x <- max_x * c(-1, 1)
  

  # Prepare plot title
  if (is.null(plot_title)) {
    title <- paste0("Volcano Plot - ", thisset_name, " - ", geneset_id)
  } else {
    title <- plot_title
  }


  # Plot data
  p <- ggplot(volcano_df_complete, aes_string(x = "log2FoldChange", y = "logTransformedpvalue")) +
       geom_point(aes_string(color = "significant", alpha = "belonging")) +
       labs(x = "log2Fold Change",
            y = "-log10 p-value",
            color = "p-value <= FDR") + 
       scale_x_continuous(limits = limit_x) +
       scale_color_manual(labels = c("significant", "not significant"), 
                          breaks = c("TRUE", "FALSE"),
                          values = c(color, "grey25")) +
       scale_alpha_manual(breaks = c("TRUE", "FALSE"),
                          values = c(1, 1/10)) +
       theme_bw() +
       theme(legend.title = element_text(size = 11, face = "bold"),
             legend.text = element_text(size = 10))  
    
  
  # adding labels to the significant points of the geneset
  p <- p + geom_text_repel(
    data = subset(volcano_df_complete, filter_info_complete),
    aes_string(label = "genes_name"),
    size = 4,
    max.overlaps = volcano_labels)

  
  # handling the title
  p <- p + ggtitle(title)
  p <- p + guides(alpha = FALSE)
  return(p)
}
