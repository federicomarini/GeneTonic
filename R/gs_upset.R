#' Create a geneset upset dataset
#' 
#' Create a data frame that can be fed to the upset function
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to see the
#' formatting requirements.
#' @param use_ids Logical - whether to use the `gs_id`entifiers as names, or the 
#' values provided as `gs_description`. Defaults to FALSE, using the full
#' descriptions
#'
#' @return A data.frame to be used in `ComplexUpset::upset()`
#' @export
#'
#' @examples
#' # TODO
create_upsetdata <- function(res_enrich, use_ids = FALSE) {
  my_gs_ids <- res_enrich$gs_id
  my_gs_desc <- res_enrich$gs_description
  my_gs_genes <- res_enrich$gs_genes
  
  gsg_list <- strsplit(my_gs_genes, ",")
  if (use_ids) {
    names(gsg_list) <- my_gs_ids
  } else {
    names(gsg_list) <- my_gs_desc
  }
  
  all_genes <- unique(unlist(gsg_list))
  upgsg <- unlist(lapply(gsg_list, function(x) {
    x <- as.vector(match(all_genes, x))
  }))
  
  upgsg[is.na(upgsg)] <- as.integer(0)
  upgsg[upgsg != 0] <- as.integer(1)
  
  upgsg <- data.frame(
    matrix(upgsg, ncol = length(gsg_list), byrow = FALSE)
  )
  
  upgsg <- upgsg[which(rowSums(upgsg) != 0), ]
  names(upgsg) <- names(gsg_list)
  rownames(upgsg) <- all_genes
  
  return(upgsg)
}




#' Upset plot for genesets
#' 
#' Create an upset plot for genesets
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be included in addition to
#' the top ones (via `n_gs`)
#' @param add_de_direction Logical, whether to add an annotation with info on the 
#' DE direction of single genes
#' @param add_de_gsgenes Logical, if set to TRUE adds an annotation with detail
#' on the single components of each defined subset
#' @param col_upDE Character, specifying the color value to be used to mark
#' upregulated genes
#' @param col_downDE Character, specifying the color value to be used to mark
#' downregulated genes
#' @param upset_geom A geom specification to be used in the upset chart. Defaults
#' sensibly to `geom_point(size = 2)`
#' @param return_upsetgsg Logical, controlling the returned value. If set to TRUE,
#' this function will not generate the plot but only create the corresponding 
#' data.frame, in case the user wants to proceed with a custom call to create an
#' upset plot.
#'
#' @return A `ggplot` object (if plotting), or alternatively a data.frame
#' @export
#'
#' @examples
#' # TODO
gs_upset <- function(res_enrich,
                     res_de,
                     annotation_obj = NULL,
                     n_gs = 10,
                     gtl = NULL,
                     gs_ids = NULL,
                     add_de_direction = FALSE,
                     add_de_gsgenes = FALSE,
                     col_upDE = "#E41A1C",
                     col_downDE = "#377EB8",
                     upset_geom = geom_point(size = 2),
                     return_upsetgsg = FALSE) {
  
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    # dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }
  
  stopifnot(is.numeric(n_gs))
  
  n_gs <- min(n_gs, nrow(res_enrich))
  
  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)], # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id] # the ones specified from the custom list
    )
  )
  
  re <- res_enrich[gs_to_use, ]
  
  upgsg <- create_upsetdata(re)
  
  if(add_de_direction) { # TODO extra checks, need DE and annot
    upgsg$logFC <- res_de[annotation_obj$gene_id[
      match(rownames(upgsg), annotation_obj$gene_name)], "log2FoldChange"]
    upgsg$logFCsign <- upgsg$logFC >= 0
    
    param_upset_baseanno <- list(
      'Intersection size' = intersection_size(
        counts = FALSE,
        mapping = aes_string(fill = "logFCsign")
      ) + 
        scale_fill_manual(
          values=c('FALSE' = col_downDE, 'TRUE' = col_upDE),
          labels = c("logFC<0", "logFC>0"),
          name = "") + 
        theme(legend.position="bottom")
    )
  } else {
    param_upset_baseanno <- "auto"
  }
  
  if (add_de_gsgenes){
    param_upset_anno <- list(
      'logFC' = (
        ggplot(mapping = aes_string(x = "intersection", y = "logFC")) +  
          geom_jitter(aes_string(color = "logFC"), na.rm = TRUE) + 
          # geom_violin(alpha = 0.5, na.rm = TRUE) + 
          theme(legend.position = "none") + 
          scale_colour_gradient2(low = col_downDE, high = col_upDE)
      )
    )
  } else {
    param_upset_anno <- NULL
  }
  
  # if desired to work on this separately
  if (return_upsetgsg) {
    return(upgsg)
  }

  # ready to plot
  ComplexUpset::upset(
    data = upgsg,
    intersect = names(upgsg)[seq_len(length(gs_to_use))],
    name = "gsg upset",
    matrix = intersection_matrix(
      geom = upset_geom
    ) + scale_y_discrete(position="right"),
    base_annotations = param_upset_baseanno,
    annotations = param_upset_anno,
    width_ratio=0.1
  )
  
}


