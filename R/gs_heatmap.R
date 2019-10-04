

#' Title
#'
#' TODO
#'
#' @param se TODO
#' @param res_de TODO
#' @param res_enrich TODO
#' @param geneset_id TODO
#' @param genelist TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#' @param annotation_obj TODO
#' @param FDR TODO
#' @param de_only TODO
#' @param cluster_rows TODO
#' @param cluster_cols TODO
#' @param center_mean TODO
#' @param scale_row TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_heatmap <- function(se,
                       res_de,
                       res_enrich,
                       geneset_id,
                       genelist,
                       genes_colname = "genes",
                       genesetname_colname = "Term",
                       genesetid_colname = "GO.ID",
                       annotation_obj = NULL,
                       FDR = 0.05,
                       de_only = FALSE,
                       cluster_rows = TRUE, # TODOTODO: options for the heatmap go on left side, as could be common to more!
                       cluster_cols = FALSE,
                       center_mean = TRUE,
                       scale_row = FALSE
                       # TODOTODO: use ellipsis for passing params to pheatmap?
                       # TODOTODO: option to just return the underlying data?s
                       # TODOTODO: options to subset to specific samples?
                       ) {

  # check that the data would ideally be a DST, so that it is not the counts/normalized?
  mydata <- assay(se)

  # if(geneset in the results)
  #   pick the genes from there
  # else
  #   option to override the geneset by providing a list
  #
  # idea: multiselect with gene names - but in the UI
  # internal matching to the IDs (in this function we use the ids already)

  rownames(res_enrich) <- res_enrich[[genesetid_colname]]
  if (geneset_id %in% res_enrich[[genesetid_colname]]) {
    thisset_name <- res_enrich[geneset_id,genesetname_colname]
    thisset_members <- unlist(strsplit(res_enrich[geneset_id,genes_colname],","))
    thisset_members_ids <- annotation_obj$gene_id[match(thisset_members, annotation_obj$gene_name)]
  } else {
    # overridable via a list
  }

  thisset_members_ids
  sig_to_keep <- (thisset_members_ids %in% rownames(se))#
  thisset_members_ids_available <- thisset_members_ids[sig_to_keep]

  mydata_sig <- mydata[thisset_members_ids_available,]

  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(mydata_sig, 1, var) == 0
  mydata_sig <- mydata_sig[!to_remove,]

  if(center_mean)
    mydata_sig <- mydata_sig - rowMeans(mydata_sig)

  if(de_only) {
    de_res <- ideal::deseqresult2DEgenes(res_de,FDR)
    de_genes <- de_res$id
    de_to_keep <- rownames(mydata_sig) %in% de_genes
    mydata_sig <- mydata_sig[de_to_keep,]
  }

  # dim(mydata_sig)

  title <- paste0("Signature heatmap - ", thisset_name)
  sample_decoration <- as.data.frame(colData(se))[,"condition",drop = FALSE]

  # if(returnData) {
  #
  # }

  pheatmap(mydata_sig,
           # annotation_col = anno_colData,
           cluster_rows = cluster_rows, cluster_cols = cluster_cols,
           scale = ifelse(scale_row,"row","none"),main = title,
           labels_row = annotation_obj[rownames(mydata_sig),]$gene_name,
           annotation_col = sample_decoration)


}
