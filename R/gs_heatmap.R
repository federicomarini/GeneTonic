

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






#' Title
#'
#' TODO
#'
#' @param se TODO
#' @param res_de TODO
#' @param res_enrich TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#' @param annotation_obj TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_scores <- function(se,
                       res_de, # maybe won't be needed?
                       res_enrich,
                       genes_colname = "genes",
                       genesetname_colname = "Term",
                       genesetid_colname = "GO.ID",
                       annotation_obj = NULL) {

  mydata <- assay(se)
  # returns a matrix, rows = genesets, cols = samples
  rownames(res_enrich) <- res_enrich[[genesetid_colname]]

  rowsd_se <- matrixStats::rowSds(mydata)
  rowavg_se <- rowMeans(mydata)
  mydata_z <- (mydata - rowavg_se) / rowsd_se

  gss_mat <- matrix(NA, nrow = nrow(res_enrich), ncol = ncol(se))
  rownames(gss_mat) <- paste0(
    res_enrich[[genesetid_colname]], "|",
    res_enrich[[genesetname_colname]])
  colnames(gss_mat) <- colnames(se)

  for(i in seq_len(nrow(res_enrich))) {

    thisset_members <- unlist(strsplit(res_enrich[i,genes_colname],","))
    thisset_members_ids <- annotation_obj$gene_id[match(thisset_members, annotation_obj$gene_name)]

    thisset_members_ids <- thisset_members_ids[thisset_members_ids %in% rownames(se)]
    thisset_zs <- mydata_z[thisset_members_ids,]
    thisset_mean_zs <- colMeans(thisset_zs)

    gss_mat[i,] <- thisset_mean_zs
  }

  return(gss_mat)
}

#' Title
#'
#' TODO
#'
#' @param mat TODO
#' @param clustering_distance_rows TODO
#' @param clustering_distance_cols TODO
#' @param cluster_rows TODO
#' @param cluster_cols TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_ggheatmap <- function(mat,
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE
) {
  d_rows <- dist(mat, method = clustering_distance_rows)
  d_cols <- dist(t(mat), method = clustering_distance_cols)

  row_tree <- hclust(d_rows)
  col_tree <- hclust(d_cols)
  score_matrix <- mat

  if(cluster_rows)
    score_matrix <- score_matrix[row_tree$order, ]

  if(cluster_cols)
    score_matrix <- score_matrix[ , col_tree$order]

  labels_rows <- factor(rownames(score_matrix),
                        levels = rev(rownames(score_matrix)))
                        # to have the top ones on top, if not clustered
  labels_cols <- factor(colnames(score_matrix),
                        levels = colnames(score_matrix))

  score_df <- expand.grid(list(labels_rows, labels_cols),
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE
  )
  scores <- data.frame(as.vector(score_matrix))
  score_df <- cbind(score_df, scores)
  colnames(score_df) <- c("GeneSet", "Sample", "Value")

  p <- ggplot(score_df, aes_string(x = "Sample", y = "GeneSet")) +
    geom_tile(aes_string(fill = "Value"), color = "white") +
  # from brewer.pal(n = 11, name = 'RdYlBu') and brewer.pal(n = 11, name = 'YlGn') - or Spectal, 11?
    scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1
      ),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 11)
  )
  return(p)
}
