
#' Title
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param res_de TODO
#' @param n_gs TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#' @param annotation_obj TODO
#' @param chars_limit TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' #TODO
enhance_table <- function(res_enrich,
                          res_de,
                          n_gs = 50,
                          genes_colname = "genes",
                          genesetname_colname = "Term",
                          genesetid_colname = "GO.ID",
                          annotation_obj,
                          chars_limit = 60) {

  # res_enrich has to have a column called containing the genes annotated to the term
  # TODOTODO

  # verify the genesets are sorted in a meaningful way?
  #TODOTODO

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_fulllist <- lapply(seq_len(n_gs), function(go) {
    genes_thisset <- res_enrich[[genes_colname]][go]
    genes_thisset <- unlist(strsplit(genes_thisset, ","))

    genesid_thisset <- annotation_obj$gene_id[match(genes_thisset, annotation_obj$gene_name)]

    res_thissubset <- res_de[genesid_thisset, ]
    res_thissubset$gene_name <- genes_thisset
    res_thissubset$goterm <- as.factor(res_enrich[[genesetname_colname]][go])
    res_thissubset$gotermshort <- substr(res_enrich[[genesetname_colname]][go], 1, 50)
    res_thissubset$goid <- res_enrich[[genesetid_colname]][go]
    return(as.data.frame(res_thissubset))
  })
  gs_fulllist <- do.call(rbind, gs_fulllist)

  this_contrast <- (sub(".*p-value: (.*)", "\\1", mcols(res_de, use.names = TRUE)["pvalue", "description"]))

  # to have first rows viewed on top
  gs_fulllist <- gs_fulllist[nrow(gs_fulllist):1, ]
  gs_fulllist$goterm <- factor(gs_fulllist$goterm, levels = rev(levels(gs_fulllist$goterm)))
  max_lfc <- max(abs(range(gs_fulllist$log2FoldChange)))

  # z score, and evtl. option to sort by that? TODOTODO

  p <- ggplot(
    gs_fulllist, aes_string(
      x = "log2FoldChange",
      y = "goterm",
      fill = "goid",
      text = "gene_name"
    )) +
    ggtitle(paste0(this_contrast, " - TODOTODO")) +
    scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
    geom_point(alpha = 0.7, shape = 21, size = 3)  +
    theme_minimal() +
    geom_vline(aes(xintercept = 0), col = "steelblue", alpha = 0.4) +
    theme(legend.position = "none") +
    scale_y_discrete(name = "",
                     labels = substr(as.character(unique(gs_fulllist$goterm)), 1, chars_limit))

  return(p)
}


#' Title TODO
#' TODO
#' @param res_enrich TODO
#' @param res_de TODO
#' @param n_gs TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#' @param annotation_obj TODO
#' @param aggrfun TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
#'
#' # TODOTODO: add some standardized column names, say "count"
get_aggrscores <- function(res_enrich,
                           res_de,
                           n_gs = 50,
                           genes_colname = "genes",
                           genesetname_colname = "Term",
                           genesetid_colname = "GO.ID",
                           annotation_obj,
                           aggrfun = mean) {

  # allgenes <- unlist(strsplit(res_enrich[[genes_colname]],","))
  gs_expanded <- tidyr::separate_rows(res_enrich, {{genes_colname}}, sep = ",")
  gs_expanded$log2FoldChange <-
    res_de[annotation_obj$gene_id[match(gs_expanded$genes, annotation_obj$gene_name)], ]$log2FoldChange

  gs_aggregated <- lapply(seq_len(nrow(res_enrich)), function(i) {
    this_gsid <- res_enrich[[genesetid_colname]][i]
    this_genesetname <- res_enrich[[genesetname_colname]][i]
    this_subset <- gs_expanded[gs_expanded[[genesetid_colname]] == this_gsid, ]

    upgenes <- sum(this_subset$log2FoldChange > 0)
    downgenes <- sum(this_subset$log2FoldChange < 0)
    z_score <- (upgenes - downgenes) / sqrt(upgenes + downgenes)

    aggr_score <- aggrfun(this_subset$log2FoldChange)
    return(c("Z_score" = z_score,
             "aggr_score" = aggr_score))
  })

  names(gs_aggregated) <- res_enrich[[genesetid_colname]]

  res_enrich$z_score <- vapply(gs_aggregated, "[", 1, FUN.VALUE = numeric(1))
  res_enrich$aggr_score <- vapply(gs_aggregated, "[", 2, FUN.VALUE = numeric(1))

  return(res_enrich)
}
