#' Title
#'
#' @param res_enrich
#' @param res_de
#' @param n_gs
#' @param annotation_obj
#' @param genes_colname
#' @param genesetname_colname
#' @param genesetid_colname
#'
#' @return
#' @export
#'
#' @examples
gs_summary_heat <- function (res_enrich,
                             res_de,
                             n_gs = 80,
                             annotation_obj,
                             genes_colname = "genes",
                             genesetname_colname = "Term",
                             genesetid_colname = "GO.ID") {

  res_enrich2 <- res_enrich[seq_len(n_gs),]

  enriched_gsids <- res_enrich2[[genesetid_colname]]
  enriched_gsnames <- res_enrich2[[genesetname_colname]]
  enriched_gsdescs <- vapply(enriched_gsids, function(arg) Definition(GOTERM[[arg]]), character(1))

  # rownames(res_enrich) <- enriched_gsids

  gs_expanded <- tidyr::separate_rows(res_enrich2, {{genes_colname}}, sep = ",")
  gs_expanded$log2FoldChange <-
    res_de[annotation_obj$gene_id[match(gs_expanded$genes, annotation_obj$gene_name)], ]$log2FoldChange

  # keep them as factor to prevent rearrangement!
  gs_expanded[[genesetid_colname]] <- factor(gs_expanded[[genesetid_colname]], levels = res_enrich2[[genesetid_colname]])
  gs_expanded[[genesetname_colname]] <- factor(gs_expanded[[genesetname_colname]], levels = res_enrich2[[genesetname_colname]])
  gs_expanded[[genes_colname]] <- factor(gs_expanded[[genes_colname]], levels = unique(gs_expanded[[genes_colname]]))

  p <- ggplot(gs_expanded,
              aes_string(genes_colname, genesetname_colname)) +
    geom_tile(aes_string(fill = "log2FoldChange"),
              color = "white") +
    scale_fill_gradient2(low = muted("deepskyblue"),
                         mid = "lightyellow",
                         high = muted("firebrick"),
                         name = "log2FoldChange") +
    xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 75, hjust = 1))

  return(p)
}

# gs_summary_heat(res_enrich,res_de,annotation_obj = anno_df, n_gs = 30)
# gs_summary_heat(res_enrich,res_de,annotation_obj = anno_df, n_gs = 30) %>% plotly::ggplotly()
