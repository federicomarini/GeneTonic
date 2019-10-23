#' Title TODO
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param res_de TODO
#' @param n_gs TODO
#' @param annotation_obj TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
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



# or something like
# TODO TODO https://www.biostars.org/p/299161/


# gs_summary_heat(res_enrich,res_de,annotation_obj = anno_df, n_gs = 30)
# gs_summary_heat(res_enrich,res_de,annotation_obj = anno_df, n_gs = 30) %>% plotly::ggplotly()
#
#
# gs_summary_dots <- function (object,
#                              x = "geneRatio",
#                              color = "p.adjust",
#                              showCategory = 10,
#                              size = NULL,
#                              split = NULL,
#                              font.size = 12,
#                              title = "",
#                              orderBy = "x",
#                              decreasing = TRUE)
# {
#   #
#   # enrich2list <- lapply(seq_len(n_gs), function(gs) { # goterm <-
#   # res_enrich$Term[gs] go_genes <- res_enrich2$genes[gs] go_genes <-
#   # strsplit(go_genes, ",") %>% unlist return(go_genes) }) names(enrich2list) <-
#   # enriched_gsids[seq_len(n_gs)]
#
#
#
#
#   colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue","Zscore","aggr_score"))
#   if (x == "geneRatio" || x == "GeneRatio") {
#     x <- "GeneRatio"
#     if (is.null(size))
#       size <- "Count"
#   }
#   else if (x == "count" || x == "Count") {
#     x <- "Count"
#     if (is.null(size))
#       size <- "GeneRatio"
#   }
#   else if (is(x, "formula")) {
#     x <- as.character(x)[2]
#     if (is.null(size))
#       size <- "Count"
#   }
#   else {
#     if (is.null(size))
#       size <- "Count"
#   }
#   df <- fortify(object, showCategory = showCategory, split = split)
#   if (!orderBy %in% colnames(df)) {
#     message("wrong orderBy parameter; set to default `orderBy = \"x\"`")
#     orderBy <- "x"
#   }
#   if (orderBy == "x") {
#     df <- dplyr::mutate(df, x = eval(parse(text = x)))
#   }
#   idx <- order(df[[orderBy]], decreasing = decreasing)
#   df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))
#   ggplot(df, aes_string(x = x, y = "Description", size = size,
#                         color = colorBy)) + geom_point() + scale_color_continuous(low = "red",
#                                                                                   high = "blue", name = color, guide = guide_colorbar(reverse = TRUE)) +
#     ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
#     scale_size(range = c(3, 8))
# }
#
