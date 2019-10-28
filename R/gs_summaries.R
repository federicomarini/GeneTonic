

#' Title TODO
#'
#' TODO
#'
#' @param res_enrich  TODO
#' @param n_gs  TODO
#' @param p_value_column  TODO
#' @param color_by  TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_summary_overview <- function(res_enrich,
                                n_gs = 20,
                                p_value_column = "p.value_elim",
                                color_by = "z_score") {
  if (!("z_score" %in% colnames(res_enrich))) {
    warning("You need to add the z_score or the aggregated score")
  } # TODO: same for aggr_score

  re <- res_enrich
  re$logp10 <- -log10(res_enrich[[p_value_column]])
  re <- re[seq_len(n_gs), ]

  # TODO: color_by to be selected

  re_sorted <- re %>%
    arrange(logp10) %>%
    mutate(Term=factor(Term, Term))
  p <- ggplot(re_sorted, ( aes(x=Term, y=logp10))) +
    geom_segment( aes(x=Term ,xend=Term, y=0, yend=logp10), color="grey") +
    geom_point(aes(col=z_score), size = 4 ) +
    scale_color_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    coord_flip() +
    theme_minimal()

  return(p)
}

#' Title  TODO
#'
#' TODO
#'
#' @param res_enrich  TODO
#' @param res_enrich2  TODO
#' @param n_gs  TODO
#' @param p_value_column  TODO
#' @param color_by  TODO
#' @param alpha_set2  TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_summary_overview_pair <- function(res_enrich,
                                     res_enrich2,
                                n_gs = 20,
                                p_value_column = "p.value_elim",
                                color_by = "z_score",
                                alpha_set2 = 0.4
                                ) {
  if (!("z_score" %in% colnames(res_enrich))) {
    warning("You need to add the z_score or the aggregated score")
  } # TODO: same for aggr_score

  # TODO: require that both res_enrich have the same terms in the table
  # identical(re1$GO.ID,re2$GO.ID) # or so

  re1 <- res_enrich

  set.seed(42)
  shuffled_ones <- sample(seq_len(nrow(re1)))
  re2 <- res_enrich
  re2$z_score <- re2$z_score[shuffled_ones]
  re2$aggr_score <- re2$aggr_score[shuffled_ones]
  re2$logp10 <- re2$logp10[shuffled_ones]

  re_both <- mutate(re1,
                    z_score_2 = re2$z_score,
                    aggr_score_2 = re2$aggr_score,
                    logp10_2 = re2$logp10)
  re_both <- re_both[seq_len(n_gs), ]

  re_both_sorted <- re_both %>%
    arrange(logp10) %>%
    mutate(Term=factor(Term, Term))

  p <- ggplot(re_both_sorted, aes(x=Term, y=logp10)) +
    geom_segment( aes(x=Term ,xend=Term, y=logp10_2, yend=logp10), color="grey") +
    geom_point(aes(col=z_score), size=4 ) +
    geom_point(aes(y = logp10_2, col = z_score_2), size = 4, alpha = alpha_set2) +
    scale_color_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    coord_flip() +
    theme_minimal()

  return(p)

}





#' Title TODO
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param res_de TODO
#' @param annotation_obj TODO
#' @param n_gs TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_summary_heat <- function(res_enrich,
                            res_de,
                            annotation_obj,
                            n_gs = 80,
                            genes_colname = "genes",
                            genesetname_colname = "Term",
                            genesetid_colname = "GO.ID") {

  res_enrich2 <- res_enrich[seq_len(n_gs), ]

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
