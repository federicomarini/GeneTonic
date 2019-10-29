

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
  p <- ggplot(re_sorted, ( aes_string(x="Term", y="logp10"))) +
    geom_segment( aes_string(x="Term" ,xend="Term", y=0, yend="logp10"), color="grey") +
    geom_point(aes_string(col="z_score"), size = 4 ) +
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
  re1$logp10 <- -log10(res_enrich[[p_value_column]])

  set.seed(42)
  shuffled_ones <- sample(seq_len(nrow(re1)))
  re2 <- res_enrich
  re2$z_score <- re2$z_score[shuffled_ones]
  re2$aggr_score <- re2$aggr_score[shuffled_ones]
  re2$logp10 <- re1$logp10[shuffled_ones]

  re_both <- mutate(re1,
                    z_score_2 = re2$z_score,
                    aggr_score_2 = re2$aggr_score,
                    logp10_2 = re2$logp10)
  re_both <- re_both[seq_len(n_gs), ]

  re_both_sorted <- re_both %>%
    arrange(logp10) %>%
    mutate(Term=factor(Term, Term))

  p <- ggplot(re_both_sorted, aes_string(x="Term", y="logp10")) +
    geom_segment( aes_string(x="Term" ,xend="Term", y="logp10_2", yend="logp10"), color="grey") +
    geom_point(aes_string(col="z_score"), size=4 ) +
    geom_point(aes_string(y = "logp10_2", col = "z_score_2"), size = 4, alpha = alpha_set2) +
    scale_color_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    coord_flip() +
    theme_minimal()

  return(p)

}

# TODO: size of points somewhat related to geneset size?

#' Title TODO
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param n_gs TODO
#' @param p_value_column TODO
#' @param color_by TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_horizon <- function(res_enrich, # TODO: should be a list of res_enrich objects!
                       n_gs = 20,
                       p_value_column = "p.value_elim",
                       color_by = "z_score") {
  # again, must be enhanced with Zscore


  # res_enrich <- get_aggrscores(topgoDE_macrophage_IFNg_vs_naive,res_macrophage_IFNg_vs_naive, annotation_obj = anno_df)
  res_enriched_1 <- res_enrich

  res_enriched_1 <- res_enriched_1[seq_len(n_gs),]
  res_enriched_1$logp10 <-  -log10(res_enriched_1[[p_value_column]])

  res_enriched_2 <-
    res_enriched_3 <-
    res_enriched_4 <- res_enriched_1

  set.seed(42)
  shuffled_r2 <- sample(seq_len(nrow(res_enriched_1)))
  shuffled_r3 <- sample(seq_len(nrow(res_enriched_1)))
  shuffled_r4 <- sample(seq_len(nrow(res_enriched_1)))

  res_enriched_2$z_score <- res_enriched_2$z_score[shuffled_r2]
  res_enriched_2$aggr_score <- res_enriched_2$aggr_score[shuffled_r2]
  res_enriched_2$logp10 <- res_enriched_2$logp10[shuffled_r2]
  res_enriched_3$z_score <- res_enriched_3$z_score[shuffled_r3]
  res_enriched_3$aggr_score <- res_enriched_3$aggr_score[shuffled_r3]
  res_enriched_3$logp10 <- res_enriched_3$logp10[shuffled_r3]
  res_enriched_4$z_score <- res_enriched_4$z_score[shuffled_r4]
  res_enriched_4$aggr_score <- res_enriched_4$aggr_score[shuffled_r4]
  res_enriched_4$logp10 <- res_enriched_4$logp10[shuffled_r4]

  res_enriched_1$scenario <- "original"
  res_enriched_2$scenario <- "shuffled_2"
  res_enriched_3$scenario <- "shuffled_3"
  res_enriched_4$scenario <- "shuffled_4"

  # to preserve the order of the terms
  res_enriched_1 <- res_enriched_1 %>%
    arrange(logp10) %>%
    mutate(Term=factor(Term, unique(Term)))

  # to preserve the sorting of scenarios
  merged_res_enh <- rbind(res_enriched_1,
                          res_enriched_2,
                          res_enriched_3,
                          res_enriched_4)
  merged_res_enh$scenario <- factor(merged_res_enh$scenario, unique(merged_res_enh$scenario))


  # if only with one...
  res_enriched_1 %>%
    # arrange(logp10) %>%
    # mutate(Term=factor(Term, unique(Term))) %>%
    ggplot(aes_string(x = "Term", y = "logp10") ) +
    geom_line(aes_string(group = "scenario", col = "scenario"), size = 3, alpha = 0.7) +
    geom_point(aes_string(fill="z_score"), size=4, pch = 21 ) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    ylim(c(0,NA)) +
    coord_flip() +
    theme_minimal()

  # sorted by category in scenario1
  merged_res_enh %>%
    mutate(Term=factor(.data$Term, unique(.data$Term))) %>%
    arrange(desc(.data$logp10)) %>%
    ggplot(aes_string(x = "Term", y = "logp10") ) +
    geom_line(aes_string(group = "scenario", col = "scenario"), size = 3, alpha = 0.7) +
    geom_point(aes_string(fill = "z_score"), size=4, pch = 21 ) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    ylim(c(0,NA)) +
    coord_flip() +
    theme_minimal()

  # with a nicer sorting - "grouped" by scenario

  nicerorder_terms <- merged_res_enh %>%
    group_by(.data$Term) %>%
    mutate(main_category = .data$scenario[which.max(.data$logp10)],
           max_value = max(.data$logp10)) %>%
    arrange(.data$main_category, desc(.data$max_value)) %>%
    dplyr::pull(.data$Term)



  merged_res_enh %>%
    # mutate(Term=factor(Term, unique(Term))) %>%
    mutate(Term=factor(.data$Term, rev(unique(nicerorder_terms)))) %>%
    arrange(desc(.data$logp10)) %>%
    ggplot(aes_string(x = "Term", y = "logp10") ) +
    geom_line(aes_string(group = "scenario", col = "scenario"), size = 3, alpha = 0.7) +
    geom_point(aes_string(fill="z_score"), size=4, pch = 21 ) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    ylim(c(0,NA)) +
    coord_flip() +
    theme_minimal()


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
