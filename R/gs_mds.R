#' Multi Dimensional Scaling plot for gene sets
#'
#' Multi Dimensional Scaling plot for gene sets, extracted from a `res_enrich`
#' object
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param similarity_measure Character, currently defaults to `kappa_matrix`, to
#' specify how to compute the similarity measure between gene sets
#' @param mds_k Integer value, number of dimensions to compute in the multi
#' dimensional scaling procedure
#' @param mds_labels Integer, defines the number of labels to be plotted on top
#' of the scatter plot for the provided gene sets.
#' @param mds_colorby Character specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults sensibly to `z_score`.
#' @param gs_labels Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be labelled.
#'
#' @return A `ggplot` object
#'
#' @seealso [create_kappa_matrix()] is used to calculate the similarity between
#' gene sets
#'
#' @export
#'
#' @examples
#' # TODO
gs_mds <- function(res_enrich,
                   res_de,
                   annotation_obj,
                   similarity_measure = "kappa_matrix",
                   mds_k = 2,
                   mds_labels = 0,
                   mds_colorby = "z_score",
                   gs_labels = NULL) { # or aggr_score

  # require res_enrich to have aggregated scores and so
  if (!("z_score" %in% colnames(res_enrich))) {
    res_enrich <- get_aggrscores(res_enrich,
                                 res_de,
                                 n_gs = nrow(res_enrich),
                                 annotation_obj = annotation_obj)
  }

  # alternative: use semantic similarity
  # library(GOSemSim)
  # mmGO <- godata('org.Mm.eg.db', ont="BP")
  # mysims <- mgoSim(res_enrich[["gs_id"]], res_enrich[["gs_id"]],
  # semData=mmGO, measure="Wang", combine=NULL)
  mds_labels <- min(mds_labels, nrow(res_enrich))


  mysets <- res_enrich[["gs_id"]]
  mysets_names <- res_enrich[["gs_description"]]

  if (similarity_measure == "kappa_matrix") {
    my_simmat <- create_kappa_matrix(res_enrich)

  } # else ... TODO

  mds_go <- cmdscale((1 - my_simmat), eig = TRUE, k = mds_k)

  # TODOTODO: cbind it to the original data frame?
  mds_go_df <- data.frame(
    dim1 = mds_go$points[, 1],
    dim2 = mds_go$points[, 2], # handle 3rd dim?
    gs_id = mysets,
    gs_name = mysets_names,
    gs_DEcount = res_enrich$DE_count,
    gs_colby = res_enrich[[mds_colorby]],
    text = paste0(mysets, ": ", mysets_names) # TODOTODO: is there a way to avoid the warning from gg?
  )

  max_z <- max(abs(range(mds_go_df$gs_colby)))
  limit <- max_z * c(-1, 1)

  p <- ggplot(mds_go_df, aes_string(x = "dim1",
                                    y = "dim2",
                                    text = "text")) +
    geom_point(aes_string(color = "gs_colby",
                          size = "gs_DEcount")) +
    scale_color_gradient2(limit = limit,
                          low = muted("deepskyblue"),
                          mid = "lightyellow",
                          high = muted("firebrick")) +
    theme_bw()

  if (mds_labels > 0) {
    p <- p + geom_label_repel(
      aes_string(label = "gs_name"), data = mds_go_df[1:mds_labels, ], size = 3)
  }

  if (!is.null(gs_labels)) {
    if (!all(gs_labels %in% res_enrich$gs_id)) {
      warning("Not all specified geneset ids were found in the `res_enrich` object")
    }
    df_gs_labels <- mds_go_df[mds_go_df$gs_id %in% gs_labels,]

    p <- p + geom_label_repel(
      aes_string(label = "gs_name"), data = df_gs_labels, size = 3)
  }



  return(p)
  ## also something to obtain clusters of terms? - well, the colors do it somehow already
}

# then ggplotly(p)
