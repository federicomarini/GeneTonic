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
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included (from the top ranked ones). Defaults to the number of rows of
#' `res_enrich`
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be included, additionally to
#' the ones specified via `n_gs`. Defaults to NULL.
#' @param similarity_measure Character, currently defaults to `kappa_matrix`, to
#' specify how to compute the similarity measure between gene sets
#' @param mds_k Integer value, number of dimensions to compute in the multi
#' dimensional scaling procedure
#' @param mds_labels Integer, defines the number of labels to be plotted on top
#' of the scatter plot for the provided gene sets.
#' @param mds_colorby Character specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults sensibly to `z_score`.
#' @param gs_labels Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be labeled.
#' @param plot_title Character string, used as title for the plot. If left `NULL`,
#' it defaults to a general description of the plot and of the DE contrast
#' @param return_data Logical, whether the function should just return the
#' data.frame of the MDS coordinates, related to the original `res_enrich`
#' object. Defaults to FALSE.
#'
#' @return A `ggplot` object
#'
#' @seealso [create_kappa_matrix()] is used to calculate the similarity between
#' gene sets
#'
#' @export
#'
#' @examples
#'
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
#' gs_mds(res_enrich,
#'        res_de,
#'        anno_df,
#'        n_gs = 200,
#'        mds_labels = 10)
gs_mds <- function(res_enrich,
                   res_de,
                   annotation_obj,
                   n_gs = nrow(res_enrich),
                   gs_ids = NULL,
                   similarity_measure = "kappa_matrix",
                   mds_k = 2,
                   mds_labels = 0,
                   mds_colorby = "z_score",
                   gs_labels = NULL,
                   plot_title = NULL,
                   return_data = FALSE) { # or aggr_score

  # TODO: match.arg on similarity matrix?
  similarity_measure <- match.arg(similarity_measure,
                                  c("kappa_matrix", "overlap_matrix"))

  # require res_enrich to have aggregated scores and so
  if (!("z_score" %in% colnames(res_enrich))) {
    res_enrich <- get_aggrscores(res_enrich,
                                 res_de,
                                 n_gs = nrow(res_enrich),
                                 annotation_obj = annotation_obj)
  }

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )


  # alternative: use semantic similarity
  # library(GOSemSim)
  # mmGO <- godata('org.Mm.eg.db', ont="BP")
  # mysims <- mgoSim(res_enrich[["gs_id"]], res_enrich[["gs_id"]],
  # semData=mmGO, measure="Wang", combine=NULL)

  if (similarity_measure == "kappa_matrix") {
    my_simmat <- create_kappa_matrix(res_enrich, n_gs = n_gs, gs_ids = gs_ids)

  } else if (similarity_measure == "overlap_matrix") {
    my_simmat <- create_jaccard_matrix(res_enrich, n_gs = n_gs, gs_ids = gs_ids, return_sym = TRUE)
  }

  # else ... TODO


  # subset here, internally
  res_enrich <- res_enrich[gs_to_use, ]
  mds_labels <- min(mds_labels, nrow(res_enrich))

  mysets <- res_enrich[["gs_id"]]
  mysets_names <- res_enrich[["gs_description"]]


  mds_gs <- cmdscale((1 - my_simmat), eig = TRUE, k = mds_k)

  # TODOTODO: cbind it to the original data frame?
  mds_gs_df <- data.frame(
    dim1 = mds_gs$points[, 1],
    dim2 = mds_gs$points[, 2], # TODO: handle 3rd dim?
    gs_id = mysets,
    gs_name = mysets_names,
    gs_DEcount = res_enrich$DE_count,
    gs_colby = res_enrich[[mds_colorby]],
    gs_text = paste0(mysets, ": ", mysets_names),  # TODOTODO: is there a way to avoid the warning from gg?
    stringsAsFactors = FALSE
  )

  if (return_data) {
    # TODO: maybe
    return(mds_gs_df)
  }

  max_z <- max(abs(range(mds_gs_df$gs_colby)))
  limit <- max_z * c(-1, 1)

  this_contrast <- (sub(".*p-value: (.*)", "\\1", mcols(res_de, use.names = TRUE)["pvalue", "description"]))

  p <- ggplot(mds_gs_df, aes_string(x = "dim1",
                                    y = "dim2",
                                    text = "gs_text")) +
    geom_point(aes_string(color = "gs_colby",
                          size = "gs_DEcount")) +
    scale_color_gradient2(limit = limit,
                          low = muted("deepskyblue"),
                          mid = "lightyellow",
                          high = muted("firebrick")) +
    labs(x = "Dimension 1",
         y = "Dimension 2",
         col = mds_colorby,
         size = "Geneset \nmembers"
    ) +
    theme_bw()

  if (mds_labels > 0) {
    label_these <- mds_gs_df$gs_id[1:mds_labels]
  } else {
    label_these <- NULL
  }

  if (!is.null(gs_labels)) {
    if (!all(gs_labels %in% res_enrich$gs_id)) {
      warning("Not all specified geneset ids were found in the `res_enrich` object")
    }
    label_those <- mds_gs_df[mds_gs_df$gs_id %in% gs_labels, "gs_id"]
  } else {
    label_those <- NULL
  }

  df_gs_labels <- mds_gs_df[mds_gs_df$gs_id %in% unique(c(label_these, label_those)), ]

  p <- p + geom_label_repel(
    aes_string(label = "gs_name"), data = df_gs_labels, size = 3)

  # handling the title
  if (is.null(plot_title)) {
    p <- p + ggtitle(paste0("Geneset MDS plot - ", this_contrast))
  } else {
    p <- p + ggtitle(plot_title)
  }

  return(p)
  ## also something to obtain clusters of terms? - well, the colors do it somehow already
}
