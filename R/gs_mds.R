#' Multi Dimensional Scaling plot for gene sets
#'
#' Multi Dimensional Scaling plot for gene sets, extracted from a `res_enrich`
#' object
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, `GeneTonic`, to see the
#' formatting requirements.
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param genes_colname Character, specifying which column of the `res_enrich`
#' object contains the genes assigned to each gene set, detected as differentially
#' expressed. Defaults to `genes`.
#' @param genesetname_colname Character, specifies which column of the `res_enrich`
#' object contains a description of the gene set. Defaults to `Term`.
#' @param genesetid_colname Character, specifies which column of the `res_enrich`
#' object contains a unique identifier of the gene set. Defaults to `GO.ID`.
#' @param genes_separator Character, specifying which separator is used in the
#' column defined by `genes_colname` to split the character of features.
#' @param similarity_measure Character, currently defaults to `kappa_matrix`, to
#' specify how to compute the similarity measure between gene sets
#' @param mds_k Integer value, number of dimensions to compute in the multi
#' dimensional scaling procedure
#' @param mds_labels Integer, defines the number of labels to be plotted on top
#' of the scatter plot for the provided gene sets.
#' @param mds_colorby Character specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults sensibly to `z_score`.
#'
#' @return A `ggplot` object
#' @export
#'
#' @examples
#' # TODO
gs_mds <- function(res_enrich,
                   res_de,
                   annotation_obj,
                   genes_colname = "genes",
                   genesetname_colname = "Term",
                   genesetid_colname = "GO.ID",
                   genes_separator = ",",
                   similarity_measure = "kappa_matrix",
                   mds_k = 2,
                   mds_labels = 10,
                   mds_colorby = "z_score") { # or aggr_score

  # require res_enrich to have aggregated scores and so
  if (!("Z_score" %in% colnames(res_enrich))) {
    res_enrich <- get_aggrscores(res_enrich,
                                 res_de,
                                 n_gs = nrow(res_enrich),
                                 genes_colname = genes_colname,
                                 genesetname_colname = genesetname_colname,
                                 genesetid_colname = genesetid_colname,
                                 annotation_obj = annotation_obj)
  }

  # alternative: use semantic similarity
  # library(GOSemSim)
  # mmGO <- godata('org.Mm.eg.db', ont="BP")
  # mysims <- mgoSim(res_enrich[[genesetid_colname]], res_enrich[[genesetid_colname]],
  # semData=mmGO, measure="Wang", combine=NULL)

  mysets <- res_enrich[[genesetid_colname]]
  mysets_names <- res_enrich[[genesetname_colname]]

  if (similarity_measure == "kappa_matrix") {
    my_simmat <- create_kappa_matrix(res_enrich,
                                     genes_colname = genes_colname,
                                     genesetname_colname = genesetname_colname,
                                     genesetid_colname = genesetid_colname,
                                     genes_separator = genes_separator)

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

  if (!is.null(mds_labels)) {
    p <- p + geom_label_repel(
      aes_string(label = "gs_name"), data = mds_go_df[1:mds_labels, ], size = 3)
  }

  return(p)
  ## also something to obtain clusters of terms? - well, the colors do it somehow already
}

# then ggplotly(p)
