#' Title TODO
#'
#' TODO
#'
#' @param res_enrich  TODO
#' @param res_de  TODO
#' @param annotation_obj TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#' @param genes_separator TODO
#' @param similarity_measure TODO
#' @param mds_k TODO
#' @param mds_labels TODO
#' @param mds_colorby TODO
#'
#' @return TODO
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
    p <- p + geom_label_repel(aes_string(label = "gs_name"), data = mds_go_df[1:mds_labels, ], size = 3)
  }

  return(p)
  ## also something to obtain clusters of terms? - well, the colors do it somehow already
}

# then ggplotly(p)
