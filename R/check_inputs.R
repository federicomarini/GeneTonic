#' Checking the input objects for GeneTonic
#'
#' Checking the input objects for GeneTonic, whether these are all set for running
#' the app
#'
#' Some suggestions on the requirements for each parameter are returned in the
#' error messages.
#'
#' @param dds A `DESeqDataSet` object, normally obtained after running your data
#' through the `DESeq2` framework.
#' @param res_de A `DESeqResults` object. As for the `dds` parameter, this is
#' also commonly used in the `DESeq2` framework.
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param annotation_obj A `data.frame` object, containing two columns, `gene_id`
#' with a set of unambiguous identifiers (e.g. ENSEMBL ids) and `gene_name`,
#' containing e.g. HGNC-based gene symbols.
#'
#' @return Invisible NULL
#' @export
#'
#' @examples
#' #
checkup_GeneTonic <- function(dds,
                              res_de,
                              res_enrich,
                              annotation_obj) {

  # checking dds
  if (!is(dds, "DESeqDataSet"))
    stop("The provided `dds` is not a DESeqDataSet object, please check your input parameters")
  # TODO: check normalization factors/size factors are in?

  # checking res_de
  if (!is(res_de, "DESeqResults"))
    stop("The provided `res_de` is not a DESeqResults object, please check your input parameters")

  # checking res_enrich
  colnames_res_enrich <- c("gs_id",
                           "gs_description",
                           "gs_pvalue",
                           "gs_genes",
                           "gs_de_count",
                           "gs_bg_count") # gs_ontology? more?
  if (!all(colnames_res_enrich %in% colnames(res_enrich)))
    stop("The provided `res_enrich` object does not respect the format required to use in GeneTonic\n",
         "e.g. all required column names have to be present.\n",
         "You might want to use one of the `shaker_*` functions to convert it.\n",
         "Required columns: ", paste(colnames_res_enrich, collapse = ", "))

  # checking annotation_obj
  colnames_annotation_obj <- c("gene_id",
                               "gene_name")
  if (!all(colnames_annotation_obj %in% colnames(annotation_obj)))
    stop("The provided `annotation_obj` object does not respect the format required to use in GeneTonic\n",
         "e.g. all required column names have to be present.\n",
         "You can use e.g. `pcaExplorer::get_annotation_orgdb()` for this purpose.\n",
         "Required columns: ", paste(colnames_annotation_obj, collapse = ", "))

  message("All set to enjoy GeneTonic!")
  invisible(NULL)
}
