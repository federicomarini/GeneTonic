#' A sample `DESeqResults` object
#'
#' A sample `DESeqResults` object, generated in the `DESeq2` framework
#'
#' @details This `DESeqResults` object on the data from the `macrophage` package
#' has been created comparing IFNg treated samples vs naive samples, accounting
#' for the different cell lines included.
#'
#' Details on how this object has been created are included in the `create_gt_data.R`
#' script, included in the `scripts` folder of the `GeneTonic` package.
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name res_macrophage_IFNg_vs_naive
#' @docType data
NULL


#' A sample `res_enrich` object
#'
#' A sample `res_enrich` object, generated with the `topGOtable` function (from
#' the `pcaExplorer` package).
#'
#' @details This `res_enrich` object on the data from the `macrophage` package
#' has been created by analyzing downstream the differentially expressed genes
#' when comparing IFNg treated samples vs naive samples, accounting
#' for the different cell lines included.
#'
#' Details on how this object has been created are included in the `create_gt_data.R`
#' script, included in the `scripts` folder of the `GeneTonic` package.
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name topgoDE_macrophage_IFNg_vs_naive
#' @docType data
NULL
