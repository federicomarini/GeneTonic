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
#' @family pathway-analysis-results
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name topgoDE_macrophage_IFNg_vs_naive
#' @docType data
NULL


#' A sample output from Enrichr
#'
#' A sample output object as created from a call to Enrichr, with the interface
#' provided by `enrichR` - using the `enrichr()` function
#'
#' @details This object has been created on the data from the `macrophage` package
#' by analyzing downstream the differentially expressed genes
#' when comparing IFNg treated samples vs naive samples, accounting
#' for the different cell lines included.
#'
#' Details on how this object has been created are included in the `create_gt_data.R`
#' script, included in the `scripts` folder of the `GeneTonic` package.
#'
#' @family pathway-analysis-results
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name enrichr_output_macrophage
#' @docType data
NULL



#' A sample output from g:Profiler
#'
#' A sample output object as created from a call to g:Profiler, with the interface
#' provided by `gprofiler2` - using the `gost()` function
#'
#' @details This object has been created on the data from the `macrophage` package
#' by analyzing downstream the differentially expressed genes
#' when comparing IFNg treated samples vs naive samples, accounting
#' for the different cell lines included.
#'
#' Details on how this object has been created are included in the `create_gt_data.R`
#' script, included in the `scripts` folder of the `GeneTonic` package.
#'
#' @family pathway-analysis-results
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name gostres_macrophage
#' @docType data
NULL


#' A sample output from fgsea
#'
#' A sample output object as created from a call to the `fgsea()` function, in
#' the `fgsea` package, as a practical framework for performing GSEA
#'
#' @details This object has been created on the data from the `macrophage` package
#' by analyzing downstream the differentially expressed genes
#' when comparing IFNg treated samples vs naive samples, accounting
#' for the different cell lines included.
#' 
#' @family pathway-analysis-results
#'
#' Details on how this object has been created are included in the `create_gt_data.R`
#' script, included in the `scripts` folder of the `GeneTonic` package.
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name fgseaRes
#' @docType data
NULL
