#' GeneTonic
#'
#' `GeneTonic` is a Bioconductor package that provides an interactive Shiny-based
#'  graphical user interface for...
#'
#' @author Federico Marini \email{marinif@@uni-mainz.de}
#'
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @importFrom DT dataTableOutput renderDataTable datatable
#' @importFrom rintrojs introjs introjsUI
#' @import igraph
#' @import visNetwork
#' @importFrom utils read.delim sessionInfo
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales alpha muted
#' @importFrom S4Vectors mcols
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @import GO.db
#' @import DESeq2
#' @importFrom AnnotationDbi Definition GOID Term Synonym Secondary
#' @importFrom plotly ggplotly plotlyOutput renderPlotly plot_ly layout add_trace
#' @importFrom tidyr separate_rows pivot_longer
#' @importFrom stats var dist hclust
#' @import pheatmap
#' @importFrom matrixStats rowSds
#' @import SummarizedExperiment
#' @importFrom methods is
#' @import pcaExplorer
#' @import shinycssloaders
#' @import bs4Dash
#' @import shinyWidgets
#' @import ggalluvial
#' @import viridis
#' @importFrom dplyr arrange desc group_by mutate
#' @importFrom ggforce geom_sina
#' @importFrom rlang .data
#'
#' @name GeneTonic-pkg
#' @docType package
NULL
