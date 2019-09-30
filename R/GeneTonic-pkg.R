#' GeneTonic
#'
#' `GeneTonic` is a Bioconductor package that provides an interactive Shiny-based
#'  graphical user interface for...
#'
#' @author Federico Marini \email{marinif@@uni-mainz.de}
#'
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @import shinydashboard
#' @importFrom DT dataTableOutput renderDataTable datatable
#' @importFrom rintrojs introjs introjsUI
#' @import igraph
#' @import visNetwork
#' @importFrom utils read.delim
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales alpha muted
#' @importFrom S4Vectors mcols
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @import GO.db
#' @importFrom AnnotationDbi Definition
#' @importFrom plotly ggplotly plotlyOutput renderPlotly
#' @importFrom tidyr separate_rows
#'
#' @name GeneTonic-pkg
#' @docType package
NULL
