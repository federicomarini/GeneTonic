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
#' @importFrom bs4Dash bs4Card bs4DashBody bs4DashControlbar
#' bs4DashFooter bs4DashNavbar bs4DashPage bs4DashSidebar
#' bs4InfoBox bs4InfoBoxOutput bs4SidebarMenu
#' bs4SidebarMenuItem bs4TabItem bs4TabItems bs4ValueBox
#' bs4ValueBoxOutput renderbs4InfoBox renderbs4ValueBox
#' bs4TabPanel bs4TabCard
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
