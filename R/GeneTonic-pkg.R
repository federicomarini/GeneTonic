#' GeneTonic
#'
#' `GeneTonic` is a Bioconductor package that provides an interactive Shiny-based
#' graphical user interface for streamlining the interpretation of RNA-seq data
#'
#' `GeneTonic` simplifies and optimizes the integration of all components of
#' Differential Expression analysis, with functional enrichment analysis and the
#' original expression quantifications.
#' It does so in a way that makes it easier to generate insightful observations
#' and hypothesis - combining the benefits of interactivity and reproducibility,
#' e.g. by capturing the features and gene sets of interest highlighted during
#' the live session, and creating an HTML report as an artifact where text,
#' code, and output coexist.
#'
#' @importFrom AnnotationDbi Definition GOID Ontology Secondary Synonym Term
#' @importFrom backbone backbone_from_projection
#' @importFrom bs4Dash bs4Card bs4DashBody bs4DashControlbar
#' bs4DashFooter bs4DashNavbar bs4DashPage bs4DashSidebar
#' bs4InfoBox bs4InfoBoxOutput bs4SidebarMenu
#' bs4SidebarMenuItem bs4TabItem bs4TabItems bs4ValueBox
#' bs4ValueBoxOutput renderbs4InfoBox renderbs4ValueBox
#' bs4TabCard
#' @importFrom circlize colorRamp2
#' @importFrom colorspace rainbow_hcl
#' @importFrom colourpicker colourInput
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom ComplexUpset upset intersection_matrix intersection_size
#' @importFrom dendextend branches_attr_by_clusters set
#' @importFrom DESeq2 vst counts estimateSizeFactors normalizationFactors sizeFactors
#' @importFrom dplyr arrange desc group_by mutate pull slice select "%>%"
#' @importFrom DT datatable dataTableOutput renderDataTable formatRound
#' formatStyle JS
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom expm "%^%"
#' @importFrom ggforce geom_sina
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom ggridges geom_density_ridges position_points_jitter
#' @import GO.db
#' @importFrom graphics par plot
#' @importFrom grDevices colorRampPalette rgb col2rgb
#' @importFrom grid gpar
#' @importFrom igraph add_edges as_adjacency_matrix as_biadjacency_matrix degree
#' delete_edges delete_vertices "%du%" E "E<-" graph_from_data_frame induced_subgraph
#' make_full_graph permute strength V "V<-" vcount as_edgelist
#' @importFrom matrixStats rowSds
#' @importFrom methods is
#' @importFrom plotly ggplotly plotlyOutput renderPlotly plot_ly layout add_trace
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rintrojs introjs introjsUI
#' @importFrom rlang .data
#' @importFrom rmarkdown render
#' @importFrom S4Vectors mcols
#' @importFrom scales alpha muted
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @importFrom shinyAce aceEditor
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyWidgets dropdownButton tooltipOptions
#' @importFrom stats var dist hclust as.dendrogram as.dist cmdscale median
#' order.dendrogram runif na.omit
#' @import SummarizedExperiment
#' @importFrom tidyr separate_rows pivot_longer
#' @importFrom tippy with_tippy
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils read.delim sessionInfo browseURL citation data write.table
#' head packageDescription
#' @importFrom visNetwork renderVisNetwork visExport visIgraph visNetworkOutput
#' visOptions
#' @importFrom viridis viridis
#' @importFrom mosdef deresult_to_df gene_plot geneinfo_to_html go_to_html
#' map_to_color styleColorBar_divergent create_link_GO create_link_NCBI
#' create_link_GeneCards create_link_GO 
#'
#' @name GeneTonic-pkg
#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription("GeneTonic", fields = "Version")
  msg <- paste0("Welcome to GeneTonic v", pkgVersion, "\n\n")
  citation <- paste0(
    "If you use GeneTonic in your work, please cite:\n\n",
    "  GeneTonic: an R/Bioconductor package for streamlining the interpretation of RNA-seq data\n",
    "  Federico Marini, Annekathrin Ludt, Jan Linke, Konstantin Strauch\n",
    "  BMC Bioinformatics, 2021 - https://doi.org/10.1186/s12859-021-04461-5\n",
    "and/or (if adopting the series of protocols as a whole)\n",
    "  Interactive and Reproducible Workflows for Exploring and Modeling RNA-seq Data with pcaExplorer, ideal, and GeneTonic\n",
    "  Annekathrin Ludt, Arsenij Ustjanzew, Harald Binder, Konstantin Strauch, Federico Marini\n",
    "  Current Protocols, 2022 - https://doi.org/10.1002/cpz1.411\n"
  )
  packageStartupMessage(paste0(msg, citation))
}
