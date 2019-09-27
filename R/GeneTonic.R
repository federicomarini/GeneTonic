#' GeneTonic
#'
#' GeneTonic, main function for the Shiny app
#'
#' @param dds dds
#' @param res_de res_de
#' @param res_enrich res_enrich
#'
#' @return A Shiny app object is returned, for interactive data exploration
#' @export
#'
#' @author Federico Marini
#'
#' @examples
#' library(GeneTonic)
#' # whatever comes next
GeneTonic <- function(dds,
                      res_de,
                      res_enrich) {
  # checks on the objects provided


  # main body of the Shiny app function!
  genetonic_ui <- shinydashboard::dashboardPage(
    skin = "black",

    header = shinydashboard::dashboardHeader(
      title = "TODOtitle",
      titleWidth = 350,
      shinydashboard::dropdownMenu(
        type = "tasks",
        icon = icon("question-circle fa-1g"),
        badgeStatus = NULL,
        headerText = "Documentation",
        shinydashboard::notificationItem(
          text = actionButton(
            "interface_overview", "Overview of the interface",
            icon("hand-o-right")
          ),
          icon = icon(""), # tricking it to not have additional icon
          status = "primary"
        )
      )
    ),

    sidebar = shinydashboard::dashboardSidebar(
      width = 250,
      shinydashboard::menuItem(
        text = "SomeSettings", icon = icon("cog"),
        startExpanded = TRUE,
        numericInput(inputId = "n_genesets",
                     label = "number of genesets",
                     value = 15, min = 1, max = 50)
      )
    ),

    body = shinydashboard::dashboardBody(
      rintrojs::introjsUI(),

      fluidRow(
        h1("The main content"),
        visNetworkOutput("mynetwork", height = "700px", width = "100%")
      )

    )
  )

  options(shiny.maxRequestSize = 15*1024^2)

  #nocov start
  genetonic_server <- function(input, output, session) {

    values <- reactiveValues()

    values$mygraph <- reactive({
      enrich2graph(res_enrich = res_enrich,
                                   res_de = res_de,
                                   n_nodes = input$n_genesets,
                                   genes_colname = "genes",
                                   genesetname_colname = "Term",
                                   genesetid_colname = "GO.ID",
                                   prettify = TRUE,
                                   geneset_graph_color = "gold",
                                   annotation_obj = example_anno_df)
    })

    output$mynetwork <- renderVisNetwork({
      # minimal example

      visNetwork::visIgraph(values$mygraph()) %>%
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), nodesIdSelection = TRUE)

    })

    observeEvent(input$interface_overview, {
      tour <- read.delim(system.file("extdata", "interface_overview.txt",
                                     package = "GeneTonic"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = "")
      rintrojs::introjs(session, options = list(steps = tour))
    })

  }
  #nocov end

  shinyApp(ui = genetonic_ui, server = genetonic_server)
}
