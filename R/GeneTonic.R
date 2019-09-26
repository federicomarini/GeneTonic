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
        startExpanded = TRUE
      )
    ),

    body = shinydashboard::dashboardBody(
      rintrojs::introjsUI(),

      fluidRow(
        "The main content"
      )

    )
  )

  options(shiny.maxRequestSize = 15*1024^2)

  #nocov start
  genetonic_server <- function(input, output, session) {

    values <- reactiveValues()

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
