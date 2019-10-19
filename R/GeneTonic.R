#' GeneTonic
#'
#' GeneTonic, main function for the Shiny app
#'
#' @param dds dds
#' @param res_de res_de
#' @param res_enrich res_enrich
#' @param annotation_obj anno_obj
#'
#' @return A Shiny app object is returned, for interactive data exploration
#' @export
#'
#' @author Federico Marini
#'
#' @examples
#' # library(GeneTonic)
#' # whatever comes next
GeneTonic <- function(dds,
                      res_de,
                      res_enrich,
                      annotation_obj) {

  options(spinner.type = 6)
  # https://projects.lukehaas.me/css-loaders/
  # or even think of https://cran.r-project.org/web/packages/shinycustomloader/README.html
  options(spinner.color = .biocgreen)


  # checks on the objects provided

  # clean up the result object, e.g. removing the NAs in the relevant columns
  res_de <- res_de[!is.na(res_de$log2FoldChange),]
  message("Removing ", sum(is.na(res_de$log2FoldChange)), " rows from the result object - logFC detected as NA")



  # UI definition -----------------------------------------------------------
  genetonic_ui <- shinydashboard::dashboardPage(
    skin = "black",

    # header definition -------------------------------------------------------
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
            icon("hand-o-right"),
            style = .actionbutton_biocstyle
          ),
          icon = icon(""), # tricking it to not have additional icon
          status = "primary"
        )
      )
    ),

    # sidebar definition ------------------------------------------------------
    sidebar = shinydashboard::dashboardSidebar(
      width = 250,
      shinydashboard::menuItem(
        text = "SomeSettings", icon = icon("cog"),
        startExpanded = TRUE,
        numericInput(inputId = "n_genesets",
                     label = "number of genesets",
                     value = 15, min = 1, max = 50),
        uiOutput("ui_exp_condition")
      )
    ),

    # body definition ---------------------------------------------------------
    body = shinydashboard::dashboardBody(
      rintrojs::introjsUI(),
      ## Define output size and style of error messages
      ## plus, define the myscrollbox div to prevent y overflow when page fills up
      tags$head(
        tags$style(
          HTML(
            ".shiny-output-error-validation {
            font-size: 15px;
            color: forestgreen;
            text-align: center;
            }

            #myScrollBox{
              overflow-y: scroll;
              .dataTables_wrapper{
                overflow-x: scroll;
              }
            }"
          )
        )
         # .icon-done {
         # color: green;
         # }

         # #myAnchorBox{}
      ),

      ## main structure of the body for the dashboard
      div(
        id = "myScrollBox", # trick to have the y direction scrollable
        tabBox(
          width=12,

          # ui panel welcome -----------------------------------------------------------
          tabPanel(
            title = "Welcome!",  icon = icon("home"), value="tab-welcome",
            h2("Whatever goes in the home/welcome page")
          ),


          # ui panel geneset-gene ---------------------------------------------------
          tabPanel(
            title = "GeneSet-Gene",  icon = icon("home"), value="tab-gsg",

            fluidRow(
              column(
                width = 9,
                withSpinner(
                  visNetworkOutput("mynetwork", height = "700px", width = "100%")
                )
              ),
              column(
                width = 3,
                h4("Genesetbox"),
                verbatimTextOutput("netnode"),
                plotOutput("net_sigheatplot"),
                plotOutput("net_geneplot"),
                # TODOTODO
              )
            )
          ),

          # ui panel enrichment map -------------------------------------------------
          tabPanel(
            title = "Enrichment map",  icon = icon("home"), value="tab-em",
            ###
            visNetworkOutput("emap_visnet", height = "700px", width = "100%")
          ),

          # ui panel de view --------------------------------------------------------
          tabPanel(
            title = "DEview!",  icon = icon("home"), value="tab-deview",
            ###
            fluidRow(
              plotOutput("go_volcano"),
              plotOutput("enriched_funcres"),
              plotlyOutput("enriched_funcres_plotly")
            )
          ),

          # ui panel about -----------------------------------------------------------
          tabPanel(
            title = "About", icon = icon("institution"), value="tab-about",

            fluidRow(
              column(
                width = 8,
                includeMarkdown(system.file("extdata", "about.md",package = "GeneTonic")),

                verbatimTextOutput("sessioninfo")
              )
            )
          )
        )
      )
      # , footer()
    )
  )

  options(shiny.maxRequestSize = 15*1024^2)

  #nocov start
  genetonic_server <- function(input, output, session) {

    # reactive objects and setup commands -------------------------------------
    values <- reactiveValues()

    myvst <- vst(dds)

    output$ui_exp_condition <- renderUI({
      poss_covars <- names(colData(dds))
      selectInput("exp_condition", label = "Group/color by: ",
                  choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
    })


    # panel GeneSet-Gene ------------------------------------------------------
    values$mygraph <- reactive({
      g <- enrich2graph(res_enrich = res_enrich,
                   res_de = res_de,
                   n_nodes = input$n_genesets,
                   genes_colname = "genes",
                   genesetname_colname = "Term",
                   genesetid_colname = "GO.ID",
                   prettify = TRUE,
                   geneset_graph_color = "gold",
                   annotation_obj = annotation_obj)
      rank_gs <- rank(V(g)$name[V(g)$nodetype == "GeneSet"])
      rank_feats <- rank(V(g)$name[V(g)$nodetype == "Feature"]) +
        length(rank_gs) # to keep the GeneSets first
      g <- permute.vertices(g, c(rank_gs, rank_feats))
      return(g)
    })


    output$mynetwork <- renderVisNetwork({
      # minimal example

      visNetwork::visIgraph(values$mygraph()) %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           degree = 1,
                                           hover = TRUE),
                   nodesIdSelection = TRUE)

    })

    output$netnode <- renderPrint({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel,V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]

      cur_gsid <- res_enrich$GO.ID[match(cur_sel,res_enrich$Term)]

      paste0("I'm selecting ",input$mynetwork_selected, ", which has index ", cur_node, " and is of type ", cur_nodetype, "this is from set", cur_gsid)

    })

    output$net_sigheatplot <- renderPlot({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel,V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "GeneSet",
                    message = "Please select a gene set."
      ))
      cur_gsid <- res_enrich$GO.ID[match(input$mynetwork_selected,res_enrich$Term)]
      gs_heatmap(myvst,
                 res_de,
                 res_enrich,
                 geneset_id = cur_gsid, # TODOTODO check that I select a gene set
                 genes_colname = "genes",
                 genesetname_colname = "Term",
                 genesetid_colname = "GO.ID",
                 annotation_obj = annotation_obj,
                 FDR = 0.05,
                 de_only = FALSE,
                 cluster_rows = TRUE, # TODOTODO: options for the heatmap go on left side, as could be common to more!
                 cluster_cols = TRUE,
                 center_mean = TRUE,
                 scale_row = TRUE
                 # TODOTODO: use ellipsis for passing params to pheatmap?
                 # TODOTODO: option to just return the underlying data?s
                 # TODOTODO: options to subset to specific samples?
      )
    })

    output$net_geneplot <- renderPlot({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel,V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "Feature",
                    message = "Please select a gene/feature."
      ))
      validate(need(input$exp_condition != "",
                    message = "Please select a group for the experimental condition."
      ))


      cur_geneid <- annotation_obj$gene_id[match(cur_sel, annotation_obj$gene_name)]
      gene_plot(dds,
                gene = cur_geneid,
                intgroup = input$exp_condition,
                annotation_obj = annotation_obj
      )
    })



    # panel DEview ------------------------------------------------------------
    output$enriched_funcres <- renderPlot({
      enhance_table(res_enrich, res_de,
                    n_gs = 50,
                    annotation_obj = annotation_obj)
    })

    output$go_volcano <- renderPlot({
      go_volcano(get_aggrscores(res_enrich,res_de,annotation_obj = annotation_obj))
    })

    output$enriched_funcres_plotly <- renderPlotly({
      ggplotly(enhance_table(res_enrich, res_de,
                    n_gs = 50,
                    annotation_obj = annotation_obj))
    })


    # panel EnrichmentMap -----------------------------------------------------
    emap_graph <- reactive({
      emg <- enrichment_map(res_enrich = res_enrich,
                            res_de = res_de,
                            annotation_obj = annotation_obj,
                            n_gs = input$n_genesets,
                            overlap_threshold = 0.1,
                            scale_edges_width = 200,
                            color_by = "p.value_elim",
                            genes_colname = "genes",
                            genesetname_colname = "Term",
                            genesetid_colname = "GO.ID")
      rank_gs <- rank(V(emg)$name)
      emg <- permute.vertices(emg, rank_gs)
      return(emg)
    })

    output$emap_visnet <- renderVisNetwork({

      visNetwork::visIgraph(emap_graph()) %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           degree = 1,
                                           hover = TRUE),
                   nodesIdSelection = TRUE)

    })

    output$sessioninfo <- renderPrint({
      sessionInfo()
    })


    # observers ---------------------------------------------------------------


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
