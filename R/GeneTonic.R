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

  # dashpage definition -----------------------------------------------------
  genetonic_ui <- bs4Dash::bs4DashPage(
    # enable_preloader = TRUE,
    title = "GeneTonic",
    sidebar_collapsed = TRUE,
    controlbar_collapsed = TRUE,

    # navbar definition -------------------------------------------------------
    navbar = bs4Dash::bs4DashNavbar(
      skin = "light",
      controlbarIcon = "gears",
      fixed = FALSE,
      rightUi = tagList(
        # actionButton(
        #   inputId = "btn_help_navbar",
        #   icon = icon("info-circle"),
        #   label = "Help", style = .actionbutton_biocstyle
        # )
        shinyWidgets::dropdownButton(
          circle = TRUE,
          status = "info",
          icon = icon("question-circle"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs !"),
          tags$h5("Documentation"),
          actionButton(
            inputId = "btn_docs_vignette",
            icon = icon("info-circle"),
            label = "Help", style = .actionbutton_biocstyle
          ),
          actionButton(
            inputId = "btn_docs_link",
            icon = icon("question-circle"),
            label = "Help", style = .actionbutton_biocstyle
          )

        ),
        shinyWidgets::dropdownButton(
          circle = TRUE,
          status = "info",
          icon = icon("info"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs !"),
          tags$h5("Additional information"),
          actionButton(
            inputId = "btn_info_session",
            icon = icon("info-circle"),
            label = "About this session", style = .actionbutton_biocstyle
          ),
          actionButton(
            inputId = "btn_info_gt",
            icon = icon("heart"),
            label = "About GeneTonic", style = .actionbutton_biocstyle
          )

        )
      )


      # title = "TODOtitle",
      # titleWidth = 350
      # bs4Dash::dropdownMenu(
      #   type = "tasks",
      #   icon = icon("question-circle fa-1g"),
      #   badgeStatus = NULL,
      #   headerText = "Documentation"
      #   # ,
      #   # bs4Dash::notificationItem(
      #   #   text = actionButton(
      #   #     "interface_overview", "Overview of the interface",
      #   #     icon("hand-o-right"),
      #   #     style = .actionbutton_biocstyle
      #   #   ),
      #   #   icon = icon(""), # tricking it to not have additional icon
      #   #   status = "primary"
      #   # )
      # )
    ),

    # sidebar definition ------------------------------------------------------
    sidebar = bs4Dash::dashboardSidebar(
      title = HTML("<small>GeneTonic</small>"),
      skin = "dark",
      status = "primary",
      brandColor = NULL,
      url = "http://bioconductor.org/",
      # src = "logos/online-learning.png",
      elevation = 1,
      opacity = 0.8,
      # width = 250,
      # bs4Dash::menuItem(
      #   text = "SomeSettings", icon = icon("cog"),
      #   startExpanded = TRUE,
      #   numericInput(inputId = "n_genesets",
      #                label = "number of genesets",
      #                value = 15, min = 1, max = 50),
      #   uiOutput("ui_exp_condition")
      # )

      bs4SidebarMenu(
        bs4SidebarMenuItem(
          "Welcome!",
          tabName = "tab_welcome",
          icon = "home"
        ),
        bs4SidebarMenuItem(
          "Gene-Geneset",
          tabName = "tab_ggs",
          icon = "share-alt-square"
        ),
        bs4SidebarMenuItem(
          "Enrichment Map",
          tabName = "tab_emap",
          icon = "map"
        ),
        bs4SidebarMenuItem(
          "DEview",
          tabName = "tab_deview",
          icon = "eye"
        ),
        bs4SidebarMenuItem(
          "About",
          tabName = "tab_about",
          icon = "institution"
        )
      )
    ),

    # body definition ---------------------------------------------------------
    body = bs4Dash::bs4DashBody(
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

      bs4TabItems(
        # Network panel
        bs4TabItem(
          tabName = "tab_welcome",
          fluidRow(
            h2("Whatever goes in the home/welcome page"),
            h3("Overview on the provided input objects")),
          fluidRow(
            bs4Dash::bs4Card(width = 6,
                             title = "Expression Matrix",
                             status = "danger",
                             solidHeader = FALSE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             DT::dataTableOutput("overview_dds")
            ),
            bs4Dash::bs4Card(width = 6,
                             title = "DE results",
                             status = "warning",
                             solidHeader = FALSE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             gradientColor = "warning",

                             DT::dataTableOutput("overview_res_de")
            ),

            bs4Dash::bs4Card(width = 6,
                             title = "Functional analysis results",
                             status = "success",
                             solidHeader = TRUE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             DT::dataTableOutput("overview_res_enrich")
            ),
            bs4Dash::bs4Card(width = 6,
                             title = "Annotation info",
                             status = "info",
                             solidHeader = TRUE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             DT::dataTableOutput("overview_annotation")
            )
          )
        ),

        bs4TabItem(
          tabName = "tab_ggs",

          fluidRow(
            column(
              width = 9,
              withSpinner(
                visNetworkOutput("mynetwork", height = "700px", width = "100%")
              )
            ),
            column(
              width = 3,
              box(
                h4("Genesetbox"),
                verbatimTextOutput("netnode"),
                plotOutput("net_sigheatplot")
              ),
              uiOutput("ui_net_geneinfo")

              # TODOTODO
            )
          )
        ),

        bs4TabItem(
          tabName = "tab_emap",
          fluidRow(
            column(
              width = 8,
              withSpinner(
                visNetworkOutput("emap_visnet", height = "700px", width = "100%")
              )
            ),
            column(
              width = 4,
              uiOutput("ui_emap_sidecontent")
            )
          )
        ),

        bs4TabItem(
          tabName = "tab_deview",
          fluidRow(
            plotOutput("gs_volcano"),
            plotOutput("enriched_funcres"),
            plotlyOutput("enriched_funcres_plotly")
          )
        ),

        bs4TabItem(
          tabName = "tab_about",

          fluidRow(
            column(
              width = 8,
              includeMarkdown(system.file("extdata", "about.md",package = "GeneTonic")),

              verbatimTextOutput("sessioninfo")
            )
          )
        )
      )

      # About section Panel



      ## main structure of the body for the dashboard
      # div(
      #   id = "myScrollBox", # trick to have the y direction scrollable
      #   bs4Dash::bs4TabCard(id = "id",
      #     width=12,
      #     # ui panel welcome -----------------------------------------------------------
      #     bs4Dash::tabPanel(
      #       tabName = "Welcome!",  icon = icon("home"), value="tab-welcome",
      #       fluidRow(
      #           h2("Whatever goes in the home/welcome page"),
      #
      #           h3("Overview on the provided input objects")),
      #           fluidRow(
      #           bs4Dash::bs4Card(width = 6,
      #                            title = "Expression Matrix",
      #                            status = "danger",
      #                            solidHeader = FALSE,
      #                            collapsible = TRUE,
      #                            collapsed = TRUE,
      #                            DT::dataTableOutput("overview_dds")
      #           ),
      #           bs4Dash::bs4Card(width = 6,
      #               title = "DE results",
      #               status = "warning",
      #               solidHeader = FALSE,
      #               collapsible = TRUE,
      #               collapsed = TRUE,
      #               gradientColor = "warning",
      #
      #               DT::dataTableOutput("overview_res_de")
      #         ),
      #         column(
      #           width = 12,
      #           bs4Dash::bs4Card(width = 6,
      #               title = "Functional analysis results",
      #               status = "success",
      #               solidHeader = TRUE,
      #               collapsible = TRUE,
      #               collapsed = TRUE,
      #               DT::dataTableOutput("overview_res_enrich")
      #           ),
      #           bs4Dash::bs4Card(width = 6,
      #               title = "Annotation info",
      #               status = "info",
      #               solidHeader = TRUE,
      #               collapsible = TRUE,
      #               collapsed = TRUE,
      #               DT::dataTableOutput("overview_annotation")
      #           )
      #         )
      #       ),
      #
      #     ),
      #
      #
      #     # ui panel geneset-gene ---------------------------------------------------
      #     bs4Dash::tabPanel(
      #       tabName = "GeneSet-Gene",  icon = icon("home"), value="tab-gsg",
      #

      #     )
      #   )
      # )
      # , footer()
    )
  )
  # genetonic_ui <- shinydashboard::dashboardPage(
  #   skin = "black",
  #
  #   # header definition -------------------------------------------------------
  #   header = shinydashboard::dashboardHeader(
  #     title = "TODOtitle",
  #     titleWidth = 350,
  #     shinydashboard::dropdownMenu(
  #       type = "tasks",
  #       icon = icon("question-circle fa-1g"),
  #       badgeStatus = NULL,
  #       headerText = "Documentation",
  #       shinydashboard::notificationItem(
  #         text = actionButton(
  #           "interface_overview", "Overview of the interface",
  #           icon("hand-o-right"),
  #           style = .actionbutton_biocstyle
  #         ),
  #         icon = icon(""), # tricking it to not have additional icon
  #         status = "primary"
  #       )
  #     )
  #   ),
  #
  #   # sidebar definition ------------------------------------------------------
  #   sidebar = shinydashboard::dashboardSidebar(
  #     width = 250,
  #     shinydashboard::menuItem(
  #       text = "SomeSettings", icon = icon("cog"),
  #       startExpanded = TRUE,
  #       numericInput(inputId = "n_genesets",
  #                    label = "number of genesets",
  #                    value = 15, min = 1, max = 50),
  #       uiOutput("ui_exp_condition")
  #     )
  #   ),
  #
  #   # body definition ---------------------------------------------------------
  #   body = shinydashboard::dashboardBody(
  #     rintrojs::introjsUI(),
  #     ## Define output size and style of error messages
  #     ## plus, define the myscrollbox div to prevent y overflow when page fills up
  #     tags$head(
  #       tags$style(
  #         HTML(
  #           ".shiny-output-error-validation {
  #           font-size: 15px;
  #           color: forestgreen;
  #           text-align: center;
  #           }
  #
  #           #myScrollBox{
  #             overflow-y: scroll;
  #             .dataTables_wrapper{
  #               overflow-x: scroll;
  #             }
  #           }"
  #         )
  #       )
  #       # .icon-done {
  #       # color: green;
  #       # }
  #
  #       # #myAnchorBox{}
  #     ),
  #
  #     ## main structure of the body for the dashboard
  #     div(
  #       id = "myScrollBox", # trick to have the y direction scrollable
  #       tabBox(
  #         width=12,
  #         # ui panel welcome -----------------------------------------------------------
  #         tabPanel(
  #           title = "Welcome!",  icon = icon("home"), value="tab-welcome",
  #           fluidRow(
  #             column(
  #               width = 12,
  #               h2("Whatever goes in the home/welcome page"),
  #
  #               h3("Overview on the provided input objects"),
  #               bs4Dash::bs4Card(width = 6,
  #                   title = "Expression Matrix",
  #                   status = "danger",
  #                   solidHeader = TRUE,
  #                   collapsible = TRUE,
  #                   collapsed = TRUE,
  #                   DT::dataTableOutput("overview_dds")
  #               ),
  #               box(width = 6,
  #                   title = "DE results",
  #                   status = "warning",
  #                   solidHeader = TRUE,
  #                   collapsible = TRUE,
  #                   collapsed = TRUE,
  #                   DT::dataTableOutput("overview_res_de")
  #               )
  #             ),
  #             column(
  #               width = 12,
  #               box(width = 6,
  #                   title = "Functional analysis results",
  #                   status = "success",
  #                   solidHeader = TRUE,
  #                   collapsible = TRUE,
  #                   collapsed = TRUE,
  #                   DT::dataTableOutput("overview_res_enrich")
  #               ),
  #               box(width = 6,
  #                   title = "Annotation info",
  #                   status = "info",
  #                   solidHeader = TRUE,
  #                   collapsible = TRUE,
  #                   collapsed = TRUE,
  #                   DT::dataTableOutput("overview_annotation")
  #               )
  #             )
  #           ),
  #
  #         ),
  #
  #
  #         # ui panel geneset-gene ---------------------------------------------------
  #         tabPanel(
  #           title = "GeneSet-Gene",  icon = icon("home"), value="tab-gsg",
  #
  #           fluidRow(
  #             column(
  #               width = 9,
  #               withSpinner(
  #                 visNetworkOutput("mynetwork", height = "700px", width = "100%")
  #               )
  #             ),
  #             column(
  #               width = 3,
  #               box(
  #                 h4("Genesetbox"),
  #                 verbatimTextOutput("netnode"),
  #                 plotOutput("net_sigheatplot")
  #               ),
  #               uiOutput("ui_net_geneinfo")
  #
  #               # TODOTODO
  #             )
  #           )
  #         ),
  #
  #         # ui panel enrichment map -------------------------------------------------
  #         tabPanel(
  #           title = "Enrichment map",  icon = icon("home"), value="tab-em",
  #           ###
  #           fluidRow(
  #             column(
  #               width = 8,
  #               withSpinner(
  #                 visNetworkOutput("emap_visnet", height = "700px", width = "100%")
  #               )
  #             ),
  #             column(
  #               width = 4,
  #               uiOutput("ui_emap_sidecontent")
  #             )
  #           )
  #
  #         ),
  #
  #         # ui panel de view --------------------------------------------------------
  #         tabPanel(
  #           title = "DEview!",  icon = icon("home"), value="tab-deview",
  #           ###
  #           fluidRow(
  #             plotOutput("gs_volcano"),
  #             plotOutput("enriched_funcres"),
  #             plotlyOutput("enriched_funcres_plotly")
  #           )
  #         ),
  #
  #         # ui panel about -----------------------------------------------------------
  #         tabPanel(
  #           title = "About", icon = icon("institution"), value="tab-about",
  #
  #           fluidRow(
  #             column(
  #               width = 8,
  #               includeMarkdown(system.file("extdata", "about.md",package = "GeneTonic")),
  #
  #               verbatimTextOutput("sessioninfo")
  #             )
  #           )
  #         )
  #       )
  #     )
  #     # , footer()
  #   )
  # )

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


    # panel Welcome -----------------------------------------------------------

    output$overview_dds <- DT::renderDataTable({
      DT::datatable(
        counts(dds),
        options = list(scrollX = TRUE,scrollY = "400px")
      )
    })
    output$overview_res_de <- DT::renderDataTable({
      DT::datatable(
        as.data.frame(res_de),
        options = list(scrollX = TRUE,scrollY = "400px")
      )
    })
    output$overview_res_enrich <- DT::renderDataTable({
      DT::datatable(
        res_enrich,
        options = list(scrollX = TRUE,scrollY = "400px")
      )
    })
    output$overview_annotation <- DT::renderDataTable({
      DT::datatable(
        annotation_obj,
        options = list(scrollX = TRUE,scrollY = "400px")
      )
    })


    # panel GeneSet-Gene ------------------------------------------------------
    values$mygraph <- reactive({
      g <- ggs_graph(res_enrich = res_enrich,
                     res_de = res_de,
                     annotation_obj = annotation_obj,
                     n_gs = input$n_genesets,
                     genes_colname = "genes",
                     genesetname_colname = "Term",
                     genesetid_colname = "GO.ID",
                     prettify = TRUE,
                     geneset_graph_color = "gold")
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
                 annotation_obj = annotation_obj,
                 geneset_id = cur_gsid, # TODOTODO check that I select a gene set
                 genes_colname = "genes",
                 genesetname_colname = "Term",
                 genesetid_colname = "GO.ID",
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

    output$ui_net_geneinfo <- renderUI({
      tagList(
        uiOutput("ggs_gene_info"),
        plotOutput("net_geneplot")
      )
    })

    output$ggs_gene_info <- renderUI({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel,V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "Feature",
                    message = "Please select a gene/feature."
      ))

      cur_geneid <- annotation_obj$gene_id[match(cur_sel, annotation_obj$gene_name)]

      # mycontent <- HTML(paste0(
      #   cur_geneid, "<br>", "<b>", cur_sel, "</b>"
      # ))

      geneinfo_2_html(cur_sel)
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
                    annotation_obj = annotation_obj,
                    n_gs = 50)
    })

    output$gs_volcano <- renderPlot({
      gs_volcano(
        get_aggrscores(res_enrich,
                       res_de,
                       annotation_obj = annotation_obj))
    })

    output$enriched_funcres_plotly <- renderPlotly({
      ggplotly(enhance_table(res_enrich,
                             res_de,
                             annotation_obj = annotation_obj,
                             n_gs = 50))
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

    output$ui_emap_sidecontent <- renderUI({
      tagList(
        plotOutput("emap_sigheatplot"),
        uiOutput("emap_geneset_info")
      )
    })

    output$emap_geneset_info <- renderUI({
      cur_gsid <- res_enrich$GO.ID[match(input$emap_visnet_selected,res_enrich$Term)]
      validate(need(!is.na(cur_gsid),
                    message = "Please select a gene set from the enrichment map."))

      # message(cur_gsid)
      # GOTERM[[cur_gsid]]
      go_2_html(cur_gsid)
    })

    output$emap_sigheatplot <- renderPlot({
      # g <- values$emap_graph()
      # cur_sel <- input$emap_visnet_selected
      # cur_node <- match(cur_sel,V(g)$name)
      # cur_nodetype <- V(g)$nodetype[cur_node]
      # validate(need(cur_nodetype == "GeneSet",
      #               message = "Please select a gene set."
      # ))
      cur_gsid <- res_enrich$GO.ID[match(input$emap_visnet_selected,res_enrich$Term)]
      validate(need(!is.na(cur_gsid),
                    message = "Please select a gene set from the enrichment map."))

      # message(cur_gsid)
      gs_heatmap(myvst,
                 res_de,
                 res_enrich,
                 annotation_obj = annotation_obj,
                 geneset_id = cur_gsid, # TODOTODO check that I select a gene set
                 genes_colname = "genes",
                 genesetname_colname = "Term",
                 genesetid_colname = "GO.ID",
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
