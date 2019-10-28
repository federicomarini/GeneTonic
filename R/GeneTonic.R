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
  res_de <- res_de[!is.na(res_de$log2FoldChange), ]
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
          "GeneSets",
          tabName = "tab_genesets",
          icon = "list-alt"
        ),
        bs4SidebarMenuItem(
          "Bookmarks",
          tabName = "tab_bookmarks",
          icon = "bookmark"
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
        # ui panel welcome -----------------------------------------------------------
        bs4TabItem(
          tabName = "tab_welcome",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_firststeps", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),
          fluidRow(
            h2("Whatever goes in the home/welcome page"),
            h3("Overview on the provided input objects")
          ),
          fluidRow(
            bs4Dash::bs4Card(width = 6,
                             title = "Expression Matrix",
                             status = "danger",
                             solidHeader = FALSE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             closable = FALSE,
                             DT::dataTableOutput("overview_dds")
            ),
            bs4Dash::bs4Card(width = 6,
                             title = "DE results",
                             status = "warning",
                             solidHeader = FALSE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             closable = FALSE,
                             DT::dataTableOutput("overview_res_de")
            ),
            bs4Dash::bs4Card(width = 6,
                             title = "Functional analysis results",
                             status = "success",
                             solidHeader = FALSE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             closable = FALSE,
                             DT::dataTableOutput("overview_res_enrich")
            ),
            bs4Dash::bs4Card(width = 6,
                             title = "Annotation info",
                             status = "info",
                             solidHeader = FALSE,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             closable = FALSE,
                             DT::dataTableOutput("overview_annotation")
            )
          ),
          uiOutput("ui_infoboxes")
        ),

        # ui panel geneset-gene ---------------------------------------------------
        bs4TabItem(
          tabName = "tab_ggs",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_ggs", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),
          fluidRow(
            column(
              width = 9,
              withSpinner(
                visNetworkOutput("mynetwork", height = "700px", width = "100%")
              ),
              actionButton("bookmark_ggs", label = "Bookmark", icon = icon("heart"),
                           style = "color: #ffffff; background-color: #ac0000; border-color: #ffffff")
            ),
            column(
              width = 3,
              bs4Card(
                width = 12,
                uiOutput("ui_ggs_genesetbox")
              ),
              # box(),
              hr(),
              uiOutput("ui_ggs_genebox")

              # TODOTODO
            )
          )
        ),

        # ui panel enrichment map -------------------------------------------------
        bs4TabItem(
          tabName = "tab_emap",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_emap", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),
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

        # ui panel de view --------------------------------------------------------
        bs4TabItem(
          tabName = "tab_deview",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_deview", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),
          fluidRow(
            plotOutput("gs_volcano"),
            plotOutput("gs_volcano_simplified"),
            plotOutput("enriched_funcres"),
            plotlyOutput("enriched_funcres_plotly")
          )
        ),

        # ui panel genesets view ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_genesets",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_genesetsview", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),
          fluidRow(
            plotOutput("gsscores_heatmap"),
            plotlyOutput("alluvial_genesets"),
            plotOutput("gs_summaryheat"),
            plotOutput("mds_genesets")
          )
        ),

        # ui panel bookmark ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_bookmarks",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_bookmarks", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              uiOutput("ui_bookmarks")
            )
          )
        ),

        # ui panel about -----------------------------------------------------------
        bs4TabItem(
          tabName = "tab_about",

          fluidRow(
            column(
              width = 8,
              includeMarkdown(system.file("extdata", "about.md", package = "GeneTonic")),

              verbatimTextOutput("sessioninfo")
            )
          )
        )
      )


    ),
    # controlbar definition ---------------------------------------------------
    controlbar = bs4Dash::bs4DashControlbar(
      numericInput(inputId = "n_genesets",
                   label = "number of genesets",
                   value = 15, min = 1, max = 50),
      uiOutput("ui_exp_condition")
    ),

    footer = bs4DashFooter(
      GeneTonic:::footer()
    )

  )

      # About section Panel



      ## main structure of the body for the dashboard
      # div(
      #   id = "myScrollBox", # trick to have the y direction scrollable
      #   bs4Dash::bs4TabCard(id = "id",
      #     width=12,
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
      #
      #     bs4Dash::tabPanel(
      #       tabName = "GeneSet-Gene",  icon = icon("home"), value="tab-gsg",
      #

      #     )
      #   )
      # )
      # , footer()

  # genetonic_ui <- shinydashboard::dashboardPage(
  #   skin = "black",
  #
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
  #
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
  #
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
  #
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

  options(shiny.maxRequestSize = 15 * 1024^2)

  #nocov start
  genetonic_server <- function(input, output, session) {

    # reactive objects and setup commands -------------------------------------
    values <- reactiveValues()

    values$mygenes <- c()
    values$mygenesets <- c()

    myvst <- vst(dds)

    output$ui_exp_condition <- renderUI({
      poss_covars <- names(colData(dds))
      selectInput("exp_condition", label = "Group/color by: ",
                  choices = c(NULL, poss_covars), selected = NULL, multiple = TRUE)
    })


    # panel Welcome -----------------------------------------------------------

    output$overview_dds <- DT::renderDataTable({
      DT::datatable(
        counts(dds),
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })
    output$overview_res_de <- DT::renderDataTable({
      DT::datatable(
        as.data.frame(res_de),
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })
    output$overview_res_enrich <- DT::renderDataTable({
      DT::datatable(
        res_enrich,
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })
    output$overview_annotation <- DT::renderDataTable({
      DT::datatable(
        annotation_obj,
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })

    output$ui_infoboxes <- renderUI({
      tagList(
        bs4ValueBoxOutput("infobox_dds"),
        bs4ValueBoxOutput("infobox_resde"),
        bs4ValueBoxOutput("infobox_resenrich"),
        bs4ValueBoxOutput("infobox_annotation")
      )
    })

    output$infobox_dds <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(nrow(dds), " x ", ncol(dds)),
        subtitle = "dds object",
        icon = "table",
        status = "danger"
      )
    })

    output$infobox_resde <- renderbs4ValueBox({
      bs4ValueBox(
        value = nrow(deseqresult2df(res_de, FDR = 0.05)), # TODO: set via widget?
        subtitle = "res object",
        icon = "vial",
        status = "warning"
      )
    })

    output$infobox_resenrich <- renderbs4ValueBox({
      bs4ValueBox(
        value = nrow(res_enrich),
        subtitle = "func enrich object",
        icon = "share-alt",
        status = "success"
      )
    })

    output$infobox_annotation <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(ncol(annotation_obj), " feature identifiers for ", nrow(dds)),
        subtitle = "annotation object",
        icon = "table",
        status = "info"
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
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]

      cur_gsid <- res_enrich$GO.ID[match(cur_sel, res_enrich$Term)]

      paste0("I'm selecting ", input$mynetwork_selected, ", which has index ", cur_node, " and is of type ", cur_nodetype, "this is from set", cur_gsid)

    })

    output$ui_ggs_genesetbox <- renderUI({
      tagList(
        h5("Genesetbox"),
        verbatimTextOutput("netnode"),
        plotOutput("net_sigheatplot"),
        uiOutput("ggs_geneset_info")
      )


    })

    output$net_sigheatplot <- renderPlot({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "GeneSet",
                    message = "Please select a gene set."
      ))
      cur_gsid <- res_enrich$GO.ID[match(input$mynetwork_selected, res_enrich$Term)]
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

    output$ggs_geneset_info <- renderUI({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "GeneSet",
                    message = "Please select a gene set."
      ))
      cur_gsid <- res_enrich$GO.ID[match(input$mynetwork_selected, res_enrich$Term)]

      # message(cur_gsid)
      # GOTERM[[cur_gsid]]
      go_2_html(cur_gsid)
    })

    output$ui_ggs_genebox <- renderUI({
      tagList(
        h5("Genebox"),
        uiOutput("ggs_gene_info"),
        plotOutput("ggs_geneplot")
      )
    })

    output$ggs_gene_info <- renderUI({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
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

    output$ggs_geneplot <- renderPlot({
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
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

    output$gs_volcano_simplified <- renderPlot({
      gs_volcano(
        get_aggrscores(gs_simplify(res_enrich, gs_overlap = 0.6),
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
      cur_gsid <- res_enrich$GO.ID[match(input$emap_visnet_selected, res_enrich$Term)]
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
      cur_gsid <- res_enrich$GO.ID[match(input$emap_visnet_selected, res_enrich$Term)]
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


    # panel genesets view -----------------------------------------------------

    gss_mat <- reactive({
      gs_scores(se = myvst,
                res_de = res_de,
                res_enrich = res_enrich,
                annotation_obj = annotation_obj,
                genes_colname = "genes",
                genesetname_colname = "Term",
                genesetid_colname = "GO.ID")
    })
    output$gsscores_heatmap <- renderPlot({
      gs_ggheatmap(gss_mat())
    })

    output$alluvial_genesets <- renderPlotly({
      gs_alluvial(res_enrich, res_de, annotation_obj, n_gs = input$n_genesets)
    })

    output$mds_genesets <- renderPlot({
      gs_mds(res_enrich, res_de, annotation_obj, mds_colorby = "z_score",
             mds_labels = input$n_genesets)
    })

    output$gs_summaryheat <- renderPlot({
      gs_summary_heat(res_enrich, res_de, annotation_obj,
                      n_gs = input$n_genesets)
    })

    # panel bookmarks ---------------------------------------------------------

    output$ui_bookmarks <- renderUI({
      tagList(
        fluidRow(
          column(
            width = 6,
            bs4InfoBoxOutput("infobox_book_genes"),
            h5("Bookmarked genes"),
            DT::dataTableOutput("bookmarks_genes")
          ),
          column(
            width = 6,
            bs4InfoBoxOutput("infobox_book_genesets"),
            h5("Bookmarked genesets"),
            DT::dataTableOutput("bookmarks_genesets")

          )
        )
      )
    })

    output$infobox_book_genes <- renderbs4InfoBox({
      bs4InfoBox(title = "Bookmarked genes",
                 value = length(values$mygenes),
                 icon = "bookmark",
                 status = "info",
                 width = 12)
    })

    output$infobox_book_genesets <- renderbs4InfoBox({
      bs4InfoBox(title = "Bookmarked genesets",
                 value = length(values$mygenesets),
                 icon = "bookmark",
                 status = "success",
                 width = 12)
    })

    output$bookmarks_genes <- DT::renderDataTable({
      datatable(data.frame(mygenes = values$mygenes))
    })
    output$bookmarks_genesets <- DT::renderDataTable({
      datatable(data.frame(mygenesets = values$mygenesets))
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


    observeEvent(input$bookmark_ggs, {
      g <- values$mygraph()
      cur_sel <- input$mynetwork_selected
      if (cur_sel == "") {
        showNotification("Select a node in the network to bookmark it", type = "warning")
      } else {
        cur_node <- match(cur_sel, V(g)$name)
        cur_nodetype <- V(g)$nodetype[cur_node]

        if (cur_nodetype == "Feature") {
          # TODO: match back to identifier and so
          values$mygenes <- c(values$mygenes, cur_sel)
          message("there go your genes... ", values$mygenes)
          showNotification(sprintf("Added %s to the bookmarked genes. The list contains now %d elements", cur_sel, length(values$mygenes)), type = "message")
        } else if (cur_nodetype == "GeneSet") {
          # TODO: match back to identifier and so
          values$mygenesets <- c(values$mygenesets, cur_sel)
          message("here are your genesets... ", values$mygenesets)
          showNotification(sprintf("Added %s to the bookmarked genesets. The list contains now %d elements", cur_sel, length(values$mygenesets)), type = "message")
        } else {
          message("bleeee")
        }
      }

    })



  }
  #nocov end

  shinyApp(ui = genetonic_ui, server = genetonic_server)
}
