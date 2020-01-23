#' GeneTonic
#'
#' GeneTonic, main function for the Shiny app
#'
#' @param dds A `DESeqDataSet` object, normally obtained after running your data
#' through the `DESeq2` framework.
#' @param res_de A `DESeqResults` object. As for the `dds` parameter, this is
#' also commonly used in the `DESeq2` framework.
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. Required columns for enjoying the full functionality of
#' [GeneTonic()] include:
#' - a gene set identifier (e.g. GeneOntology id, `gs_id`) and its term description
#' (`gs_description`)
#' - a numeric value for the significance of the enrichment (`gs_pvalue`)
#' - a column named `gs_genes` containing a comma separated vector of the gene names
#' associated to the term, one for each term
#' - the number of genes in the geneset of interest detected as differentially
#' expressed (`gs_de_count`), or in the background set of genes (`gs_bg_count`)
#' See [shake_topGOtableResult()] or [shake_enrichResult()] for examples of such
#' formatting helpers
#' @param annotation_obj A `data.frame` object, containing two columns, `gene_id`
#' with a set of unambiguous identifiers (e.g. ENSEMBL ids) and `gene_name`,
#' containing e.g. HGNC-based gene symbols. This object can be constructed via
#' the `org.eg.XX.db` packages, e.g. with convenience functions such as
#' [pcaExplorer::get_annotation_orgdb()].
#' @param project_id  A character string, which can be considered as an identifier
#' for the set/session, and will be e.g. used in the title of the report created
#' via [happy_hour()]
#'
#' @return A Shiny app object is returned, for interactive data exploration
#' @export
#'
#' @author Federico Marini
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'                      keys = rownames(dds_macrophage),
#'                      column = "SYMBOL",
#'                      keytype = "ENSEMBL"),
#'   stringsAsFactors = FALSE,
#'   row.names = rownames(dds_macrophage)
#' )
#'
#' # res object
#' data(res_de_macrophage, package = "GeneTonic")
#' res_de <- res_macrophage_IFNg_vs_naive
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#' res_enrich <- get_aggrscores(res_enrich, res_de, anno_df)
#'
#' # now everything is in place to launch the app
#' if (interactive())
#'   GeneTonic(dds = dds_macrophage,
#'             res_de = res_de,
#'             res_enrich = res_enrich,
#'             annotation_obj = anno_df,
#'             project_id = "myexample")
#'
GeneTonic <- function(dds,
                      res_de,
                      res_enrich,
                      annotation_obj,
                      project_id = "") {

  options(spinner.type = 6)
  # https://projects.lukehaas.me/css-loaders/
  # or even think of https://cran.r-project.org/web/packages/shinycustomloader/README.html
  options(spinner.color = .biocgreen)

  usage_mode <- "shiny_mode"


  # checks on the objects provided
  checkup_GeneTonic(dds,
                    res_de,
                    res_enrich,
                    annotation_obj)

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
      skin = "dark",
      controlbarIcon = "gears",
      fixed = TRUE,
      leftUi = tagList(
        tags$code(tags$h3("GeneTonic")),
        actionButton("bookmarker", label = "Bookmark", icon = icon("heart"),
                     style = "color: #ffffff; background-color: #ac0000; border-color: #ffffff", class = "ml-5")
      ),
      rightUi = tagList(
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_docs",
          circle = FALSE,
          status = "info",
          icon = icon("book"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "More documentation"),
          tags$h5("Documentation"),
          actionButton(
            inputId = "btn_docs_vignette",
            icon = icon("book-open"),
            label = "Open GeneTonic Vignette", style = .actionbutton_biocstyle,
            onclick = ifelse(system.file("doc", "GeneTonic_manual.html", package = "GeneTonic") != "",
                             "",
                             "window.open('https://federicomarini.github.io/GeneTonic/articles/GeneTonic_manual.html', '_blank')")
                             # sprintf("window.open('http://bioconductor.org/packages/%s/bioc/vignettes/GeneTonic/inst/doc/GeneTonic_manual.html', '_blank')",
                             #         ifelse(unlist(packageVersion("GeneTonic"))[2] %% 2L==0L, "release", "devel")
                             # )
            # )
          ),
          actionButton(
            inputId = "btn_docs_link",
            icon = icon("question-circle"),
            label = "Help", style = .actionbutton_biocstyle
          )

        ),
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_info",
          circle = FALSE,
          status = "info",
          icon = icon("info"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "More info!"),
          tags$h5("Additional information"),
          actionButton(
            inputId = "btn_info_session",
            icon = icon("info-circle"),
            label = "About this session", style = .actionbutton_biocstyle
          ),
          actionButton(
            inputId = "btn_info_genetonic",
            icon = icon("heart"),
            label = "About GeneTonic", style = .actionbutton_biocstyle
          )
        )
      )
    ),

    # sidebar definition ------------------------------------------------------
    sidebar = bs4Dash::bs4DashSidebar(
      title = HTML("<small>GeneTonic</small>"),
      src = "GeneTonic/GeneTonic.png",
      skin = "dark",
      status = "primary",
      brandColor = NULL,
      url = "http://bioconductor.org/",
      # src = "logos/online-learning.png",
      elevation = 1,
      opacity = 0.8,

      bs4SidebarMenu(
        id = "gt_tabs",
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
          icon = "project-diagram" # hubspot? map?
        ),
        bs4SidebarMenuItem(
          "Overview",
          tabName = "tab_overview",
          icon = "eye"
        ),
        bs4SidebarMenuItem(
          "GSViz",
          tabName = "tab_gsviz",
          icon = "images"
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
            }
            "
          )
        )
      ),

      tags$head(
        tags$style(
          ".biocdlbutton{background-color:#0092AC;} .biocdlbutton{color: #ffffff;}"
        )
      ),

      tags$script(HTML("$(function(){
      $(document).keyup(function(e) {
      if (e.which == 17) {
        $('#bookmarker').click()
      }
      });
      })")),

      # 27: esc, works
      # 60, <, works NOT1
      # 17, ctrl left, works

      # see more here:
      # https://stackoverflow.com/questions/41675059/keyboard-shortcuts-to-trigger-reactive-flows-in-r-shiny
      # https://stackoverflow.com/questions/10655202/detect-multiple-keys-on-single-keypress-event-in-jquery
      # http://keycode.info/

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
            h2("Overview on the provided input")
          ),
          fluidRow(
            bs4Dash::bs4Card(
              width = 6,
              inputId = "card_em",
              title = "Expression Matrix",
              status = "danger",
              solidHeader = FALSE,
              collapsible = TRUE,
              collapsed = TRUE,
              closable = FALSE,
              DT::dataTableOutput("overview_dds"
              )
            ),
            bs4Dash::bs4Card(
              width = 6,
              inputId = "card_de",
              title = "DE results",
              status = "warning",
              solidHeader = FALSE,
              collapsible = TRUE,
              collapsed = TRUE,
              closable = FALSE,
              DT::dataTableOutput("overview_res_de")
            ),
            bs4Dash::bs4Card(
              width = 6,
              inputId = "card_enrich",
              title = "Functional analysis results",
              status = "success",
              solidHeader = FALSE,
              collapsible = TRUE,
              collapsed = TRUE,
              closable = FALSE,
              DT::dataTableOutput("overview_res_enrich")
            ),
            bs4Dash::bs4Card(
              width = 6,
              inputId = "card_anno",
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
              width = 8,
              withSpinner(
                visNetworkOutput("ggsnetwork", height = "700px", width = "100%")
              )
            ),
            column(
              width = 4,
              bs4Card(
                width = 12,
                uiOutput("ui_ggs_genesetbox")
              ),
              # box(),
              hr(),
              bs4Card(
                width = 12,
                uiOutput("ui_ggs_genebox")
              )
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

        # ui panel overview --------------------------------------------------------
        bs4TabItem(
          tabName = "tab_overview",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_overview", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),

          fluidRow(
            bs4Dash::column(
              width = 10,
              offset = 1,
              bs4Dash::bs4TabCard(
                id = "tabcard_deview",
                title = "Different views",
                side = "right",
                elevation = 2,
                width = 12,
                bs4TabPanel(
                  tabName = "Tab 1",
                  active = TRUE,
                  withSpinner(plotOutput("gs_volcano"))
                ),
                bs4TabPanel(
                  tabName = "Tab 2",
                  active = FALSE,
                  withSpinner(plotOutput("gs_volcano_simplified"))
                ),
                bs4TabPanel(
                  tabName = "Tab 3",
                  active = FALSE,
                  withSpinner(plotOutput("enriched_funcres"))
                ),
                bs4TabPanel(
                  tabName = "Tab 4",
                  active = FALSE,
                  withSpinner(plotlyOutput("enriched_funcres_plotly"))
                )
              )
            )
          )
        ),

        # ui panel GSViz view ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_gsviz",
          fluidRow(
            column(
              width = 11
            ),
            column(
              width = 1,
              actionButton(
                "tour_gsviz", label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              )
            )
          ),

          fluidRow(
            bs4Dash::column(
              width = 10,
              offset = 1,
              bs4Dash::bs4TabCard(
                id = "tabcard_genesets",
                title = "Different views",
                side = "right",
                elevation = 2,
                width = 12,
                bs4TabPanel(
                  tabName = "Tab 1",
                  active = TRUE,
                  withSpinner(plotOutput("gsscores_heatmap"))

                ),
                bs4TabPanel(
                  tabName = "Tab 2",
                  active = FALSE,
                  withSpinner(plotlyOutput("alluvial_genesets"))
                ),
                bs4TabPanel(
                  tabName = "Tab 3",
                  active = FALSE,
                  withSpinner(plotOutput("gs_summaryheat"))
                ),
                bs4TabPanel(
                  tabName = "Tab 4",
                  active = FALSE,
                  withSpinner(plotOutput("mds_genesets"))
                ),
                bs4TabPanel(
                  tabName = "Tab 5",
                  active = FALSE,
                  withSpinner(plotOutput("gs_summaryoverview"))
                ),
                bs4TabPanel(
                  tabName = "Tab 6",
                  active = FALSE,
                  plotOutput("gs_summary_overview_pair")
                ),
                bs4TabPanel(
                  tabName = "Tab 7",
                  active = FALSE,
                  withSpinner(plotOutput("gs_summaryhorizon"))
                ),
                bs4TabPanel(
                  tabName = "Tab 8",
                  active = FALSE,
                  withSpinner(plotlyOutput("gs_summaryradar"))
                ),
                bs4TabPanel(
                  tabName = "Tab 9",
                  active = FALSE,
                  withSpinner(plotOutput("gs_dendro"))
                )

              )
            )
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
          ),
          fluidRow(
            bs4Dash::column(
              width = 8,
              offset = 2,
              gt_downloadButton(
                "start_happyhour",
                "Start the happy hour!",
                class = "biocdlbutton",
                icon = "cocktail") # magic?
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
      numericInput(inputId = "de_fdr",
                   label = "False Discovery Rate (FDR) for DE",
                   value = 0.05, min = 0.0001, max = 1, step = 0.01),
      numericInput(inputId = "n_genesets",
                   label = "Number of genesets",
                   value = 15, min = 1, max = 50),
      selectInput("exp_condition", label = "Group/color by: ",
                  choices = c(NULL, names(colData(dds))), selected = NULL, multiple = TRUE)
    ),

    # footer definition -------------------------------------------------------
    footer = bs4DashFooter(
      GeneTonic_footer,
      right_text = NULL
    )

  )

  options(shiny.maxRequestSize = 15 * 1024^2)

  #nocov start
  genetonic_server <- function(input, output, session) {

    # reactive objects and setup commands -------------------------------------
    reactive_values <- reactiveValues()

    reactive_values$mygenes <- c()
    reactive_values$mygenesets <- c()

    myvst <- vst(dds)

    res_enhanced <- get_aggrscores(res_enrich = res_enrich,
                                   res_de = res_de,
                                   annotation_obj = annotation_obj)

    # output$ui_exp_condition <- renderUI({
      # selectInput("exp_condition", label = "Group/color by: ",
                  # choices = c(NULL, poss_covars), selected = NULL, multiple = TRUE)
    # })


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
        options = list(
          scrollX = TRUE,
          scrollY = "400px",
          pageLength = 25,
          columnDefs = list(
            list(className = "dt-center", targets = "_all")
          )
        )
      ) %>%
        formatRound(columns = c("log2FoldChange"), digits = 3) %>%
        formatStyle(
          "log2FoldChange",
          background = styleColorBar_divergent(as.data.frame(res_de)$log2FoldChange,
                                               scales::alpha("navyblue", 0.4),
                                               scales::alpha("darkred", 0.4)),
          backgroundSize = "100% 90%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
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
        fluidRow(
          column(
            width = 12,
            bs4ValueBoxOutput("infobox_dds"),
            bs4ValueBoxOutput("infobox_resde"),
            bs4ValueBoxOutput("infobox_resenrich"),
            bs4ValueBoxOutput("infobox_annotation")
          )
          # ,
          # column(
          #   width = 6,
          #   img(src = "GeneTonic/GeneTonic.png", height = "350px")
          # )
        )
      )
    })

    output$infobox_dds <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(nrow(dds), " genes x ", ncol(dds), " samples"),
        subtitle = "dds object",
        icon = "table",
        status = "danger",
        width = NULL
      )
    })

    output$infobox_resde <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(
          nrow(deseqresult2df(res_de, FDR = input$de_fdr)),
          " DE genes"
        ),
        subtitle = "res object",
        icon = "vial",
        status = "warning",
        width = NULL
      )
    })

    output$infobox_resenrich <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(
          nrow(res_enrich),
          " functional categories"
        ),
        subtitle = "func enrich object",
        icon = "share-alt",
        status = "success",
        width = NULL
      )
    })

    output$infobox_annotation <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(ncol(annotation_obj), " feature identifiers for ", nrow(dds), " features"),
        subtitle = "annotation object",
        icon = "table",
        status = "info",
        width = NULL
      )
    })

    # panel GeneSet-Gene ------------------------------------------------------
    reactive_values$ggs_graph <- reactive({
      g <- ggs_graph(
        res_enrich = res_enrich,
        res_de = res_de,
        annotation_obj = annotation_obj,
        n_gs = input$n_genesets,
        prettify = TRUE,
        geneset_graph_color = "gold"
      )
      # rank_gs <- rank(V(g)$name[V(g)$nodetype == "GeneSet"])
      # rank_feats <- rank(V(g)$name[V(g)$nodetype == "Feature"]) +
      #   length(rank_gs) # to keep the GeneSets first
      # g <- permute.vertices(g, c(rank_gs, rank_feats))
      # return(g)
    })

    output$ggsnetwork <- renderVisNetwork({
      # minimal example

      visNetwork::visIgraph(reactive_values$ggs_graph()) %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           degree = 1,
                                           hover = TRUE),
                   nodesIdSelection = TRUE)

    })

    output$netnode <- renderPrint({
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]

      cur_gsid <- res_enrich$gs_id[match(cur_sel, res_enrich$gs_description)]

      paste0("I'm selecting ", input$ggsnetwork_selected, ", which has index ", cur_node, " and is of type ", cur_nodetype, "this is from set", cur_gsid)
    })

    output$ui_ggs_genesetbox <- renderUI({
      tagList(
        h5("Genesetbox"),
        # verbatimTextOutput("netnode"),
        plotOutput("net_sigheatplot"),
        uiOutput("ggs_geneset_info")
      )
    })

    output$net_sigheatplot <- renderPlot({
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "GeneSet",
                    message = "Please select a gene set."
      ))
      cur_gsid <- res_enrich$gs_id[match(input$ggsnetwork_selected, res_enrich$gs_description)]

      if (!is.null(input$exp_condition)) {
        gs_heatmap(
          myvst,
          res_de,
          res_enrich,
          annotation_obj = annotation_obj,
          geneset_id = cur_gsid,
          FDR = input$de_fdr,
          de_only = FALSE,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          center_mean = TRUE,
          scale_row = TRUE,
          anno_col_info = input$exp_condition
        )
      } else {
        gs_heatmap(
          myvst,
          res_de,
          res_enrich,
          annotation_obj = annotation_obj,
          geneset_id = cur_gsid,
          FDR = input$de_fdr,
          de_only = FALSE,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          center_mean = TRUE,
          scale_row = TRUE
        )
      }
    })

    output$ggs_geneset_info <- renderUI({
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "GeneSet",
                    message = "" # "Please select a gene set."
      ))
      cur_gsid <- res_enrich$gs_id[match(input$ggsnetwork_selected, res_enrich$gs_description)]

      go_2_html(cur_gsid, res_enrich)
    })

    output$ui_ggs_genebox <- renderUI({
      tagList(
        h5("Genebox"),
        uiOutput("ggs_gene_info"),
        plotOutput("ggs_geneplot")
      )
    })

    output$ggs_gene_info <- renderUI({
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
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
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "Feature",
                    message = "" # "Please select a gene/feature."
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



    # panel EnrichmentMap -----------------------------------------------------
    emap_graph <- reactive({
      emg <- enrichment_map(
        res_enrich = res_enrich,
        res_de = res_de,
        annotation_obj = annotation_obj,
        n_gs = input$n_genesets,
        overlap_threshold = 0.1,
        scale_edges_width = 200,
        color_by = "gs_pvalue"
      )
      # rank_gs <- rank(V(emg)$name)
      # emg <- permute.vertices(emg, rank_gs)
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
      cur_gsid <- res_enrich$gs_id[match(input$emap_visnet_selected, res_enrich$gs_description)]
      validate(need(!is.na(cur_gsid),
                    message = "Please select a gene set from the enrichment map."))

      # message(cur_gsid)
      # GOTERM[[cur_gsid]]
      go_2_html(cur_gsid, res_enrich)
    })

    output$emap_sigheatplot <- renderPlot({
      # g <- reactive_values$emap_graph()
      # cur_sel <- input$emap_visnet_selected
      # cur_node <- match(cur_sel,V(g)$name)
      # cur_nodetype <- V(g)$nodetype[cur_node]
      # validate(need(cur_nodetype == "GeneSet",
      #               message = "Please select a gene set."
      # ))
      cur_gsid <- res_enrich$gs_id[match(input$emap_visnet_selected, res_enrich$gs_description)]
      validate(need(!is.na(cur_gsid),
                    message = "" # "Please select a gene set from the enrichment map."
                    ))


      if (!is.null(input$exp_condition)) {
        # message(cur_gsid)
        gs_heatmap(
          myvst,
          res_de,
          res_enrich,
          annotation_obj = annotation_obj,
          geneset_id = cur_gsid,
          FDR = input$de_fdr,
          de_only = FALSE,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          center_mean = TRUE,
          scale_row = TRUE,
          anno_col_info = input$exp_condition
        )
      } else {
        gs_heatmap(
          myvst,
          res_de,
          res_enrich,
          annotation_obj = annotation_obj,
          geneset_id = cur_gsid,
          FDR = input$de_fdr,
          de_only = FALSE,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          center_mean = TRUE,
          scale_row = TRUE
        )
      }
    })

    # panel DEview ------------------------------------------------------------

    output$enriched_funcres <- renderPlot({
      enhance_table(res_enrich, res_de,
                    annotation_obj = annotation_obj,
                    n_gs = input$n_genesets)
    })

    output$gs_volcano <- renderPlot({
      gs_volcano(
        get_aggrscores(res_enrich,
                       res_de,
                       annotation_obj = annotation_obj),
        volcano_labels = input$n_genesets
      )
    })

    output$gs_volcano_simplified <- renderPlot({
      gs_volcano(
        get_aggrscores(gs_simplify(res_enrich, gs_overlap = 0.6),
                       res_de,
                       annotation_obj = annotation_obj),
        volcano_labels = input$n_genesets
      )
    })

    output$enriched_funcres_plotly <- renderPlotly({
      ggplotly(enhance_table(res_enrich,
                             res_de,
                             annotation_obj = annotation_obj,
                             n_gs = input$n_genesets))
    })


    # panel genesets view -----------------------------------------------------

    gss_mat <- reactive({
      gs_scores(se = myvst,
                res_de = res_de,
                res_enrich = res_enrich,
                annotation_obj = annotation_obj)
    })
    output$gsscores_heatmap <- renderPlot({
      gs_scoresheat(
        gss_mat(),
        n_gs = input$n_genesets)
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

    output$gs_summaryoverview <- renderPlot({
      gs_summary_overview(res_enrich = res_enhanced,
                          n_gs = input$n_genesets)
    })

    output$gs_summaryoverview_pair <- renderPlot({
      gs_summary_overview_pair(res_enrich = res_enhanced,
                               n_gs = input$n_genesets)
    })

    output$gs_summaryhorizon <- renderPlot({
      gs_horizon(res_enrich = res_enhanced,
                 n_gs = input$n_genesets)
    })

    output$gs_summaryradar <- renderPlotly({
      gs_radar(res_enrich = res_enhanced,
               n_gs = input$n_genesets)
    })

    output$gs_dendro <- renderPlot({
      gs_dendro(res_enrich = res_enhanced,
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
            DT::dataTableOutput("bookmarks_genes"),
            downloadButton("btn_export_genes", label = "", class = "biocdlbutton")
            # ideally completed by a function/param to upload them

          ),
          column(
            width = 6,
            bs4InfoBoxOutput("infobox_book_genesets"),
            h5("Bookmarked genesets"),
            DT::dataTableOutput("bookmarks_genesets"),
            downloadButton("btn_export_genesets", label = "", class = "biocdlbutton")

          )
        )
      )
    })

    output$infobox_book_genes <- renderbs4InfoBox({
      bs4InfoBox(
        title = "Bookmarked genes",
        value = length(reactive_values$mygenes),
        icon = "bookmark",
        status = "info",
        width = 12
      )
    })

    output$infobox_book_genesets <- renderbs4InfoBox({
      bs4InfoBox(
        title = "Bookmarked genesets",
        value = length(reactive_values$mygenesets),
        icon = "bookmark",
        status = "success",
        width = 12
      )
    })

    output$bookmarks_genes <- DT::renderDataTable({
      book_df_genes <- annotation_obj[reactive_values$mygenes, ]
      datatable(book_df_genes, rownames = FALSE)
    })
    output$bookmarks_genesets <- DT::renderDataTable({
      book_df_genesets <- res_enrich[reactive_values$mygenesets, 1:2]
      datatable(book_df_genesets, rownames = FALSE)
    })

    output$btn_export_genes <- downloadHandler(
      filename = function() {
        paste0("GeneTonicBookmarks_genes_", project_id, "_", gsub(" ", "_", gsub("-", "", gsub(":", "-", as.character(Sys.time())))), ".txt")
      }, content = function(file) {
        writeLines(text = reactive_values$mygenes,
                   con = file)
      }
    )

    output$btn_export_genesets <- downloadHandler(
      filename = function() {
        paste0("GeneTonicBookmarks_genesets_", project_id, "_", gsub(" ", "_", gsub("-", "", gsub(":", "-", as.character(Sys.time())))), ".txt")
      }, content = function(file) {
        writeLines(text = reactive_values$mygenesets,
                   con = file)
      }
    )

    output$start_happyhour <- downloadHandler(
      filename = paste0(
        Sys.Date(),
        "_", round(runif(1) * 100), # for not having all w the same name
        "_GeneTonicReport.html"),
      content = function(file) {
        # temporarily switch to the temp dir, in case you do not have write permission to the current working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        # cat(tmp_content,file="GeneTonic_tempreport.Rmd",sep="\n")
        withProgress(rmarkdown::render(
          input = system.file("extdata", "cocktail_recipe.Rmd", package = "GeneTonic"),
          output_file = file,
          # fragment.only = TRUE,
          quiet = TRUE),
          message = "Generating the html report",
          detail = "This can take some time")
      }
    )

    output$sessioninfo <- renderPrint({
      sessionInfo()
    })


    # observers ---------------------------------------------------------------

    observeEvent(input$tour_firststeps, {
      tour <- read.delim(system.file("extdata", "tour_welcome.txt", package = "GeneTonic"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = "")
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_ggs, {
      tour <- read.delim(system.file("extdata", "tour_ggs.txt", package = "GeneTonic"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = "")
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_emap, {
      tour <- read.delim(system.file("extdata", "tour_emap.txt", package = "GeneTonic"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = "")
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_overview, {
      tour <- read.delim(system.file("extdata", "tour_overview.txt", package = "GeneTonic"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = "")
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_gsviz, {
      tour <- read.delim(system.file("extdata", "tour_gsviz.txt", package = "GeneTonic"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = "")
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_bookmarks, {
      tour <- read.delim(system.file("extdata", "tour_bookmarks.txt", package = "GeneTonic"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = "")
      rintrojs::introjs(session, options = list(steps = tour))
    })


    # observe({
    #   print(input$gt_tabs)
    # })

    observeEvent(input$btn_docs_vignette, {
      path <- system.file("doc", "GeneTonic_manual.html", package = "GeneTonic")
      if (path == "") {
        showNotification("This vignette has not been built on this system - Opening the online documentation. Please note that the versions might not be coincident!", type = "warning")
      } else {
        browseURL(path)
      }
    })


    observeEvent(input$btn_info_session, {
      showModal(modalDialog(
        title = "Session information", size = "l", fade = TRUE,
        footer = NULL, easyClose = TRUE,
        tagList(renderPrint({
          sessionInfo()
        }))
      ))
    })

    observeEvent(input$btn_info_genetonic, {
      showModal(modalDialog(
        title = "About GeneTonic", size = "l", fade = TRUE,
        footer = NULL, easyClose = TRUE,
        tagList(
          "GeneTonic_info", br(), br(),
          HTML("If you use this package, please use the following citation information:"),
          renderPrint({
            citation("GeneTonic")
          })
        )
      ))
    })


    # bookmarker --------------------------------------------------------------
    observeEvent(input$bookmarker, {
      if (input$gt_tabs == "tab_welcome")
        showNotification("welcome on board!")
      else if (input$gt_tabs == "tab_ggs") {
        showNotification("ggs baby")
        g <- reactive_values$ggs_graph()
        cur_sel <- input$ggsnetwork_selected
        if (cur_sel == "") {
          showNotification("Select a node in the network to bookmark it", type = "warning")
        } else {
          cur_node <- match(cur_sel, V(g)$name)
          cur_nodetype <- V(g)$nodetype[cur_node]

          if (cur_nodetype == "Feature") {
            cur_sel_id <- annotation_obj$gene_id[match(cur_sel, annotation_obj$gene_name)]
            if (cur_sel_id %in% reactive_values$mygenes) {
              showNotification(sprintf("The selected gene %s (%s) is already in the set of the bookmarked genes.", cur_sel, cur_sel_id), type = "default")
            } else {
              reactive_values$mygenes <- unique(c(reactive_values$mygenes, cur_sel_id))
              # message("there go your genes... ", reactive_values$mygenes)
              showNotification(sprintf("Added %s (%s) to the bookmarked genes. The list contains now %d elements", cur_sel, cur_sel_id, length(reactive_values$mygenes)), type = "message")

            }
          } else if (cur_nodetype == "GeneSet") {
            cur_sel_id <- res_enrich$gs_id[match(cur_sel, res_enrich$gs_description)]
            if (cur_sel_id %in% reactive_values$mygenesets) {
              showNotification(sprintf("The selected gene set %s (%s) is already in the set of the bookmarked genesets.", cur_sel, cur_sel_id), type = "default")
            } else {
              reactive_values$mygenesets <- unique(c(reactive_values$mygenesets, cur_sel_id))
              # message("here are your genesets... ", reactive_values$mygenesets)
              showNotification(sprintf("Added %s (%s) to the bookmarked genesets. The list contains now %d elements", cur_sel, cur_sel_id, length(reactive_values$mygenesets)), type = "message")
            }

          } else {
            message("bleeee")
          }
        }
      }
      else if (input$gt_tabs == "tab_emap") {
        showNotification("maaap maaap")
        g <- reactive_values$ggs_graph()
        cur_sel <- input$emap_visnet_selected
        cur_sel_id <- res_enrich$gs_id[match(cur_sel, res_enrich$gs_description)]

        if (cur_sel == "") {
          showNotification("Select a node in the network to bookmark it", type = "warning")
        } else {
          if (cur_sel_id %in% reactive_values$mygenesets) {
            showNotification(sprintf("The selected gene set %s (%s) is already in the set of the bookmarked genesets.", cur_sel, cur_sel_id), type = "default")
          } else {
            reactive_values$mygenesets <- unique(c(reactive_values$mygenesets, cur_sel_id))
            message("here are your genesets... ", reactive_values$mygenesets)
            showNotification(sprintf("Added %s (%s) to the bookmarked genesets. The list contains now %d elements", cur_sel, cur_sel_id, length(reactive_values$mygenesets)), type = "message")
          }
        }
      }
      else if (input$gt_tabs == "tab_bookmarks")
        showNotification("catching up and wrapping up...")
      else if (input$gt_tabs == "tab_about")
        showNotification("youwannaknowwhatitsallabout")

    })


    # observeEvent(input$start_happyhour, {
    #   showNotification("The happy hour is officially on!",
    #                    type = "message")
    # })



  }
  #nocov end

  shinyApp(ui = genetonic_ui, server = genetonic_server)
}
