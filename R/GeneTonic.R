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
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param project_id  A character string, which can be considered as an identifier
#' for the set/session, and will be e.g. used in the title of the report created
#' via [happy_hour()]
#' @param size_gtl Numeric value, specifying the maximal size in MB for the 
#' accepted GeneTonicList object - this applies when uploading the dataset at 
#' runtime
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
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'     keys = rownames(dds_macrophage),
#'     column = "SYMBOL",
#'     keytype = "ENSEMBL"
#'   ),
#'   stringsAsFactors = FALSE,
#'   row.names = rownames(dds_macrophage)
#' )
#'
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
#' if (interactive()) {
#'   GeneTonic(
#'     dds = dds_macrophage,
#'     res_de = res_de,
#'     res_enrich = res_enrich,
#'     annotation_obj = anno_df,
#'     project_id = "myexample"
#'   )
#' }
#' # alternatively...
#' gtl_macrophage <- GeneTonic_list(
#'   dds = dds_macrophage,
#'   res_de = res_de,
#'   res_enrich = res_enrich,
#'   annotation_obj = anno_df
#' )
#' # GeneTonic(gtl = gtl_macrophage)
GeneTonic <- function(dds = NULL,
                      res_de = NULL,
                      res_enrich = NULL,
                      annotation_obj = NULL,
                      gtl = NULL,
                      project_id = "",
                      size_gtl = 50) {

  # https://projects.lukehaas.me/css-loaders/
  # or even think of https://cran.r-project.org/web/packages/shinycustomloader/README.html
  oopt <- options(spinner.type = 6, spinner.color = .biocgreen)
  # play nice with other previously chosen options
  on.exit(options(oopt))

  usage_mode <- "shiny_mode"
  
  all_components_provided <- 
    !is.null(dds) & !is.null(res_de) & !is.null(res_enrich) & !is.null(annotation_obj) 
  
  if (all_components_provided) {
    message("Checks, all components TODO")
    checkup_GeneTonic(
      dds,
      res_de,
      res_enrich,
      annotation_obj
    )
  }
  
  if (!is.null(gtl)) {
    message("Checks, gtl TODO")
    
    checkup_gtl(gtl)
  }

  # UI definition -----------------------------------------------------------

  # dashpage definition -----------------------------------------------------
  genetonic_ui <- bs4Dash::bs4DashPage(
    # enable_preloader = TRUE,
    title = "GeneTonic",
    dark = NULL,
    # navbar definition -------------------------------------------------------
    header = bs4Dash::bs4DashNavbar(
      # older leftUi elements
      tagList(
        tags$code(tags$h3("GeneTonic")),
        actionButton("bookmarker",
          label = "Bookmark", icon = icon("heart"),
          style = "color: #ffffff; background-color: #ac0000; border-color: #ffffff", class = "ml-5"
        )
      ),
      # older rightUi elements
      # ideally to be spaced on the full right side...
      tags$span(style = "display:inline-block; width: 30%"),
      tagList(
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
              "window.open('https://federicomarini.github.io/GeneTonic/articles/GeneTonic_manual.html', '_blank')"
            )
            # sprintf("window.open('http://bioconductor.org/packages/%s/bioc/vignettes/GeneTonic/inst/doc/GeneTonic_manual.html', '_blank')",
            #         ifelse(unlist(packageVersion("GeneTonic"))[2] %% 2L==0L, "release", "devel")
            # )
            # )
          ),
          actionButton(
            inputId = "btn_first_help",
            icon = icon("question-circle"),
            label = "First Help", style = .actionbutton_biocstyle
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
        ),
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_savegtl",
          circle = FALSE,
          status = "info",
          icon = icon("save"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "Save the gtl object!"),
          tags$h5("Save data as gtl object"),
          # actionButton(
          #   inputId = "btn_save_gtl",
          #   icon = icon("gift"),
          #   label = "Save as gtl serialized object TODO", style = .actionbutton_biocstyle
          # )
          uiOutput("ui_gtl_download")
        )
      ),
      title = bs4Dash::bs4DashBrand(
        title = HTML("<small>GeneTonic</small>"),
        href = "https://bioconductor.org/packages/GeneTonic",
        # color = "info",
        image = "GeneTonic/GeneTonic.png"
      ),
      skin = "dark",
      status = "gray-dark",
      border = FALSE,
      controlbarIcon = icon("cogs"),
      fixed = TRUE # ,
      # leftUi =
      # rightUi =
    ),

    # sidebar definition ------------------------------------------------------
    sidebar = bs4Dash::bs4DashSidebar(
      title = HTML("<small>GeneTonic</small>"),
      src = "GeneTonic/GeneTonic.png",
      skin = "dark",
      status = "primary",
      brandColor = NULL,
      url = "https://bioconductor.org/packages/GeneTonic",
      collapsed = TRUE,
      # src = "logos/online-learning.png",
      elevation = 1,
      opacity = 0.8,
      bs4SidebarMenu(
        id = "gt_tabs",
        bs4SidebarMenuItem(
          "Welcome!",
          tabName = "tab_welcome",
          icon = icon("home")
        ),
        bs4SidebarMenuItem(
          "Gene-Geneset",
          tabName = "tab_ggs",
          icon = icon("share-alt-square")
        ),
        bs4SidebarMenuItem(
          "Enrichment Map",
          tabName = "tab_emap",
          icon = icon("project-diagram") # hubspot? map?
        ),
        bs4SidebarMenuItem(
          "Overview",
          tabName = "tab_overview",
          icon = icon("eye")
        ),
        bs4SidebarMenuItem(
          "GSViz",
          tabName = "tab_gsviz",
          icon = icon("images")
        ),
        bs4SidebarMenuItem(
          "Bookmarks",
          tabName = "tab_bookmarks",
          icon = icon("bookmark")
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
          tagList(
            uiOutput("ui_uploadgtl"),
            uiOutput("ui_panel_welcome")
          )
          # ,
          # verbatimTextOutput("gtl_described")
        ),

        # ui panel geneset-gene ---------------------------------------------------
        bs4TabItem(
          tabName = "tab_ggs",
          uiOutput("ui_panel_ggs")
        ),

        # ui panel enrichment map -------------------------------------------------
        bs4TabItem(
          tabName = "tab_emap",
          uiOutput("ui_panel_em")
        ),

        # ui panel overview --------------------------------------------------------
        bs4TabItem(
          tabName = "tab_overview",
          uiOutput("ui_panel_overview")
        ),

        # ui panel GSViz ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_gsviz",
          uiOutput("ui_panel_gsviz")
        ),

        # ui panel bookmark ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_bookmarks",
          uiOutput("ui_panel_bookmarks")
        )
      )
    ),
    # controlbar definition ---------------------------------------------------
    controlbar = bs4Dash::bs4DashControlbar(
      collapsed = TRUE,
      uiOutput("ui_controlbar")
    ),

    # footer definition -------------------------------------------------------
    footer = bs4DashFooter(
      left = GeneTonic_footer,
      right = NULL
    )
  )

  options(shiny.maxRequestSize = size_gtl * 1024^2)

  # nocov start
  genetonic_server <- function(input, output, session) {
    
    # initializing reactives --------------------------------------------------
    reactive_values <- reactiveValues()
    
    # for the usage in the bookmarks
    reactive_values$mygenes <- c()
    reactive_values$mygenesets <- c()
    
    if (!is.null(gtl)) {
      message("GeneTonic info: gtl object provided")
      
      checkup_gtl(gtl)
      
      reactive_values$gtl <- gtl
      reactive_values$dds <- gtl$dds
      reactive_values$res_de <- gtl$res_de
      reactive_values$res_enrich <- gtl$res_enrich
      reactive_values$annotation_obj <- gtl$annotation_obj
      
      # clean up the result object, e.g. removing the NAs in the relevant columns
      removed_genes <- is.na(gtl$res_de$log2FoldChange)
      message(
        "Removing ", sum(removed_genes),
        "/", nrow(gtl$res_de), " rows from the DE `res_de` object - log2FC values detected as NA"
      )
      reactive_values$res_de <- gtl$res_de[!removed_genes, ]
      
      reactive_values$myvst <- reactive({
        vst(reactive_values$dds)
      })
      
      reactive_values$res_enhanced <- reactive({
        get_aggrscores(
          res_enrich = reactive_values$res_enrich,
          res_de = reactive_values$res_de,
          annotation_obj = reactive_values$annotation_obj
        )
      })
      
    } else if (all_components_provided){
      message("GeneTonic info: all components provided")
      
      checkup_GeneTonic(dds = dds,
                        res_de = res_de,
                        res_enrich= res_enrich,
                        annotation_obj = annotation_obj)
      
      reactive_values$dds <- dds
      reactive_values$res_de <- res_de
      reactive_values$res_enrich <- res_enrich
      reactive_values$annotation_obj <- annotation_obj
      # also creating the gtl reactive object
      reactive_values$gtl <- GeneTonic_list(
        dds = dds,
        res_de = res_de,
        res_enrich = res_enrich,
        annotation_obj = annotation_obj
      )
      
      # clean up the result object, e.g. removing the NAs in the relevant columns
      removed_genes <- is.na(res_de$log2FoldChange)
      message(
        "Removing ", sum(removed_genes),
        "/", nrow(res_de), " rows from the DE `res_de` object - log2FC values detected as NA"
      )
      reactive_values$res_de <- res_de[!removed_genes, ]
      
      reactive_values$myvst <- reactive({
        vst(reactive_values$dds)
      })
      
      reactive_values$res_enhanced <- reactive({
        get_aggrscores(
          res_enrich = reactive_values$res_enrich,
          res_de = reactive_values$res_de,
          annotation_obj = reactive_values$annotation_obj
        )
      })
      
    } else {
      message("GeneTonic info: no input data provided, upload at runtime expected")

      reactive_values$dds <- NULL
      reactive_values$res_de <- NULL
      reactive_values$res_enrich <- NULL
      reactive_values$annotation_obj <- NULL
      
      reactive_values$upload_active <- TRUE
    }
    
    
    # TODO: defining the logic of the data provided
    # 
    # if( gtl is provided) {
    #   assign gtl components to the dds, res_de, res_enrich, and annotation
    #   ready to go
    #   don't display upload gtl button
    # } else if (all components provided) {
    #   ready to go
    #   don't display upload gtl button
    # } else {
    #   display upload gtl button
    #   upon uploading, check and assign gtl components to the dds, res_de, res_enrich, and annotation
    #   ready to go
    # 
    # }

    
    # panel Welcome -----------------------------------------------------------

    output$ui_uploadgtl <- renderUI({
      validate(
        need(reactive_values$upload_active, message = "")
      )
      
      tagList(
        fluidRow(
          column(
            width = 10,
            offset = 1,
            bs4Dash::box(
              title = "Upload your data as a GeneTonicList",
              id = "box_upload",
              width = 12,
              collapsible = TRUE,
              collapsed = FALSE,
              closable = FALSE,
              fluidRow(
                column(
                  width = 6,
                  img(src = "GeneTonic/GeneTonic.png", height = "150px"),
                  fileInput(inputId = "uploadgtl", 
                            label = "Upload a GeneTonicList serialized object")
                ),
                column(
                  width = 6,
                  tags$details(
                    tags$summary("What is a GeneTonicList object?"),
                    includeMarkdown(system.file("extdata", "help_gtl.md", package = "GeneTonic")),
                    h4("Help on format, markdown"),
                    p("TODO some description and/or a link to how this format works")
                  ),
                  uiOutput("ui_describegtl")
                )
              )
            )
          )
        )
        
      )
    })
    
    output$ui_describegtl <- renderUI({
      validate(
        need(
          !is.null(reactive_values$in_gtl), message = ""
        )
      )
      
      tags$details(
        tags$summary("What's in my GeneTonicList object?"),
        p("This is the output of the `describe_gtl` function on the provided object:"),
        verbatimTextOutput("gtl_described")
      )
      
      # TODO: once not required, remove this and replace with a collapsible
      # element that contains the text output
      # nrow(reactive_values$in_gtl$res_enrich)
    })
    
    
    # TODO: maybe do it for ANY gtl, either provided as param or uploaded
    output$gtl_described <- renderText({
      validate(
        need(
          !is.null(reactive_values$in_gtl), message = "No gtl file uploaded"
        )
      )
      
      describe_gtl(reactive_values$in_gtl)
    })
    
    output$ui_panel_welcome <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl) , 
             message = "Please provide a GeneTonicList (e.g. via the 'Upload your data as a GeneTonicList' box) or its components (possible if you are calling GeneTonic from the command line).\n\n\nAll the content of this tab will be displayed upon providing this object")
      )
      tagList(
        fluidRow(
          column(
            width = 11
          ),
          column(
            width = 1,
            actionButton(
              "tour_firststeps",
              label = "", icon = icon("question-circle"),
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
            id = "card_em",
            title = "Expression Matrix",
            status = "danger",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            DT::dataTableOutput("overview_dds")
          ),
          bs4Dash::bs4Card(
            width = 6,
            id = "card_de",
            title = "DE results",
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            DT::dataTableOutput("overview_res_de")
          ),
          bs4Dash::bs4Card(
            width = 6,
            id = "card_enrich",
            title = "Functional analysis results",
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            DT::dataTableOutput("overview_res_enrich")
          ),
          bs4Dash::bs4Card(
            width = 6,
            id = "card_anno",
            title = "Annotation info",
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            DT::dataTableOutput("overview_annotation")
          )
        ),
        uiOutput("ui_infoboxes")
      )
    })
      
    output$overview_dds <- DT::renderDataTable({
      DT::datatable(
        counts(reactive_values$dds),
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })
    output$overview_res_de <- DT::renderDataTable({
      res_de <- reactive_values$res_de
      res_de <- res_de[order(res_de$padj), ]
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
          background = styleColorBar_divergent(
            as.data.frame(res_de)$log2FoldChange,
            scales::alpha("navyblue", 0.4),
            scales::alpha("darkred", 0.4)
          ),
          backgroundSize = "100% 90%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    output$overview_res_enrich <- DT::renderDataTable({
      DT::datatable(
        reactive_values$res_enrich,
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })
    output$overview_annotation <- DT::renderDataTable({
      DT::datatable(
        reactive_values$annotation_obj,
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
        value = paste0(nrow(reactive_values$dds), " genes x ", ncol(reactive_values$dds), " samples"),
        subtitle = "dds object",
        icon = icon("table"),
        color = "danger",
        width = NULL
      )
    })

    output$infobox_resde <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(
          nrow(deseqresult2df(reactive_values$res_de, FDR = input$de_fdr)),
          " DE genes"
        ),
        subtitle = "res object",
        icon = icon("vial"),
        color = "warning",
        width = NULL
      )
    })

    output$infobox_resenrich <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(
          nrow(reactive_values$res_enrich),
          " functional categories"
        ),
        subtitle = "func enrich object",
        icon = icon("share-alt"),
        color = "success",
        width = NULL
      )
    })

    output$infobox_annotation <- renderbs4ValueBox({
      bs4ValueBox(
        value = paste0(ncol(reactive_values$annotation_obj), " feature identifiers for ", nrow(reactive_values$dds), " features"),
        subtitle = "annotation object",
        icon = icon("table"),
        color = "info",
        width = NULL
      )
    })

    # panel GeneSet-Gene ------------------------------------------------------
    
    output$ui_panel_ggs <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl) , 
             message = "Gene-Geneset Panel\n\n\n\nPlease provide a GeneTonicList or its components to display the content of this tab.")
      )
      
      tagList(
        fluidRow(
          column(
            width = 11
          ),
          column(
            width = 1,
            actionButton(
              "tour_ggs",
              label = "", icon = icon("question-circle"),
              style = .helpbutton_biocstyle
            )
          )
        ),
        fluidRow(
          column(
            width = 8,
            withSpinner(
              visNetworkOutput("ggsnetwork",
                               height = "700px",
                               width = "100%"
              )
            ),
            fluidRow(
              column(
                width = 1,
                tippy::tippy(
                  actionButton(inputId = "coder_ggsnetwork",
                               label = "",
                               style = .helpbutton_biocstyle,
                               icon = icon("user-edit")),
                  "Show the code for this plot",
                  placement = "right"
                )
              )
            )
          ),
          column(
            width = 4,
            bs4Card(
              id = "box_geneset",
              title = "Geneset Box",
              width = 12,
              closable = FALSE,
              uiOutput("ui_ggs_genesetbox")
            ),
            # box(),
            hr(),
            bs4Card(
              id = "box_gene",
              title = "Gene Box",
              width = 12,
              closable = FALSE,
              uiOutput("ui_ggs_genebox")
            )
          )
        ),
        fluidRow(
          bs4Dash::bs4Card(
            width = 12,
            id = "card_ggsbackbone",
            title = "Gene-geneset graph summaries",
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            fluidRow(
              column(
                width = 8,
                uiOutput("ui_backbone_launch"),
                radioButtons(
                  inputId = "backbone_on",
                  label = "Compute the backbone on",
                  choices = c("genesets", "features"),
                  selected = "genesets",
                  inline = TRUE
                ),
                withSpinner(
                  visNetworkOutput("backbone_graph")
                )
              ),
              column(
                width = 4,
                uiOutput("ui_graph_summary")
              )
            )
          )
        )
      )
    })
    
    reactive_values$ggs_graph <- reactive({
      g <- ggs_graph(
        res_enrich = reactive_values$res_enrich,
        res_de = reactive_values$res_de,
        annotation_obj = reactive_values$annotation_obj,
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
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE
        ) %>%
        visExport(
          name = "ggs_network",
          type = "png",
          label = "Save ggs graph"
        )
    })

    output$netnode <- renderPrint({
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]

      cur_gsid <- reactive_values$res_enrich$gs_id[match(cur_sel, reactive_values$res_enrich$gs_description)]

      paste0("I'm selecting ", input$ggsnetwork_selected, ", which has index ", cur_node, " and is of type ", cur_nodetype, "this is from set", cur_gsid)
    })

    output$ui_graph_summary <- renderUI({
      tagList(
        # TODO: if other UI elements should be in, we can place them here
        h4("Highly connected genes"),
        DT::dataTableOutput("table_graph_summary")
      )
    })

    output$table_graph_summary <- DT::renderDataTable({
      g <- reactive_values$ggs_graph()

      node_degrees <- summarize_ggs_hubgenes(g)
      DT::datatable(node_degrees, escape = FALSE)
    })

    output$ui_ggs_genesetbox <- renderUI({
      tagList(
        # verbatimTextOutput("netnode"),
        plotOutput("net_sigheatplot"),
        plotOutput("sig_volcano"),
        uiOutput("ggs_geneset_info")
      )
    })

    output$net_sigheatplot <- renderPlot({
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "GeneSet",
        message = "Please select a gene set from the Gene-Geneset Graph."
      ))
      cur_gsid <- reactive_values$res_enrich$gs_id[match(input$ggsnetwork_selected, reactive_values$res_enrich$gs_description)]
      
      if (!is.null(input$exp_condition)) {
        gs_heatmap(
          reactive_values$myvst(),
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
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
          reactive_values$myvst(),
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
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

    output$sig_volcano <- renderPlot({
      g <- reactive_values$ggs_graph()
      cur_sel <- input$ggsnetwork_selected
      cur_node <- match(cur_sel, V(g)$name)
      cur_nodetype <- V(g)$nodetype[cur_node]
      validate(need(cur_nodetype == "GeneSet",
        message = ""
      ))
      cur_gsid <- reactive_values$res_enrich$gs_id[match(input$ggsnetwork_selected, reactive_values$res_enrich$gs_description)]
      
      # if (input$box_geneset$collapsed) 
      #   bs4Dash::updatebs4Card("box_geneset", action = "toggle")
      # if (!input$box_gene$collapsed) 
      #   bs4Dash::updatebs4Card("box_gene", action = "toggle")

      if (input$labels) {
        signature_volcano(
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
          geneset_id = cur_gsid,
          FDR = input$de_fdr,
          color = input$col,
          volcano_labels = Inf
        )
      } else {
        signature_volcano(
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
          geneset_id = cur_gsid,
          FDR = input$de_fdr,
          color = input$col
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
      cur_gsid <- reactive_values$res_enrich$gs_id[match(input$ggsnetwork_selected, reactive_values$res_enrich$gs_description)]

      go_2_html(cur_gsid, res_enrich)
    })

    output$ui_ggs_genebox <- renderUI({
      tagList(
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

      cur_geneid <- reactive_values$annotation_obj$gene_id[match(cur_sel, reactive_values$annotation_obj$gene_name)]
      
      # if (input$box_gene$collapsed) 
      #   bs4Dash::updatebs4Card("box_gene", action = "toggle")
      # if (!input$box_geneset$collapsed) 
      #   bs4Dash::updatebs4Card("box_geneset", action = "toggle")
      
      geneinfo_2_html(cur_sel, reactive_values$res_de)
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

      cur_geneid <- reactive_values$annotation_obj$gene_id[match(cur_sel, reactive_values$annotation_obj$gene_name)]
      gene_plot(reactive_values$dds,
        gene = cur_geneid,
        intgroup = input$exp_condition,
        annotation_obj = reactive_values$annotation_obj
      )
    })

    reactive_values$backbone_graph <- reactive({
      not_msg <- sprintf("Computing backbone on %s of the current gene-geneset graph, please hold on...", input$backbone_on)
      showNotification(not_msg)

      bbg <- ggs_backbone(
        res_enrich = reactive_values$res_enrich,
        res_de = reactive_values$res_de,
        annotation_obj = reactive_values$annotation_obj,
        n_gs = input$n_genesets,
        bb_on = input$backbone_on
        # prettify = TRUE,
        # geneset_graph_color = "gold"
      )
      return(bbg)
    })

    output$backbone_graph <- renderVisNetwork({
      # minimal example
      bbg <- reactive_values$backbone_graph()
      validate(
        need(
          {
            igraph::vcount(bbg) > 0
          },
          message = "Graph has no nodes, try increasing the number of sets to include..."
        )
      )

      visNetwork::visIgraph(bbg) %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE
        ) %>%
        visExport(
          name = "backbone_network",
          type = "png",
          label = "Save backbone graph"
        )
    })


    # panel EnrichmentMap -----------------------------------------------------
    
    output$ui_panel_em <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl) , 
             message = "Enrichment Map Panel\n\n\n\nPlease provide a GeneTonicList or its components to display the content of this tab.")
      )
      
      tagList(
        fluidRow(
          column(
            width = 11
          ),
          column(
            width = 1,
            actionButton(
              "tour_emap",
              label = "", icon = icon("question-circle"),
              style = .helpbutton_biocstyle
            )
          )
        ),
        fluidRow(
          column(
            width = 8,
            withSpinner(
              tagList(
                selectInput(
                  inputId = "emap_colorby",
                  label = "Color enrichment map by",
                  choices = colnames(reactive_values$res_enrich)[unlist(lapply(reactive_values$res_enrich, is.numeric))],
                  selected = "gs_pvalue"
                ),
                visNetworkOutput("emap_visnet",
                                 height = "700px",
                                 width = "100%"
                ),
                fluidRow(
                  column(
                    width = 1,
                    tippy::tippy(
                      actionButton(inputId = "coder_emap_visnet",
                                   label = "",
                                   style = .helpbutton_biocstyle,
                                   icon = icon("user-edit")),
                      "Show the code for this plot",
                      placement = "right"
                    )
                  )
                )
              )
            )
          ),
          column(
            width = 4,
            bs4Card(
              title = "Geneset Box",
              width = 12,
              closable = FALSE,
              uiOutput("ui_emap_sidecontent")
            )
          )
        )
        
        ,
        fluidRow(
          bs4Dash::bs4Card(
            width = 12,
            id = "card_distillery",
            title = "Geneset distillery",
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            fluidRow(
              column(
                width = 8,
                withSpinner(
                  DT::dataTableOutput("dt_distill")
                ),
                uiOutput("distill_launch"),
                numericInput(
                  inputId = "n_genesets_distill",
                  label = "Number of genesets",
                  value = min(50, nrow(reactive_values$res_enrich)), min = 1, max = nrow(reactive_values$res_enrich)
                )
              ),
              column(
                width = 4,
                bs4Card(
                  title = "Meta-geneset Box",
                  width = 12,
                  closable = FALSE,
                  uiOutput("ui_metags_sidecontent")
                )
              )
            )
          )
        )
      )
    })
        
    emap_graph <- reactive({
      emg <- enrichment_map(
        res_enrich = reactive_values$res_enrich,
        res_de = reactive_values$res_de,
        annotation_obj = reactive_values$annotation_obj,
        n_gs = input$n_genesets,
        overlap_threshold = 0.1,
        scale_edges_width = 200,
        color_by = input$emap_colorby
      )
      # rank_gs <- rank(V(emg)$name)
      # emg <- permute.vertices(emg, rank_gs)
      return(emg)
    })

    output$emap_visnet <- renderVisNetwork({
      visNetwork::visIgraph(emap_graph()) %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE
        ) %>%
        visExport(
          name = "emap_network",
          type = "png",
          label = "Save enrichment map"
        )
    })

    output$ui_emap_sidecontent <- renderUI({
      tagList(
        plotOutput("emap_sigheatplot"),
        uiOutput("emap_geneset_info")
      )
    })

    output$emap_geneset_info <- renderUI({
      cur_gsid <- reactive_values$res_enrich$gs_id[match(input$emap_visnet_selected, reactive_values$res_enrich$gs_description)]
      validate(need(!is.na(cur_gsid),
        message = ""
      ))

      # message(cur_gsid)
      # GOTERM[[cur_gsid]]
      go_2_html(cur_gsid, reactive_values$res_enrich)
    })

    output$emap_sigheatplot <- renderPlot({
      # g <- reactive_values$emap_graph()
      # cur_sel <- input$emap_visnet_selected
      # cur_node <- match(cur_sel,V(g)$name)
      # cur_nodetype <- V(g)$nodetype[cur_node]
      # validate(need(cur_nodetype == "GeneSet",
      #               message = "Please select a gene set."
      # ))
      cur_gsid <- reactive_values$res_enrich$gs_id[match(input$emap_visnet_selected, reactive_values$res_enrich$gs_description)]
      validate(need(!is.na(cur_gsid),
        message = "Please select a gene set from the Enrichment Map."
      ))

      if (!is.null(input$exp_condition)) {
        # message(cur_gsid)
        gs_heatmap(
          reactive_values$myvst(),
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
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
          reactive_values$myvst(),
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
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


    # geneset distillery
    reactive_values$distillat <- reactive({
      distillat <- distill_enrichment(
        res_enrich = reactive_values$res_enrich,
        res_de = reactive_values$res_de,
        annotation_obj = reactive_values$annotation_obj,
        n_gs = input$n_genesets_distill
      )
      return(distillat)
    })

    output$dt_distill <- DT::renderDataTable({
      dist_table <- reactive_values$distillat()$distilled_table

      # TODO: reorder the columns from the distilled table

      DT::datatable(
        dist_table[, 1:4],
        selection = "single",
        rownames = FALSE,
        options = list(
          pageLength = 50,
          scrollX = TRUE,
          scrollY = "400px"
        )
      )
    })

    output$ui_metags_sidecontent <- renderUI({
      tagList(
        plotOutput("distill_heatmap"),
        uiOutput("distill_info")
      )
    })

    output$distill_info <- renderUI({
      "Meta-geneset information"
      # TODO: structure up the content, a la gene2html
    })

    output$distill_launch <- renderUI({
      tagList(
        actionButton(
          inputId = "btn_show_emap_distilled",
          icon = icon("hubspot"),
          label = "Distill emap", style = .actionbutton_biocstyle
        )
      )
    })

    output$distill_heatmap <- renderPlot({
      dist_table <- reactive_values$distillat()$distilled_table
      s <- input$dt_distill_rows_selected

      validate(need(length(s) > 0,
        message = "Please select a meta-geneset from the table"
      ))

      selrow <- dist_table[s, ]$metags_msgs

      sel_genes <- strsplit(dist_table[s, ]$metags_genes, ",")[[1]]
      # message(length(sel_genes))
      sel_genes_id <- reactive_values$annotation_obj$gene_id[match(sel_genes, reactive_values$annotation_obj$gene_name)]
      # message(length(sel_genes_id))

      if (!is.null(input$exp_condition)) {
        gs_heatmap(
          reactive_values$myvst(),
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
          genelist = sel_genes_id,
          FDR = input$de_fdr,
          de_only = FALSE,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          center_mean = TRUE,
          scale_row = TRUE,
          anno_col_info = input$exp_condition,
          plot_title = selrow
        )
      } else {
        gs_heatmap(
          reactive_values$myvst(),
          reactive_values$res_de,
          reactive_values$res_enrich,
          annotation_obj = reactive_values$annotation_obj,
          genelist = sel_genes_id,
          FDR = input$de_fdr,
          de_only = FALSE,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          center_mean = TRUE,
          scale_row = TRUE,
          plot_title = selrow
        )
      }
    })

    # output$distill_graph <- renderPlot({
    #   plot(reactive_values$distillat()$distilled_em)
    # })

    output$distill_graph <- renderVisNetwork({
      ig <- reactive_values$distillat()$distilled_em
      # TODO: define color palette
      colpal <- colorspace::rainbow_hcl(length(unique(V(ig)$color)))[V(ig)$color]
      V(ig)$color.background <- scales::alpha(colpal, alpha = 0.8)
      V(ig)$color.highlight <- scales::alpha(colpal, alpha = 1)
      V(ig)$color.hover <- scales::alpha(colpal, alpha = 0.5)

      V(ig)$color.border <- "black"

      visNetwork::visIgraph(ig) %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE,
          selectedBy = "membership"
        )
    })



    # panel Overview ------------------------------------------------------------

    output$ui_panel_overview <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl) , 
             message = "Overview Panel\n\n\n\nPlease provide a GeneTonicList or its components to display the content of this tab.")
      )
      
      tagList(
        fluidRow(
          column(
            width = 11
          ),
          column(
            width = 1,
            actionButton(
              "tour_overview",
              label = "", icon = icon("question-circle"),
              style = .helpbutton_biocstyle
            )
          )
        ),
        fluidRow(
          bs4Dash::column(
            width = 11,
            offset = 0,
            bs4Dash::bs4Card(
              id = "tabcard_deview",
              # selected = "Geneset Volcano",
              title = "Overview",
              # side = "right",
              elevation = 2,
              width = 12,
              closable = FALSE,
              bs4Dash::tabsetPanel(
                id = "tsp1",
                type = "pills",
                selected = "Geneset Volcano",
                side = "right",
                tabPanel(
                  title = "Geneset Volcano",
                  withSpinner(
                    plotOutput("gs_volcano",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_volcano",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Geneset Volcano - simplified",
                  numericInput(
                    inputId = "gs_overlap",
                    label = "Gene Set overlap",
                    value = 0.6, min = 0, max = 1, step = 0.05,
                    width = "30%"
                  ),
                  withSpinner(
                    plotOutput("gs_volcano_simplified",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_volcano_simplified",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Enhanced Table",
                  withSpinner(
                    plotOutput("enriched_funcres",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_enhancedtable",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Enhanced Table - interactive",
                  withSpinner(
                    plotlyOutput("enriched_funcres_plotly",
                                 height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_enhancedtableinteractive",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    })
    
    output$enriched_funcres <- renderPlot({
      enhance_table(reactive_values$res_enrich, 
                    reactive_values$res_de,
                    annotation_obj = reactive_values$annotation_obj,
                    n_gs = input$n_genesets
      )
    })

    output$gs_volcano <- renderPlot({
      gs_volcano(
        get_aggrscores(reactive_values$res_enrich,
                       reactive_values$res_de,
                       annotation_obj = reactive_values$annotation_obj
        ),
        volcano_labels = input$n_genesets
      )
    })

    output$gs_volcano_simplified <- renderPlot({
      gs_volcano(
        get_aggrscores(gs_simplify(reactive_values$res_enrich, gs_overlap = input$gs_overlap),
                       reactive_values$res_de,
                       annotation_obj = reactive_values$annotation_obj
        ),
        volcano_labels = input$n_genesets
      )
    })

    output$enriched_funcres_plotly <- renderPlotly({
      ggplotly(enhance_table(reactive_values$res_enrich,
                             reactive_values$res_de,
                             annotation_obj = reactive_values$annotation_obj,
                             n_gs = input$n_genesets
      ))
    })


    # panel GSViz -----------------------------------------------------

    output$ui_panel_gsviz <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl) , 
             message = "GSViz Panel\n\n\n\nPlease provide a GeneTonicList or its components to display the content of this tab.")
      )
      
      tagList(
        fluidRow(
          column(
            width = 11
          ),
          column(
            width = 1,
            actionButton(
              "tour_gsviz",
              label = "", icon = icon("question-circle"),
              style = .helpbutton_biocstyle
            )
          )
        ),
        fluidRow(
          bs4Dash::column(
            width = 11,
            offset = 0,
            bs4Dash::bs4Card(
              id = "tabcard_genesets",
              title = "GSViz",
              # side = "right",
              elevation = 2,
              width = 12,
              closable = FALSE,
              bs4Dash::tabsetPanel(
                id = "tsp2",
                type = "pills",
                selected = "Scores Heatmap",
                side = "right",
                shiny::tabPanel(
                  title = "Scores Heatmap",
                  withSpinner(
                    plotOutput("gsscores_heatmap",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gsscores_heatmap",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Alluvial Plot",
                  withSpinner(
                    plotlyOutput("alluvial_genesets",
                                 height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_alluvial_genesets",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Summary Heatmap",
                  withSpinner(
                    plotOutput("gs_summaryheat",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_summaryheat",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Geneset MDS",
                  withSpinner(
                    plotOutput("mds_genesets",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_mds_genesets",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Summary Overview",
                  withSpinner(
                    plotOutput("gs_summaryoverview",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_summaryoverview",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Geneset Radar",
                  withSpinner(
                    plotlyOutput("gs_summaryradar",
                                 height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_summaryradar",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Geneset Dendrogram",
                  withSpinner(
                    plotOutput("gs_dendro",
                               height = "650px")
                  ),
                  fluidRow(
                    column(
                      width = 1,
                      tippy::tippy(
                        actionButton(inputId = "coder_gs_dendro",
                                     label = "",
                                     style = .helpbutton_biocstyle,
                                     icon = icon("user-edit")),
                        "Show the code for this plot",
                        placement = "right"
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    })
    
    gss_mat <- reactive({
      gs_scores(
        se = reactive_values$myvst(),
        res_de = reactive_values$res_de,
        res_enrich = reactive_values$res_enrich,
        annotation_obj = reactive_values$annotation_obj
      )
    })
    output$gsscores_heatmap <- renderPlot({
      gs_scoresheat(
        gss_mat(),
        n_gs = input$n_genesets
      )
    })

    output$alluvial_genesets <- renderPlotly({
      gs_alluvial(reactive_values$res_enrich, 
                  reactive_values$res_de, 
                  reactive_values$annotation_obj, 
                  n_gs = input$n_genesets)
    })

    output$mds_genesets <- renderPlot({
      gs_mds(reactive_values$res_enrich, 
             reactive_values$res_de, reactive_values$annotation_obj,
             mds_colorby = "z_score",
             mds_labels = input$n_genesets
      )
    })

    output$gs_summaryheat <- renderPlot({
      gs_summary_heat(reactive_values$res_enrich, 
                      reactive_values$res_de, 
                      reactive_values$annotation_obj,
        n_gs = input$n_genesets
      )
    })

    output$gs_summaryoverview <- renderPlot({
      gs_summary_overview(
        res_enrich = reactive_values$res_enhanced,
        n_gs = input$n_genesets
      )
    })

    output$gs_summaryoverview_pair <- renderPlot({
      gs_summary_overview_pair(
        res_enrich = reactive_values$res_enhanced,
        n_gs = input$n_genesets
      )
    })

    output$gs_summaryhorizon <- renderPlot({
      gs_horizon(
        res_enrich = reactive_values$res_enhanced,
        n_gs = input$n_genesets
      )
    })

    output$gs_summaryradar <- renderPlotly({
      gs_radar(
        res_enrich = reactive_values$res_enhanced,
        n_gs = input$n_genesets
      )
    })

    output$gs_dendro <- renderPlot({
      gs_dendro(
        res_enrich = reactive_values$res_enhanced,
        n_gs = input$n_genesets
      )
    })



    # panel Bookmarks ---------------------------------------------------------

    output$ui_panel_bookmarks <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl) , 
             message = "Bookmarks Panel\n\n\n\nPlease provide a GeneTonicList or its components to display the content of this tab.")
      )
      
      tagList(
        fluidRow(
          column(
            width = 11
          ),
          column(
            width = 1,
            actionButton(
              "tour_bookmarks",
              label = "", icon = icon("question-circle"),
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
            br(), br(),
            gt_downloadButton(
              "start_happyhour",
              "Start the happy hour!",
              class = "biocdlbutton",
              icon = "cocktail"
            ) # magic?
          )
        ),
        hr(),
        br(), br(), br(),
        fluidRow(
          column(
            width = 4,
            textInput(
              "se_export_name",
              label = "Choose a filename for the serialized .rds object",
              value = "se_GeneTonic_toiSEE.rds"
            )
          ),
          column(
            width = 4,
            gt_downloadButton(
              "button_iSEEexport",
              label = "Export as serialized SummarizedExperiment",
              class = "biocdlbutton",
              icon = "glasses"
            )
          )
        )
      )
    })
    
    output$ui_bookmarks <- renderUI({
      tagList(
        fluidRow(
          column(
            width = 6,
            bs4InfoBoxOutput("infobox_book_genes",
              width = 6
            ),
            h5("Bookmarked genes"),
            DT::dataTableOutput("bookmarks_genes"),
            downloadButton("btn_export_genes", label = "", class = "biocdlbutton")
            # ideally completed by a function/param to upload them
          ),
          column(
            width = 6,
            bs4InfoBoxOutput("infobox_book_genesets",
              width = 6
            ),
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
        icon = icon("bookmark"),
        iconElevation = 3,
        color = "info",
        fill = TRUE,
        width = 12
      )
    })

    output$infobox_book_genesets <- renderbs4InfoBox({
      bs4InfoBox(
        title = "Bookmarked genesets",
        value = length(reactive_values$mygenesets),
        icon = icon("bookmark"),
        iconElevation = 3,
        color = "success",
        fill = TRUE,
        width = 12
      )
    })

    output$bookmarks_genes <- DT::renderDataTable({
      book_df_genes <- reactive_values$annotation_obj[reactive_values$mygenes, ]
      datatable(book_df_genes, rownames = FALSE)
    })
    output$bookmarks_genesets <- DT::renderDataTable({
      book_df_genesets <- reactive_values$res_enrich[reactive_values$mygenesets, c("gs_id", "gs_description")]
      datatable(book_df_genesets, rownames = FALSE)
    })

    output$btn_export_genes <- downloadHandler(
      filename = function() {
        paste0("GeneTonicBookmarks_genes_", project_id, "_", gsub(" ", "_", gsub("-", "", gsub(":", "-", as.character(Sys.time())))), ".txt")
      }, content = function(file) {
        writeLines(
          text = reactive_values$mygenes,
          con = file
        )
      }
    )

    output$btn_export_genesets <- downloadHandler(
      filename = function() {
        paste0("GeneTonicBookmarks_genesets_", project_id, "_", gsub(" ", "_", gsub("-", "", gsub(":", "-", as.character(Sys.time())))), ".txt")
      }, content = function(file) {
        writeLines(
          text = reactive_values$mygenesets,
          con = file
        )
      }
    )

    output$start_happyhour <- downloadHandler(
      filename = paste0(
        Sys.Date(),
        "_", round(runif(1) * 100), # for not having all w the same name
        "_GeneTonicReport.html"
      ),
      content = function(file) {
        # temporarily switch to the temp dir, in case you do not have write permission to the current working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        # cat(tmp_content,file="GeneTonic_tempreport.Rmd",sep="\n")
        withProgress(rmarkdown::render(
          input = system.file("extdata", "cocktail_recipe.Rmd", package = "GeneTonic"),
          output_file = file,
          # fragment.only = TRUE,
          quiet = TRUE
        ),
        message = "Generating the html report",
        detail = "This can take some time"
        )
      }
    )

    output$button_iSEEexport <- downloadHandler(
      filename = function() {
        input$se_export_name
      }, content = function(file) {
        se <- export_for_iSEE(reactive_values$dds, reactive_values$res_de)
        saveRDS(se, file = file)
      }
    )
    
    output$ui_gtl_download <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl), "provide GTL TODO")
      )
      gt_downloadButton(
        "btn_download_gtl",
        "Save as gtl serialized object TODO",
        class = "biocdlbutton",
        icon = "gift"
      )
    })
    
    output$btn_download_gtl <- downloadHandler(
      filename = function() {
        "exported_gtl.RDS"
      }, content = function(file) {
        showNotification("Saving... TODO", type = "default")
        saveRDS(reactive_values$gtl, file = file)
        showNotification("Done! TODO", type = "message")
      }
    )

    output$sessioninfo <- renderPrint({
      sessionInfo()
    })

    # controlbar --------------------------------------------------------------
    
    output$ui_controlbar <- renderUI({
      validate(
        need(!is.null(reactive_values$gtl) , 
             message = "\n\n\nPlease provide a GeneTonicList or its components to display the content of the control bar.")
      )
      
      # message(nrow(reactive_values$res_enrich))
      
      tagList(
        numericInput(
          inputId = "de_fdr",
          label = "False Discovery Rate (FDR) for DE",
          value = 0.05, min = 0.0001, max = 1, step = 0.01
        ),
        numericInput(
          inputId = "n_genesets",
          label = "Number of genesets",
          value = 15, min = 1, max = nrow(reactive_values$res_enrich)
        ),
        selectInput("exp_condition",
                    label = "Group/color by: ",
                    choices = c(NULL, names(colData(reactive_values$dds))), selected = NULL, multiple = TRUE
        ),
        colourInput("col", "Select colour for volcano plot", "#1a81c2",
                    returnName = TRUE,
                    allowTransparent = TRUE
        ),
        checkboxInput("labels", label = "Display all labels", value = FALSE)
      )
    })
    
    outputOptions(output, "ui_controlbar", suspendWhenHidden = FALSE)

    # observers ---------------------------------------------------------------

    observeEvent(input$tour_firststeps, {
      tour <- read.delim(system.file("extdata", "tour_welcome.txt", package = "GeneTonic"),
        sep = ";", stringsAsFactors = FALSE,
        row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_ggs, {
      tour <- read.delim(system.file("extdata", "tour_ggs.txt", package = "GeneTonic"),
        sep = ";", stringsAsFactors = FALSE,
        row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_emap, {
      tour <- read.delim(system.file("extdata", "tour_emap.txt", package = "GeneTonic"),
        sep = ";", stringsAsFactors = FALSE,
        row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_overview, {
      tour <- read.delim(system.file("extdata", "tour_overview.txt", package = "GeneTonic"),
        sep = ";", stringsAsFactors = FALSE,
        row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_gsviz, {
      tour <- read.delim(system.file("extdata", "tour_gsviz.txt", package = "GeneTonic"),
        sep = ";", stringsAsFactors = FALSE,
        row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_bookmarks, {
      tour <- read.delim(system.file("extdata", "tour_bookmarks.txt", package = "GeneTonic"),
        sep = ";", stringsAsFactors = FALSE,
        row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })


    # observe({
    #   print(input$gt_tabs)
    # })

    observeEvent(input$btn_first_help, {
      showModal(
        modalDialog(
          title = "First Help Info", size = "l", fade = TRUE,
          footer = NULL, easyClose = TRUE,
          tagList(
            includeMarkdown(system.file("extdata", "GeneTonic101.md", package = "GeneTonic")),
          )
        )
      )
    })

    observeEvent(input$btn_docs_vignette, {
      path <- system.file("doc", "GeneTonic_manual.html", package = "GeneTonic")
      if (path == "") {
        showNotification("This vignette has not been built on this system - Opening the online documentation. Please note that the versions might not be coincident!", type = "warning")
      } else {
        browseURL(path)
      }
    })


    observeEvent(input$btn_info_session, {
      showModal(
        modalDialog(
          title = "Session information", size = "l", fade = TRUE,
          footer = NULL, easyClose = TRUE,
          tagList(
            tags$code("> sessionInfo()"),
            renderPrint({
              sessionInfo()
            })
          )
        )
      )
    })

    observeEvent(input$btn_info_genetonic, {
      showModal(
        modalDialog(
          title = "About GeneTonic", size = "l", fade = TRUE,
          footer = NULL, easyClose = TRUE,
          tagList(
            includeMarkdown(system.file("extdata", "about.md", package = "GeneTonic")),
            renderPrint({
              citation("GeneTonic")
            })
          )
        )
      )
    })

    observeEvent(input$btn_show_emap_distilled, {
      showModal(
        modalDialog(
          title = "distillery", size = "l", fade = TRUE,
          footer = NULL, easyClose = TRUE,
          visNetworkOutput("distill_graph")
        )
      )
    })
    
    observeEvent(input$uploadgtl, {
      reactive_values$in_gtl <- readRDS(input$uploadgtl$datapath)
      
      if (is.null(checkup_gtl(reactive_values$in_gtl))) {
        reactive_values$gtl <- reactive_values$in_gtl
        reactive_values$dds <- reactive_values$in_gtl$dds
        reactive_values$res_de <- reactive_values$in_gtl$res_de
        reactive_values$res_enrich <- reactive_values$in_gtl$res_enrich
        reactive_values$annotation_obj <- reactive_values$in_gtl$annotation_obj
        
        showNotification(
          ui = "Upload complete! Re-open the 'File Upload' collapsible box if you want to upload another gtl object.",
          type = "message"
        )
        
        
        removed_genes <- is.na(reactive_values$res_de$log2FoldChange)
        message(
          "Removing ", sum(removed_genes),
          "/", nrow(reactive_values$res_de), " rows from the DE `res_de` object - log2FC values detected as NA"
        )
        reactive_values$res_de <- reactive_values$res_de[!removed_genes, ]
        
        
        reactive_values$myvst <- reactive({
          vst(reactive_values$dds)
        })
        
        reactive_values$res_enhanced <- reactive({
          get_aggrscores(
            res_enrich = reactive_values$res_enrich,
            res_de = reactive_values$res_de,
            annotation_obj = reactive_values$annotation_obj
          )
        })
        
        bs4Dash::updateBox("box_upload", action = "toggle")
      } else {
        showNotification("wrong file format, TODO!")
      }
      
    })


    # bookmarker --------------------------------------------------------------
    observeEvent(input$bookmarker, {
      if (input$gt_tabs == "tab_welcome") {
        showNotification("Welcome to GeneTonic!")
      } else if (input$gt_tabs == "tab_ggs") {
        g <- reactive_values$ggs_graph()
        cur_sel <- input$ggsnetwork_selected
        if (cur_sel == "") {
          showNotification("Select a node in the network to bookmark it", type = "warning")
        } else {
          cur_node <- match(cur_sel, V(g)$name)
          cur_nodetype <- V(g)$nodetype[cur_node]

          if (cur_nodetype == "Feature") {
            cur_sel_id <- reactive_values$annotation_obj$gene_id[match(cur_sel, reactive_values$annotation_obj$gene_name)]
            if (cur_sel_id %in% reactive_values$mygenes) {
              showNotification(sprintf("The selected gene %s (%s) is already in the set of the bookmarked genes.", cur_sel, cur_sel_id), type = "default")
            } else {
              reactive_values$mygenes <- unique(c(reactive_values$mygenes, cur_sel_id))
              # message("there go your genes... ", reactive_values$mygenes)
              showNotification(sprintf("Added %s (%s) to the bookmarked genes. The list contains now %d elements", cur_sel, cur_sel_id, length(reactive_values$mygenes)), type = "message")
            }
          } else if (cur_nodetype == "GeneSet") {
            cur_sel_id <- reactive_values$res_enrich$gs_id[match(cur_sel, reactive_values$res_enrich$gs_description)]
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
        g <- reactive_values$ggs_graph()
        cur_sel <- input$emap_visnet_selected
        cur_sel_id <- reactive_values$res_enrich$gs_id[match(cur_sel, reactive_values$res_enrich$gs_description)]

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
      else if (input$gt_tabs == "tab_bookmarks") {
        showNotification("You are already in the Bookmarks tab...")
      }
    })


    # code popups -------------------------------------------------------------
    observeEvent(input$coder_gs_volcano, {
      mycode <- c(
        .gt_code_setup,
        "gs_volcano(",
        "  get_aggrscores(res_enrich = res_enrich,",
        "                 res_de = res_de,",
        "                 annotation_obj = annotation_obj",
        "  ),",
        paste0("  volcano_labels = ", input$n_genesets),
        ")",
        .gt_code_closeup
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_volcano",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
            # renderPrint({writeLines(mycode)})
          )
        )
      )
    })

    observeEvent(input$coder_gs_volcano_simplified, {
      mycode <- c(
        .gt_code_setup,

        "gs_volcano(",
        "  get_aggrscores(",
        "    res_enrich = gs_simplify(res_enrich = res_enrich, ",
        paste0("                             gs_overlap = ", input$gs_overlap, "),"),
        "    res_de = res_de,",
        "    annotation_obj = annotation_obj",
        "  ),",
        paste0("  volcano_labels = ", input$n_genesets),
        ")",

        .gt_code_closeup
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_volcano_simplified",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_gs_enhancedtable, {
      mycode <- c(
        .gt_code_setup,
        "enhance_table(res_enrich = res_enrich,",
        "              res_de = res_de,",
        "              annotation_obj = annotation_obj,",
        paste0("              n_gs = ", input$n_genesets, ")"),
        .gt_code_closeup
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_enhancedtable",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_gs_enhancedtableinteractive, {
      mycode <- c(
        .gt_code_setup,
        "library('plotly')",
        "ggplotly(",
        "  enhance_table(res_enrich = res_enrich,",
        "                res_de = res_de,",
        "                annotation_obj = annotation_obj,",
        paste0("                n_gs = ", input$n_genesets, ")"),
        ")"

      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_enhancedtableinteractive",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_gsscores_heatmap, {
      mycode <- c(
        .gt_code_setup,

        "myvst <- vst(dds)",
        "gss_mat <- gs_scores(",
        "  se = myvst,",
        "  res_de = res_de,",
        "  res_enrich = res_enrich,",
        "  annotation_obj = annotation_obj",
        ")",
        "",
        "gs_scoresheat(",
        "  gss_mat,",
        paste0("  n_gs = ", input$n_genesets),
        ")",

        .gt_code_closeup
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gsscores_heatmap",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_alluvial_genesets, {
      mycode <- c(
        .gt_code_setup,

        "gs_alluvial(res_enrich = res_enrich,",
        "            res_de = res_de, ",
        "            annotation_obj = annotation_obj, ",
        paste0("            n_gs = ", input$n_genesets, ")")

      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-alluvial_genesets",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_gs_summaryheat, {
      mycode <- c(
        .gt_code_setup,

        "gs_summary_heat(",
        "  res_enrich = res_enrich, ",
        "  res_de = res_de, ",
        "  annotation_obj = annotation_obj,",
        paste0("  n_gs = ", input$n_genesets),
        ")",

        .gt_code_closeup
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_summaryheat",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_mds_genesets, {
      mycode <- c(
        .gt_code_setup,

        "gs_mds(res_enrich = res_enrich, ",
        "       res_de = res_de, ",
        "       annotation_obj = annotation_obj,",
        "       mds_colorby = 'z_score',",
        paste0("       mds_labels = ", input$n_genesets),
        ")",

        .gt_code_closeup
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-mds_genesets",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_gs_summaryoverview, {
      mycode <- c(
        .gt_code_setup,

        "res_enhanced <- get_aggrscores(",
        "  res_enrich = res_enrich,",
        "  res_de = res_de,",
        "  annotation_obj = annotation_obj",
        ")",
        "",
        "gs_summary_overview(",
        "  res_enrich = res_enhanced,",
        paste0("  n_gs = ", input$n_genesets),
        ")",

        .gt_code_closeup
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_summaryoverview",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_gs_summaryradar, {
      mycode <- c(
        .gt_code_setup,

        "res_enhanced <- get_aggrscores(",
        "  res_enrich = res_enrich,",
        "  res_de = res_de,",
        "  annotation_obj = annotation_obj",
        ")",
        "",
        "gs_radar(",
        "  res_enrich = res_enhanced,",
        paste0("  n_gs = ", input$n_genesets),
        ")"

      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_summaryradar",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_gs_dendro, {
      mycode <- c(
        .gt_code_setup,

        "res_enhanced <- get_aggrscores(",
        "  res_enrich = res_enrich,",
        "  res_de = res_de,",
        "  annotation_obj = annotation_obj",
        ")",
        "",
        "gs_dendro(",
        "  res_enrich = res_enhanced,",
        paste0("  n_gs = ", input$n_genesets),
        ")"
      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_dendro",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_ggsnetwork, {
      mycode <- c(
        .gt_code_setup,

        "g <- ggs_graph(",
        "  res_enrich = res_enrich,",
        "  res_de = res_de,",
        "  annotation_obj = annotation_obj,",
        paste0("  n_gs = ", input$n_genesets, ","),
        "  prettify = TRUE,",
        "  geneset_graph_color = 'gold'",
        ")",
        "",
        "library('visNetwork')",
        "library('magrittr')",
        "visNetwork::visIgraph(g) %>%",
        "  visOptions(",
        "    highlightNearest = list(",
        "      enabled = TRUE,",
        "      degree = 1,",
        "      hover = TRUE",
        "    ),",
        "    nodesIdSelection = TRUE",
        "  ) %>%",
        "  visExport(",
        "    name = 'ggs_network',",
        "    type = 'png',",
        "    label = 'Save ggs graph'",
        "  )"

      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-gs_ggsnetwork",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })

    observeEvent(input$coder_emap_visnet, {
      mycode <- c(
        .gt_code_setup,

        "emg <- enrichment_map(",
        "  res_enrich = res_enrich,",
        "  res_de = res_de,",
        "  annotation_obj = annotation_obj,",
        paste0("  n_gs = ", input$n_genesets, ","),
        "  overlap_threshold = 0.1,",
        "  scale_edges_width = 200,",
        paste0("  color_by = '", input$emap_colorby, "'"),
        ")",
        "",
        "library('visNetwork')",
        "library('magrittr')",
        "visNetwork::visIgraph(emg) %>%",
        "  visOptions(",
        "    highlightNearest = list(",
        "      enabled = TRUE,",
        "      degree = 1,",
        "      hover = TRUE",
        "    ),",
        "    nodesIdSelection = TRUE",
        "  ) %>%",
        "  visExport(",
        "    name = 'emap_network',",
        "    type = 'png',",
        "    label = 'Save enrichment map'",
        "  )"

      )

      showModal(
        modalDialog(
          title = "About this plot", size = "l", fade = TRUE, footer = NULL, easyClose = TRUE,
          tagList(
            shinyAce::aceEditor(
              "editor-emap_visnet",
              theme = "solarized_light",
              readOnly = TRUE,
              value = paste(mycode, collapse = "\n")
            )
          )
        )
      )
    })



    # observeEvent(input$start_happyhour, {
    # showNotification("The happy hour is on! Please wait for the report to be fully compiled",
    # type = "message")
    # })
  }
  # nocov end

  shinyApp(ui = genetonic_ui, server = genetonic_server)
}
