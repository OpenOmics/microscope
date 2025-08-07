## app.R
##
## This Shiny application implements an interactive dashboard for exploring
## diversity and taxonomy data.  The UI structure is based on the layout
## patterns found in the `bioNPS` project (see https://github.com/abenedetti/bioNPS).
## Specifically, it uses shinydashboard to provide a persistent sidebar for
## navigating between sections, while each tab presents its own set of
## controls arranged horizontally using `fluidRow()` and `column()` calls.
## Download buttons and input controls are grouped and spaced using
## `div(style=...)` wrappers to ensure consistent margins.

library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(plotly)
library(bslib)
library(tibble)
library(ape)
library(vegan)
library(usedist)
library(scales)

## Load helper functions from utils.R.  These functions define:
## * load_diversity_metrics() - loads alpha/beta diversity from files
## * tax_collapse()           - collapses count tables to a taxonomic level
## * plot_alpha_diversity()   - returns a ggplot for alpha diversity
## * plot_beta_pcoa()         - returns a ggplot for beta diversity PCoA
## * plot_tax_stacked_bar()   - returns a ggplot for stacked bar plots
source("../R/utils.R")

## Define a colour palette for discrete groups.  The palette is reversed
## relative to the original specification so that more common taxa appear
## earlier in the legend.  You may adjust or replace this palette as
## desired.
col21 <- rev(c(
  "tomato1", "darkblue", "turquoise1", "lightblue", "darkred",
  "mediumblue", "purple", "bisque", "greenyellow", "yellow",
  "violetred2", "darkgreen", "darkgoldenrod1", "deeppink3",
  "cadetblue4", "orchid2", "seagreen3", "purple4",
  "dodgerblue2", "red", "gray27"
))

# User interface ------------------
ui <- dashboardPage(
  title = "Microscope",
  dashboardHeader(
    title = tags$div("Microscope",
                     style = "font-size:32px; font-style:italic; font-weight:bold; text-align:left;")
  ),
  dashboardSidebar(width = 250,
    tags$div(
         style = "text-align: center; padding: 50px 0px 50px 0px;",
         tags$img(src = "logo.microscope.nobackground.png", height = "220px")
    ),
    sidebarMenu(
      menuItem("Home", tabName = "home_tab", icon = icon("home")),
      menuItem("Metadata", tabName = "metadata", icon = icon("table")),
      menuItem("Alpha Diversity", tabName = "alpha", icon = icon("chart-bar")),
      menuItem("Beta Diversity", tabName = "beta", icon = icon("chart-line")),
      menuItem("Taxonomy", tabName = "taxo", icon = icon("layer-group"))
    ),
    div(
      style = "position: absolute; bottom: 20px; width: 100%; display: flex; justify-content: center; gap: 25px;",
      tags$img(src = "logo.rtb.png", height = "60px", style = "margin-bottom:10px;"),
      tags$img(src = "logo.openomics.png", height = "60px")
    )
  ),
  dashboardBody(
    tabItems(
      ## Home tab --------------------------------------------------
      tabItem(tabName = "home_tab",
              fluidRow(
                column(width = 8, offset = 1,
                       tags$hr(),
                       # First line: large, italic, centered
                       tags$p("Welcome to Microscope",
                              style = "font-size:40px; font-style:italic; text-align:left;"),
                       tags$hr(),
                       tags$hr(),
                       
                       # All other text: larger font, left-aligned
                       tags$p("This dashboard allows you to explore microbiome diversity and taxonomy data.",
                              style = "font-size:25px; text-align:left;"),
                       
                       tags$p("Use the sidebar to navigate to Alpha Diversity, Beta Diversity, or Taxonomy.",
                              style = "font-size:25px; text-align:left;"),
                       
                       tags$hr(),
                       
                       # Final message: larger, centered
                       tags$p("Upload your metadata and diversity metrics to get started!",
                              style = "font-size:25px; text-align:left;"),
                  
                )
              )
      ),
      ## Metadata tab ------------------------------------------------------
      tabItem(tabName = "metadata",
        fluidRow(
          column(6,
            div(style = "margin-bottom:15px;",
                fileInput("metadata", "Upload Metadata CSV", accept = ".csv")
            )
          ),
        ),
        fluidRow(
          column(12,
            card(card_header("Metadata Table"), tableOutput("metadata_preview"))
          )
        )
      ),
      ## Alpha diversity tab ----------------------------------------------
      tabItem(tabName = "alpha",
        fluidRow(
          column(4,
            div(style = "margin-bottom:15px;",
                fileInput("diversity_zip", "Upload Diversity Files (.zip)", accept = ".zip")
            )
          ),
          column(4,
            uiOutput("group_column_ui")
          )
        ),
        fluidRow(
          column(12,
            card(card_header("Alpha Diversity Plot"), plotlyOutput("alpha_plot", height = "600px"))
          )
        ),
        fluidRow(
          uiOutput('download_alpha_plot_ui'),
        )
      ),
      ## Beta diversity tab -----------------------------------------------
      tabItem(tabName = "beta",
        fluidRow(
          column(4, uiOutput("beta_matrix_ui")),
          column(4, uiOutput("beta_group_ui"))
        ),
        fluidRow(
          column(12,
            card(card_header("PCoA Plot"), plotOutput("beta_plot", height = "600px"))
          )
        ),
        fluidRow(
          uiOutput('download_beta_plot_ui')
        )
      ),
      ## Taxonomy bar plot tab -------------------------------------------
      tabItem(tabName = "taxo",
        fluidRow(
          column(3,
            div(style = "margin-bottom:15px;",
                fileInput("count_matrix", "Upload Count Matrix (CSV/TSV)", accept = c(".csv", ".tsv"))
            )
          ),
          column(3, uiOutput("tax_level_ui")),
          column(3, uiOutput("split_by_1_ui")),
          column(3, uiOutput("split_by_2_ui"))
        ),
        fluidRow(
          column(3,
            div(style = "margin-bottom:15px;",
                sliderInput("top_n", "Number of Top Taxa:", min = 1, max = 20, value = 9)
            )
          ),
          column(3,
            checkboxInput("show_names", "Show Sample Names", FALSE)
          )
        ),
        fluidRow(
          column(12,
            card(card_header("Stacked Bar Plot"), plotlyOutput("tax_bar_plot", height = "800px"))
          )
        ),
        fluidRow(
          column(3,
                 uiOutput('download_tax_plot_ui'),
                 uiOutput('download_tax_data_ui')
          )
        )
      )
    )
  )
)

# Server ------------------
server <- function(input, output, session) {
  ## Reactive values to store loaded data ----------
  alpha_data <- reactiveVal(NULL)
  beta_data  <- reactiveVal(NULL)
  pheno_data <- reactiveVal(NULL)
  count_matrix_data <- reactiveVal(NULL)
  tax_data   <- reactiveVal(NULL)

  ## Handle metadata upload ----------
  observeEvent(input$metadata, {
    req(input$metadata)
    meta <- readr::read_csv(input$metadata$datapath)
    pheno_data(meta)
    output$metadata_preview <- renderTable({ meta })
  })

  ## Automatically load diversity metrics when the folder changes ----------
  observeEvent(input$diversity_zip, {
    req(input$metadata, input$diversity_zip)
    # Temporary directory for extracted files
    extract_dir <- file.path(tempdir(), "diversity_data")
    unlink(extract_dir, recursive = TRUE) # clean old extractions
    dir.create(extract_dir, showWarnings = FALSE)
    # Extract uploaded ZIP file
    unzip(input$diversity_zip$datapath, exdir = extract_dir,junkpaths = T)
    # Load diversity metrics
    diversity <- load_diversity_metrics(input$metadata$datapath, extract_dir)
    alpha_data(diversity$alpha)
    beta_data(diversity)
  })

  ## UI for alpha diversity group selection ----------
  output$group_column_ui <- renderUI({
    req(alpha_data())
    selectInput("group_by", "Group by Column:",
      choices = setdiff(colnames(alpha_data()), c("sampleID", "Metric", "Value")))
  })

  ## Render alpha diversity plot ----------
  output$alpha_plot <- renderPlotly({
    req(alpha_data(), input$group_by)
    p = plot_alpha_diversity(alpha_data(), input$group_by)
    plotly::ggplotly(p)
  })

  ## UI for beta diversity matrix selection ----------
  output$beta_matrix_ui <- renderUI({
    req(beta_data())
    selectInput("beta_matrix", "Distance Matrix:", choices = names(beta_data()$beta))
  })

  ## UI for beta diversity grouping variable ----------
  output$beta_group_ui <- renderUI({
    req(beta_data())
    selectInput("beta_group", "Grouping Variable:",
      choices = setdiff(colnames(beta_data()$metadata), "sampleID"))
  })

  ## Render beta diversity plot ----------
  output$beta_plot <- renderPlot({
    req(beta_data(), input$beta_matrix, input$beta_group)
    dm <- beta_data()$beta[[input$beta_matrix]]
    metadata <- beta_data()$metadata
    plot_beta_pcoa(dm = dm, metatable = metadata,
                   grp = input$beta_group, dm_name = input$beta_matrix,
                   permtest = TRUE, grp_continuous = FALSE)
  })

  ## Handle count matrix upload for taxonomy bar plot ----------
  observeEvent(input$count_matrix, {
    req(input$count_matrix)
    file_ext <- tools::file_ext(input$count_matrix$name)
    cm <- if (file_ext == "tsv") {
      readr::read_tsv(input$count_matrix$datapath)
    } else {
      readr::read_csv(input$count_matrix$datapath)
    }
    count_matrix_data(cm)
    ## Render taxonomic level selection once the matrix is loaded 
    output$tax_level_ui <- renderUI({
      selectInput("tax_level", "Taxonomic Level:",
        choices = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
        selected = "genus")
    })
  })

  ## UI for split_by selectors ----------
  output$split_by_1_ui <- renderUI({
    req(pheno_data())
    selectInput("split_by_1", "Facet by:", choices = c("None", colnames(pheno_data())), selected = "None")
  })
  output$split_by_2_ui <- renderUI({
    req(pheno_data())
    selectInput("split_by_2", "Second Facet:", choices = c("None", colnames(pheno_data())), selected = "None")
  })

  ## Render taxonomy bar plot ----------
  output$tax_bar_plot <- renderPlotly({
    req(count_matrix_data(), pheno_data(), input$tax_level)
    tax_result <- tax_collapse(count_matrix_data(), tax = input$tax_level, top_n = input$top_n)
    split_1 <- if (input$split_by_1 == "None") NULL else input$split_by_1
    split_2 <- if (input$split_by_2 == "None") NULL else input$split_by_2
    p = plot_tax_stacked_bar(exps = tax_result[[3]], pheno = pheno_data(),
                             split_by_1 = split_1, split_by_2 = split_2,
                             show_names = input$show_names)
    plotly::ggplotly(p)
  })
  
  ## Show download buttons only if files are available ---------------
  output$download_alpha_plot_ui <- renderUI({
    req(alpha_data(), input$group_by) # Show button only when data exists
    downloadButton("download_alpha_plot", "Download Plot")
  })
  output$download_beta_plot_ui <- renderUI({
    req(beta_data(), input$group_by) # Show button only when data exists
    downloadButton("download_beta_plot", "Download Plot")
  })
  output$download_tax_plot_ui <- renderUI({
    req(count_matrix_data(), pheno_data()) # Show button only when data exists
    downloadButton("download_tax_plot", "Download Plot")
  })
  output$download_tax_data_ui <- renderUI({
    req(count_matrix_data(), pheno_data(),input$tax_level, input$top_n) # Show button only when data exists
    downloadButton("download_tax_data", "Download Data")
  })

  ## Download handlers -------------
  output$download_alpha_plot <- downloadHandler(
    filename = function() { "alpha_diversity_plot.png" },
    content = function(file) {
      req(alpha_data(), input$group_by)
      p <- plot_alpha_diversity(alpha_data(), input$group_by)
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  output$download_beta_plot <- downloadHandler(
    filename = function() { "beta_diversity_plot.png" },
    content = function(file) {
      req(beta_data(), input$beta_matrix, input$beta_group)
      dm <- beta_data()$beta[[input$beta_matrix]]
      metadata <- beta_data()$metadata
      p <- plot_beta_pcoa(dm, metadata, input$beta_group, input$beta_matrix)
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  output$download_tax_plot <- downloadHandler(
    filename = function() { 
      req(count_matrix_data(), pheno_data(),input$tax_level, input$top_n)
      paste0("taxonomydata_", input$tax_level, "_top", input$top_n, ".png")
    },
    content = function(file) {
      req(count_matrix_data(), pheno_data())
      tax_result <- tax_collapse(count_matrix_data(), tax = input$tax_level, top_n = input$top_n)
      split_1 <- if (input$split_by_1 == "None") NULL else input$split_by_1
      split_2 <- if (input$split_by_2 == "None") NULL else input$split_by_2
      p <- plot_tax_stacked_bar(exps = tax_result[[3]], pheno = pheno_data(),
                                split_by_1 = split_1, split_by_2 = split_2,
                                show_names = input$show_names)
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )
  output$download_tax_data <- downloadHandler(
    filename = function() { 
      req(count_matrix_data(), pheno_data(),input$tax_level, input$top_n)
      paste0("taxonomydata_", input$tax_level, "_top", input$top_n, ".csv")
    },
    content = function(file) {
      req(count_matrix_data(), pheno_data())
      tax_result <- tax_collapse(count_matrix_data(), tax = input$tax_level, top_n = input$top_n)[[2]]
      write_csv(tax_result, file)
    }
  )
}

## Run the application
shinyApp(ui, server)