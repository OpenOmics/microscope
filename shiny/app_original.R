# app.R

packages_to_load = c(
  "shiny", "ggplot2", "ggpubr", "dplyr", "tidyr", "readr", "stringr",
  "plotly", "bslib", "tibble", "ape", "vegan", "usedist", "scales", "shinyFiles"
)

# Load packages quietly
invisible(suppressPackageStartupMessages({
  lapply(packages_to_load, library, character.only = TRUE)
}))

# Load helper functions
source("utils.R")

# Color palette for discrete groups
col21 <- rev(c(
  "tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
  "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
  "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"
))

ui <- fluidPage(
  titlePanel("Microscope"),
  theme = bs_theme(version = 5),
  
  tabsetPanel(
    # Metadata Preview Tab
    tabPanel("Metadata Preview",
             sidebarLayout(
               sidebarPanel(
                 div(style="margin-bottom:15px;",
                     fileInput("metadata", "Upload Metadata CSV", accept = ".csv")
                 ),
                 div(style="margin-bottom:10px;",
                     downloadButton("download_metadata", "Download Metadata CSV")
                 )
               ),
               mainPanel(
                 card(card_header("Metadata Table"), tableOutput("metadata_preview"))
               )
             )
    ),
    
    # Alpha Diversity Tab
    tabPanel("Alpha Diversity",
             sidebarLayout(
               sidebarPanel(
                 div(style="margin-bottom:15px;",
                     shinyDirButton("div_dir", "Select Diversity Folder", "Select a folder")
                 ),
                 uiOutput("group_column_ui"),
                 br(),
                 div(style="margin-bottom:10px;", downloadButton("download_alpha_plot", "Download Alpha Plot")),
                 div(style="margin-bottom:10px;", downloadButton("download_alpha_data", "Download Alpha Data"))
               ),
               mainPanel(
                 card(card_header("Alpha Diversity Plot"), plotlyOutput("alpha_plot", height = "600px"))
               )
             )
    ),
    
    # Beta Diversity Tab
    tabPanel("Beta Diversity",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("beta_matrix_ui"),
                 uiOutput("beta_group_ui"),
                 br(),
                 div(style="margin-bottom:10px;", downloadButton("download_beta_plot", "Download Beta Plot")),
                 div(style="margin-bottom:10px;", downloadButton("download_beta_data", "Download PCoA Data"))
               ),
               mainPanel(
                 card(card_header("PCoA Plot"), plotOutput("beta_plot", height = "600px"))
               )
             )
    ),
    
    # Taxonomy Bar Plot Tab
    tabPanel("Taxonomy Bar Plot",
             sidebarLayout(
               sidebarPanel(
                 div(style="margin-bottom:15px;",
                     fileInput("count_matrix", "Upload Count Matrix (CSV/TSV)", accept = c(".csv", ".tsv"))
                 ),
                 uiOutput("tax_level_ui"),
                 uiOutput("split_by_1_ui"),
                 uiOutput("split_by_2_ui"),
                 div(style="margin-bottom:15px;",
                     sliderInput("top_n", "Number of Top Taxa:", min = 5, max = 30, value = 9)
                 ),
                 checkboxInput("show_names", "Show Sample Names", FALSE),
                 br(),
                 div(style="margin-bottom:10px;", downloadButton("download_tax_plot", "Download Taxonomy Plot")),
                 div(style="margin-bottom:10px;", downloadButton("download_tax_data", "Download Taxonomy Data"))
               ),
               mainPanel(
                 card(card_header("Stacked Bar Plot"), plotOutput("tax_bar_plot", height = "800px"))
               )
             )
    )
  )
)

server <- function(input, output, session) {
  alpha_data <- reactiveVal(NULL)
  beta_data <- reactiveVal(NULL)
  pheno_data <- reactiveVal(NULL)
  count_matrix_data <- reactiveVal(NULL)
  tax_data <- reactiveVal(NULL)
  
  volumes <- c(Home = fs::path_home(), "Root" = "/")
  shinyDirChoose(input, "div_dir", roots = volumes, session = session)
  
  # ---- Metadata ----
  observeEvent(input$metadata, {
    meta <- read_csv(input$metadata$datapath)
    pheno_data(meta)
    output$metadata_preview <- renderTable({ head(meta, 10) })
  })
  
  # ---- Auto-load diversity when folder selected ----
  observeEvent(input$div_dir, {
    req(input$metadata)
    dir_path <- parseDirPath(volumes, input$div_dir)
    req(dir_path)
    diversity <- load_diversity_metrics(input$metadata$datapath, dir_path)
    alpha_data(diversity$alpha)
    beta_data(diversity)
  })
  
  # ---- Alpha Diversity ----
  output$group_column_ui <- renderUI({
    req(alpha_data())
    selectInput("group_by", "Group by Column:",
                choices = setdiff(colnames(alpha_data()), c("sampleID", "Metric", "Value")))
  })
  
  output$alpha_plot <- renderPlotly({
    req(alpha_data(), input$group_by)
    plot_alpha_diversity(alpha_data(), input$group_by)
  })
  
  # ---- Beta Diversity ----
  output$beta_matrix_ui <- renderUI({
    req(beta_data())
    selectInput("beta_matrix", "Distance Matrix:",
                choices = names(beta_data()$beta))
  })
  
  output$beta_group_ui <- renderUI({
    req(beta_data())
    selectInput("beta_group", "Grouping Variable:",
                choices = setdiff(colnames(beta_data()$metadata), "sampleID"))
  })
  
  output$beta_plot <- renderPlot({
    req(beta_data(), input$beta_matrix, input$beta_group)
    dm <- beta_data()$beta[[input$beta_matrix]]
    metadata <- beta_data()$metadata
    plot_beta_pcoa(dm = dm, metatable = metadata,
                   grp = input$beta_group, dm_name = input$beta_matrix,
                   permtest = TRUE, grp_continuous = FALSE)
  })
  
  # ---- Taxonomy Bar Plot ----
  observeEvent(input$count_matrix, {
    req(input$count_matrix)
    file_ext <- tools::file_ext(input$count_matrix$name)
    cm <- if (file_ext == "tsv") read_tsv(input$count_matrix$datapath) else read_csv(input$count_matrix$datapath)
    count_matrix_data(cm)
    
    output$tax_level_ui <- renderUI({
      selectInput("tax_level", "Taxonomic Level:",
                  choices = c("kingdom","phylum","class","order","family","genus","species"),
                  selected = "genus")
    })
  })
  
  output$split_by_1_ui <- renderUI({
    req(pheno_data())
    selectInput("split_by_1", "Facet by:", choices = c("None", colnames(pheno_data())), selected = "None")
  })
  
  output$split_by_2_ui <- renderUI({
    req(pheno_data())
    selectInput("split_by_2", "Second Facet:", choices = c("None", colnames(pheno_data())), selected = "None")
  })
  
  output$tax_bar_plot <- renderPlot({
    req(count_matrix_data(), pheno_data(), input$tax_level)
    tax_result <- tax_collapse(count_matrix_data(), tax = input$tax_level, top_n = input$top_n)
    split_1 <- if (input$split_by_1 == "None") NULL else input$split_by_1
    split_2 <- if (input$split_by_2 == "None") NULL else input$split_by_2
    plot_tax_stacked_bar(exps = tax_result[[3]], pheno = pheno_data(),
                         split_by_1 = split_1, split_by_2 = split_2,
                         show_names = input$show_names)
  })
  
  # ---- Downloads ----
  output$download_metadata <- downloadHandler(
    filename = function() { "metadata.csv" },
    content = function(file) {
      req(pheno_data())
      write_csv(pheno_data(), file)
    }
  )
  
  output$download_alpha_plot <- downloadHandler(
    filename = function() { "alpha_diversity_plot.png" },
    content = function(file) {
      req(alpha_data(), input$group_by)
      p <- plot_alpha_diversity(alpha_data(), input$group_by)
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  
  output$download_alpha_data <- downloadHandler(
    filename = function() { "alpha_diversity_data.csv" },
    content = function(file) {
      req(alpha_data())
      write_csv(alpha_data(), file)
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
  
  output$download_beta_data <- downloadHandler(
    filename = function() { "beta_diversity_data.csv" },
    content = function(file) {
      req(beta_data())
      write_csv(beta_data()$metadata, file)
    }
  )
  
  output$download_tax_plot <- downloadHandler(
    filename = function() { "taxonomy_bar_plot.png" },
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
    filename = function() { "taxonomy_data.csv" },
    content = function(file) {
      req(count_matrix_data())
      write_csv(count_matrix_data(), file)
    }
  )
}

shinyApp(ui, server)
