library(ggplot2)
library(tidyverse)
library(tidyquant)
library(purrr)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(reticulate)
library(shinythemes)
library(zoo)
library(shiny)
library(shinydashboard)
library(shinyjs)

print(paste("ulimit (in bytes) is", Cstack_info()["size"]))
print(paste("R version is", R.version.string))

# Source the code needed to make plots
source("plotting.R")

# Source the code needed to read from HDF
source("read_hdf.R")

data_folder = "../data/"

# Set up the user interface
ui <- dashboardPage(
  dashboardHeader(title = "GLAM"),
  dashboardSidebar(
    useShinyjs(),
    # Dropdown menu to select the dataset of interest
    selectInput(
      "dataset",
      "Dataset:",
      find_datasets(data_folder)
    )
  ),
  dashboardBody(
    useShinyjs(),
    fluidRow(
      box(
        width = 12,
        title="Dataset Summary",
        collapsible = TRUE,
        column(
          width = 4,
          DT::DTOutput("dataset_summary_DT")
        ),
        column(
          width = 8,
          div(
            style = "padding-right: 20px",
            fluidRow(
              tabsetPanel(
                type = "tabs",
                tabPanel(
                  "Total Reads",
                  plotOutput("total_reads_plot")
                ),
                tabPanel(
                  "Aligned Reads",
                  plotOutput("aligned_reads_plot")
                ),
                tabPanel(
                  "Prop. Aligned (Hist)",
                  plotOutput("proportion_aligned_plot")
                ),
                tabPanel(
                  "Prop. Aligned (Scatter)",
                  plotOutput("aligned_scatter_plot")
                ),
                tabPanel(
                  "PCA - Samples",
                  fluidRow(
                    uiOutput("color_pca_by")
                  ),
                  fluidRow(
                    plotOutput("pca_samples_plot")
                  )
                ),
                tabPanel(
                  "PCA - CAGs",
                  plotOutput("pca_cags_plot")
                )
              )
            ),
            fluidRow(
              fluidRow(
                div(
                  downloadButton("dataset_summary_pdf", label="PDF"),
                  downloadButton("dataset_summary_csv", label="CSV"),
                  style="text-align:right"
                )
              )
            )
          )
        )
      )
    ),
    fluidRow(
      box(
        width = 12,
        title="Statistical Analysis Results",
        collapsible = TRUE,
        column(
          fluidRow(uiOutput("parameter_select")),
          fluidRow(DT::DTOutput("corncob_results_DT")),
          width = 6
        ),
        column(
          fluidRow(plotOutput("corncob_results_plot")),
          fluidRow(
            div(
              downloadButton("corncob_results_pdf", label="PDF"),
              downloadButton("corncob_results_csv", label="CSV"),
              style="text-align:right"
            )
          ),
          width = 6
        )
      )
    ),
    fluidRow(
      box(
        width = 12,
        title="CAGs in Dataset",
        collapsible = TRUE,
        column(
          fluidRow(DT::DTOutput("cag_summary_DT")),
          width = 6
        ),
        column(
          fluidRow(
            column(
              width = 4,
              uiOutput("cag_summary_x_select")
            ),
            column(
              width = 4,
              uiOutput("cag_summary_y_select")
            ),
            column(
              width = 4,
              uiOutput("cag_summary_hue_select")
            )
          ),
          fluidRow(plotOutput("cag_summary_plot")),
          fluidRow(
            div(
              downloadButton("cag_summary_pdf", label="PDF"),
              downloadButton("cag_summary_csv", label="CSV"),
              style="text-align:right"
            )
          ),
          width = 6
        )
      )
    )
  )
)

panel_elements <- c("_DT", "_plot", "_pdf", "_csv")
show_panel <- function(panel_name){
  for(suffix in panel_elements){
    show(paste(panel_name, suffix, sep=""))
  }
}
hide_panel <- function(panel_name){
  for(suffix in panel_elements){
    hide(paste(panel_name, suffix, sep=""))
  }
}

server <- function(input, output) {
  
  # Watch for when the user changes the dataset
  observeEvent(input$dataset, {
    
    # If the new dataset contains corncob results, show the corncob menu and display
    if(read_hdf_has_corncob(input$dataset, data_folder)){
      
      show("parameter_select")
      show_panel("corncob_results")

    } else {
      
      # Otherwise, hide those elements
      hide("parameter_select")
      hide_panel("corncob_results")

    }
    
  })
  
  # Reactive element reading in the manifest for this dataset
  manifest_df <- reactive({read_hdf_manifest(input$dataset, data_folder)})
  
  # Reactive element with the list of parameters
  unique_parameters <- reactive({
    l <- unique(corncob_results_df()$parameter)
    return(rev(l))
  })
  
  # If the source HDF5 has corncob results, display the parameter list in the sidebar
  output$parameter_select <- renderUI({
    selectInput(
      "parameter",
      "Parameter:",
      unique_parameters()
    )
  })
  
  # If the source HDF5 has corncob results, display the parameter list in the sidebar
  output$color_pca_by <- renderUI({
    selectInput(
      "color_pca_by",
      "Color By:",
      manifest_df() %>%
        colnames
    )
  })

  # Reactive element recording whether the HDF5 contains corncob results
  has_corncob_results <- reactive({
    return(read_hdf_has_corncob(input$dataset, data_folder))
  })

  ########################
  # REACTIVE DATA TABLES #
  ########################
  
  # Reactive element with a table summarizing the overall dataset
  dataset_summary_df <- reactive({
    # Set up the table
    summary_df <- data.frame(Dataset=input$dataset)
    
    # Add the number of samples and CAGs
    summary_df["Number of Samples"] <- ncol(cag_abund_df())
    summary_df["Number of CAGs"] <- nrow(cag_abund_df())
    
    # Add the number of genes
    summary_df["Number of Genes"] <- nrow(gene_annotations_df())
    
    # Format the table for plotting
    summary_df <- t(summary_df)
    colnames(summary_df) <- c("Value")
    return(summary_df)
  })
  
  # Reactive element with the number of reads for this dataset
  readcounts_df <- reactive({
    read_hdf_readcounts(input$dataset, data_folder)
  })
  
  # Reactive element with the CAG abundances for the whole dataset
  cag_abund_df <- reactive({
    read_cag_abundances(input$dataset, data_folder)
  })
  
  # Reactive element with the annotations for every gene in this dataset
  gene_annotations_df <- reactive({
    read_gene_annotations(input$dataset, data_folder)
  })
  
  # Reactive element with the table of corncob results
  corncob_results_df <- reactive({
    if(has_corncob_results()){
      read_hdf_corncob_results(
        input$dataset, 
        data_folder
      )
    }else{
      data.frame(parameter=c('empty'), estimate=c('empty'), stdev=c('empty'))
    }
  })
  
  # Filter the corncob results by the selected parameter
  corncob_results_filtered_df <- reactive({
    if(is.null(input$parameter)){
      return(corncob_results_df())
    }else{
      return(corncob_results_df() %>% 
        filter(parameter == input$parameter) %>%
        select( -parameter ))
    }
  })
  
  # Manifest for the dataset
  manifest_df <- reactive({
    read_hdf_manifest(input$dataset, data_folder)
  })
  
  # Summary of all CAGs in this dataset
  cag_summary_df <- reactive({
    make_cag_summary(
      gene_annotations_df(),
      cag_abund_df(),
      corncob_results_filtered_df()
    )
  })
  
  ###################
  # DATASET SUMMARY #
  ###################
  
  # Render the dataset summary as a DataTable
  output$dataset_summary_DT <- DT::renderDT({
    dataset_summary_df()
  }, options = list(dom="t"))
  
  # Render the plots in the dataset summary tabset
  output$total_reads_plot <- renderPlot(plot_readcounts(readcounts_df(), "total"))
  output$aligned_reads_plot <- renderPlot(plot_readcounts(readcounts_df(), "aligned"))
  output$proportion_aligned_plot <- renderPlot(plot_readcounts(readcounts_df(), "proportion"))
  output$aligned_scatter_plot <- renderPlot(plot_readcounts(readcounts_df(), "scatter"))
  output$pca_samples_plot <- renderPlot(
    plot_dataset_pca(
      cag_abund_df(), manifest_df(), "Samples", input$color_pca_by
    )
  )
  output$pca_cags_plot <- renderPlot(
    plot_dataset_pca(
      cag_abund_df(), manifest_df(), "CAGs", input$color_pca_by
    )
  )
  
  # Make the dataset summary PDF available for download
  output$dataset_summary_pdf <- downloadHandler(
    filename = paste(
      input$dataset,
      "summary.pdf",
      sep="."
    ),
    content = function(file) {
      pdf(file)
      print(plot_readcounts(readcounts_df(), "total"))
      print(plot_readcounts(readcounts_df(), "aligned"))
      print(plot_readcounts(readcounts_df(), "proportion"))
      print(plot_readcounts(readcounts_df(), "scatter"))
      for(color_by in colnames(manifest_df())){
        print(plot_dataset_pca(cag_abund_df(), manifest_df(), "Samples", color_by))
      }
      print(plot_dataset_pca(cag_abund_df(), manifest_df(), "CAGs", input$color_pca_by))
      
      dev.off()
    }
  )
  
  # Make the dataset summary CSV available for download
  output$dataset_summary_csv <- downloadHandler(
    filename = paste(input$dataset, "summary.csv", sep="."),
    content = function(file) {
      write_csv(readcounts_df(), path = file)
    },
    contentType = "text/csv"
  )
    
  ###################
  # CORNCOB RESULTS #
  ###################
  
  # Render the corncob results as a DataTable
  output$corncob_results_DT <- DT::renderDT({
    corncob_results_filtered_df()
  }, rownames = FALSE)
  
  # Render the corncob results as a scatter plot
  output$corncob_results_plot <- renderPlot(
    plot_corncob_results(corncob_results_df(), input$parameter)
  )
  
  # Make the corncob results PDF available for download
  output$corncob_results_pdf <- downloadHandler(
    filename = paste(input$dataset, "corncob.pdf", sep="."),
    content = function(file) {
      ggsave(
        file, 
        plot = plot_corncob_results(corncob_results_df(), input$parameter), 
        device="pdf"
      )
    }
  )
  
  # Make the corncob results CSV available for download
  output$corncob_results_csv <- downloadHandler(
    filename = paste(input$dataset, "corncob.csv", sep="."),
    content = function(file) {
      write_csv(corncob_results_df(), path = file)
    },
    contentType = "text/csv"
  )
  
  ###############
  # CAG SUMMARY #
  ###############
  
  # Render the CAG summary as a DataTable
  output$cag_summary_DT <- DT::renderDT({
    cag_summary_df()
  }, rownames = FALSE)
  
  # Provide menus to let the user select what kind of plot to make
  output$cag_summary_x_select <- renderUI({
    selectInput(
      "cag_summary_x_select",
      "Horizontal Axis:",
      cag_summary_df() %>%
        colnames,
      selected = "Max Abund"
    )
  })
  output$cag_summary_y_select <- renderUI({
    selectInput(
      "cag_summary_y_select",
      "Vertical Axis:",
      cag_summary_df() %>%
        colnames,
      selected = "Prevalence"
    )
  })
  output$cag_summary_hue_select <- renderUI({
    selectInput(
      "cag_summary_hue_select",
      "Color:",
      cag_summary_df() %>%
        colnames,
      selected = "Number of Genes"
    )
  })
  
  # Render the CAG summary as a scatter plot
  output$cag_summary_plot <- renderPlot(
    plot_cag_summary(
      cag_summary_df(), 
      input$cag_summary_x_select, 
      input$cag_summary_y_select, 
      input$cag_summary_hue_select
    )
  )
  
  # Make the CAG summary PDF available for download
  output$cag_summary_pdf <- downloadHandler(
    filename = paste(input$dataset, "CAG.summary.pdf", sep="."),
    content = function(file) {
      ggsave(
        file, 
        plot = plot_cag_summary(
          cag_summary_df(), 
          input$cag_summary_x_select, 
          input$cag_summary_y_select, 
          input$cag_summary_hue_select
        ), 
        device="pdf"
      )
    }
  )
  
  # Make the CAG summary CSV available for download
  output$cag_summary_csv <- downloadHandler(
    filename = paste(input$dataset, "CAG.summary.csv", sep="."),
    content = function(file) {
      write_csv(cag_summary_df(), path = file)
    },
    contentType = "text/csv"
  )
}

shinyApp(ui, server)

