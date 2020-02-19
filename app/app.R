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
                  "Size",
                  plotOutput("cag_size_histogram")
                ),
                tabPanel(
                  "Prevalence",
                  plotOutput("cag_prevalence_histogram")
                ),
                tabPanel(
                  "Combined",
                  plotOutput("cag_size_prevalence_distribution")
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
        title="Sample Ordination",
        collapsible = TRUE,
        column(
          width = 4,
          uiOutput("color_ordination_by")
        ),
        column(
          width = 8,
          div(
            style = "padding-right: 20px",
            fluidRow(
              tabsetPanel(
                type = "tabs",
                tabPanel(
                  "PCA",
                  plotOutput("ordination_pca_scatter")
                ),
                tabPanel(
                  "t-SNE",
                  plotOutput("ordination_tsne_scatter")
                )
              )
            ),
            fluidRow(
              fluidRow(
                div(
                  downloadButton("ordination_pdf", label="PDF"),
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
        title="Statistical Analysis",
        collapsible = TRUE,
        column(
          fluidRow(uiOutput("parameter_select")),
          width = 4
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
          width = 8
        )
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "CAGs Summary Table",
        collapsible = TRUE,
        fluidRow(div(
          style = "padding-right: 20px; padding-left: 20px",
          DT::DTOutput("cag_summary_DT")
        )),
        fluidRow(div(
          downloadButton("cag_summary_csv", label="CSV"),
          style="text-align:right; padding-right: 20px; padding-top: 10px"
        ))
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
    l <- l[l != "(Intercept)"]
    return(l)
  })
  
  # If the source HDF5 has corncob results, display the parameter list in the sidebar
  output$parameter_select <- renderUI({
    selectInput(
      "parameter",
      "Parameter:",
      unique_parameters()
    )
  })
  
  # Display the list of options for coloring samples in the ordination plot
  output$color_ordination_by <- renderUI({
    selectInput(
      "color_ordination_by",
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
    summary_df["Number of Samples"] <- nrow(manifest_df())
    summary_df["Number of CAGs"] <- nrow(cag_summary_df())
    
    # Add the number of genes
    summary_df["Number of Genes"] <- sum(cag_summary_df()["size"])
    
    # Format the table for plotting
    summary_df <- t(summary_df)
    colnames(summary_df) <- c("Value")
    return(summary_df)
  })
  
  # Reactive element with the number of reads for this dataset
  readcounts_df <- reactive({
    read_hdf_readcounts(input$dataset, data_folder)
  })
  
  # Reactive element with the PCA results
  pca_df <- reactive({
    read_hdf_pca(input$dataset, data_folder)
  })
  
  # Reactive element with the t-SNE results
  tsne_df <- reactive({
    read_hdf_tsne(input$dataset, data_folder)
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
    read_hdf_cag_summary(input$dataset, data_folder)
  })
  
  # Extended CAG summary table containing corncob result metrics
  cag_extended_summary_df <- reactive({
    if(has_corncob_results() && is.null(input$parameter) == FALSE){
      cag_summary_df() %>%
        full_join(corncob_results_filtered_df()) %>%
        rename(
          `Mean Abundance` = mean_abundance,
          Prevalence = prevalence,
          `Number of Genes` = size,
          `Estimated Coefficient` = estimate,
          `Std. Error` = std_error,
          `p-value` = p_value
        )
    } else {
      cag_summary_df() %>%
        rename(
          `Mean Abundance` = mean_abundance,
          Prevalence = prevalence,
          `Number of Genes` = size
        )
    }
  })
  
  ###################
  # DATASET SUMMARY #
  ###################
  
  # Render the dataset summary as a DataTable
  output$dataset_summary_DT <- DT::renderDT({
    dataset_summary_df()
  }, options = list(dom="t"))
  
  # Render the plots in the dataset summary tabset
  output$cag_size_histogram <- renderPlot(plot_cag_hist(cag_summary_df(), "size"))
  output$cag_prevalence_histogram <- renderPlot(plot_cag_hist(cag_summary_df(), "prevalence"))
  output$cag_size_prevalence_distribution <- renderPlot(plot_cag_size_prevalence_distribution(cag_summary_df()))

  # Make the dataset summary PDF available for download
  output$dataset_summary_pdf <- downloadHandler(
    filename = paste(
      input$dataset,
      "summary.pdf",
      sep="."
    ),
    content = function(file) {
      pdf(file)
      print(plot_cag_hist(cag_summary_df(), "size"))
      print(plot_cag_hist(cag_summary_df(), "prevalence"))
      print(plot_cag_size_prevalence_distribution(cag_summary_df()))
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
  
  ##############
  # ORDINATION #
  ##############
  
  output$ordination_pca_scatter <- renderPlot(plot_ordination_scatter(
    pca_df(), manifest_df(), input$color_ordination_by, "PCA"
  ))
  output$ordination_tsne_scatter <- renderPlot(plot_ordination_scatter(
    tsne_df(), manifest_df(), input$color_ordination_by, "t-SNE"
  ))
  output$ordination_pdf <- downloadHandler(
    filename = paste(
      input$dataset,
      "ordination.pdf",
      sep="."
    ),
    content = function(file) {
      pdf(file)
      for(color_by in colnames(manifest_df())){
        print(plot_ordination_scatter(
          pca_df(), manifest_df(), color_by, "PCA"
        ))
      }
      for(color_by in colnames(manifest_df())){
        print(plot_ordination_scatter(
          tsne_df(), manifest_df(), color_by, "t-SNE"
        ))
      }
      dev.off()
    }
  )

  ###################
  # CORNCOB RESULTS #
  ###################


  # Render the corncob results as a scatter plot
  output$corncob_results_plot <- renderPlot(
    plot_corncob_results(corncob_results_df(), input$parameter)
  )

  # Make the corncob results PDF available for download
  output$corncob_results_pdf <- downloadHandler(
    filename = paste(input$dataset, "corncob.pdf", sep="."),
    content = function(file) {
      pdf(file)
      for(parameter_name in unique_parameters()){
        print(plot_corncob_results(corncob_results_df(), parameter_name))
      }
      dev.off()
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
    cag_extended_summary_df()
  }, rownames = FALSE)
  
  # Make the CAG summary CSV available for download
  output$cag_summary_csv <- downloadHandler(
    filename = paste(input$dataset, "CAG.summary.csv", sep="."),
    content = function(file) {
      write_csv(cag_extended_summary_df(), path = file)
    },
    contentType = "text/csv"
  )
}

shinyApp(ui, server)

