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
source_python("read_hdf.py")

data_folder = Sys.getenv("DATA_DIR")

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
                  "Abundance",
                  plotOutput("cag_abundance_histogram")
                ),
                tabPanel(
                  "Combined",
                  plotOutput("cag_size_abundance_distribution")
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
        title = "CAG Summary Table",
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
    ),
    fluidRow(
      box(
        width = 12,
        title = "Selected CAG Details",
        collapsible = TRUE,
        fluidRow(div(
          style = "padding-right: 20px; padding-left: 20px",
          column(
            width = 8,
            DT::DTOutput("cag_details_DT")
          ),
          column(
            width = 4,
            align = "center",
            plotOutput("cag_details_tax_bars", width = "300px", height = "600px")
          )
        )),
        fluidRow(div(
          downloadButton("cag_details_csv", label="CSV"),
          downloadButton("cag_details_pdf", label="PDF"),
          style="text-align:right; padding-right: 20px; padding-top: 10px"
        ))
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "Selected CAG Abundance",
        collapsible = TRUE,
        column(
          width = 4,
          uiOutput("cag_abundance_x"),
          uiOutput("cag_abundance_hue"),
          uiOutput("cag_abundance_col"),
          uiOutput("cag_abundance_group"),
          uiOutput("cag_abundance_geom"),
          checkboxInput(
            "cag_abundance_geom_logy", 
            "Log Scale",
            value = TRUE
          )
        ),
        column(
          width = 8,
          fluidRow(div(
            plotOutput("cag_abundance_plot")
          )),
          fluidRow(div(
            downloadButton("cag_abundance_csv", label="CSV"),
            downloadButton("cag_abundance_pdf", label="PDF"),
            style="text-align:right; padding-right: 20px; padding-top: 10px"
          ))
        )
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "CAG Abundance Heatmap",
        collapsible = TRUE,
        column(
          width = 4,
          selectizeInput(
            "heatmap_cag_selector",
            "CAGs",
            NULL,
            multiple = TRUE
          ),
          uiOutput("heatmap_color_by_primary"),
          uiOutput("heatmap_color_by_secondary"),
          uiOutput("heatmap_color_by_tertiary"),
          fluidRow(
            column(
              width = 6,
              checkboxInput(
                "heatmap_checkbox_sample_name", 
                "Sample Names",
                value = TRUE
              )
            ),
            column(
              width = 6,
              checkboxInput(
                "heatmap_checkbox_cluster", 
                "Cluster"
              )
            )
          )
        ),
        column(
          width = 8,
          fluidRow(div(
            plotOutput("heatmap_cag_plot")
          )),
          fluidRow(div(
            downloadButton("heatmap_cag_csv", label="CSV"),
            downloadButton("heatmap_cag_pdf", label="PDF"),
            style="text-align:right; padding-right: 20px; padding-top: 10px"
          ))
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

server <- function(input, output, session) {
  
  # Watch for when the user changes the dataset
  observeEvent(input$dataset, {
    
    # Read all of the data objects while showing status updates
    withProgress(message = 'Reading data', value = 0, {
      
      # Number of objects to read
      n = 8
      # Read in each table in turn
      dummy <- manifest_df()
      # And update the progress bar for each one
      incProgress(1/n, detail = "Read in metadata table")
      # etc.
      dummy <- corncob_results_df()
      incProgress(1/n, detail = "Read in corncob results")
      dummy <- dataset_summary_df()
      incProgress(1/n, detail = "Read in dataset summary table")
      dummy <- readcounts_df()
      incProgress(1/n, detail = "Read in aligned reads table")
      dummy <- pca_df()
      incProgress(1/n, detail = "Read in PCA table")
      dummy <- tsne_df()
      incProgress(1/n, detail = "Read in t-SNE table")
      dummy <- cag_summary_df()
      dummy <- cag_extended_summary_df()
      incProgress(1/n, detail = "Read in CAG summary table")
      dummy <- cag_details_df()
      incProgress(1/n, detail = "Read in CAG detail table")
    })

    # If the new dataset contains corncob results, show the corncob menu and display
    if(read_hdf_has_corncob(input$dataset, data_folder)){

      show("parameter_select")
      show_panel("corncob_results")

    } else {

      # Otherwise, hide those elements
      hide("parameter_select")
      hide_panel("corncob_results")

    }
    
    # Update the multiple CAG selector with the CAGs for this dataset
    updateSelectizeInput(
      session, 
      'heatmap_cag_selector', 
      choices = cag_summary_df()$CAG, 
      server = TRUE
    )
    
    # If the dataset has more than 15 samples, don't plot the sample names by default in the heatmap
    updateCheckboxInput(
      session, 
      "heatmap_checkbox_sample_name", 
      value = nrow(manifest_df()) <= 15
    )
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
    # Read the summary table directly from the HDF5
    withProgress(message = "Reading dataset summary", {
      return(
        read_hdf_summary(input$dataset, data_folder) %>%
          remove_rownames %>%
          column_to_rownames("variable")
      )
    })
  })
  
  # Reactive element with the number of reads for this dataset
  readcounts_df <- reactive({
    withProgress(
      message = "Reading alignment metrics",
      {read_hdf_readcounts(input$dataset, data_folder)}
    )
  })
  
  # Reactive element with the PCA results
  pca_df <- reactive({
    withProgress(
      message = "Reading PCA",
      {read_hdf_pca(input$dataset, data_folder)}
    )
  })
  
  # Reactive element with the t-SNE results
  tsne_df <- reactive({
    withProgress(
      message = "Reading t-SNE",
      {read_hdf_tsne(input$dataset, data_folder)}
    )
  })

  # Reactive element with the table of corncob results
  corncob_results_df <- reactive({
    if(has_corncob_results()){
      withProgress(
        message = "Reading statistical analysis results",
        {read_hdf_corncob_results(
          input$dataset, 
          data_folder
        )}
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
    withProgress(
      message = "Reading manifest",
      {read_hdf_manifest(input$dataset, data_folder)}
    )
  })
  
  # Summary of all CAGs in this dataset
  cag_summary_df <- reactive({
    withProgress(
      message = "Reading CAG summary",
      {read_hdf_cag_summary(input$dataset, data_folder)}
    )
  })
  
  # Extended CAG summary table containing corncob result metrics
  cag_extended_summary_df <- reactive({
    if(has_corncob_results() && is.null(input$parameter) == FALSE){
      cag_summary_df() %>%
        full_join(corncob_results_filtered_df()) %>%
        arrange(
          p_value
        ) %>%
        drop_na(p_value) %>% # Drop CAGs with no p-values computed
        mutate( # Truncate significant digits
          estimate = estimate %>% signif(4),
          p_value = p_value %>% signif(4),
          std_error = std_error %>% signif(4)
        ) %>%
        rename(
          `Mean Abundance` = mean_abundance,
          `Mean Abundance (log10)` = mean_abundance_log10,
          Prevalence = prevalence,
          `Number of Genes` = size,
          `Number of Genes (log10)` = size_log10,
          `Estimated Coefficient` = estimate,
          `Std. Error` = std_error,
          `p-value` = p_value
        )
    } else {
      cag_summary_df() %>%
        rename(
          `Mean Abundance` = mean_abundance,
          `Mean Abundance (log10)` = mean_abundance_log10,
          Prevalence = prevalence,
          `Number of Genes` = size,
          `Number of Genes (log10)` = size_log10
        )
    }
  })
  
  # Summary table for the single selected CAG
  cag_details_df <- reactive({
    if(is.null(input$cag_summary_DT_rows_selected)){
      return(data.frame())
    } else {
      # Get the ID for the CAG which has been selected
      cag_id <- cag_extended_summary_df()$CAG[input$cag_summary_DT_rows_selected]
      return(
        read_cag_details(
          file.path(
            data_folder, 
            paste(input$dataset, ".hdf5", sep="")
          ), 
          cag_id
        )
      )
    }
  })
  
  # Relative abundance table for the single selected CAG
  cag_abundance_df <- reactive({
    if(is.null(input$cag_summary_DT_rows_selected)){
      return(data.frame())
    } else {
      # Get the ID for the CAG which has been selected
      cag_id <- cag_extended_summary_df()$CAG[input$cag_summary_DT_rows_selected]
      return(
        read_cag_abundance(
          file.path(
            data_folder, 
            paste(input$dataset, ".hdf5", sep="")
          ), 
          cag_id
        ) %>% 
          as_tibble
        ) %>% full_join(
          manifest_df()
        )
      }
  })
  
  # Relative abundance table for multiple selected CAGs
  multiple_cag_abundance_df <- reactive({
    if(length(input$heatmap_cag_selector) == 0){
      return(data.frame())
    } else {
      return(
        read_multiple_cag_abundances(
          file.path(
            data_folder, 
            paste(input$dataset, ".hdf5", sep="")
          ), 
          input$heatmap_cag_selector
        ) %>% 
          as_tibble
      ) %>% full_join(
        manifest_df()
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
  output$cag_abundance_histogram <- renderPlot(plot_cag_hist(cag_summary_df(), "abundance"))
  output$cag_size_abundance_distribution <- renderPlot(plot_cag_size_abundance_distribution(cag_summary_df()))

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
      print(plot_cag_hist(cag_summary_df(), "abundance"))
      print(plot_cag_size_abundance_distribution(cag_summary_df()))
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
  output$cag_summary_DT <- DT::renderDT(
    {cag_extended_summary_df()}, 
    rownames = FALSE,
    server = TRUE,
    selection = list(
      mode = 'single',
      selected = 1
    )
  )
  
  # Make the CAG summary CSV available for download
  output$cag_summary_csv <- downloadHandler(
    filename = paste(input$dataset, "CAG.summary.csv", sep="."),
    content = function(file) {
      write_csv(cag_extended_summary_df(), path = file)
    },
    contentType = "text/csv"
  )
  
  ###############
  # CAG DETAILS #
  ###############
  
  # Render the details for a single CAG as a DataTable
  output$cag_details_DT <- DT::renderDT(
    {cag_details_df()},
    rownames = FALSE,
    server = TRUE,
    selection = 'none'
  )
  
  # Render the assigned taxa for a CAG as a bargraph
  output$cag_details_tax_bars <- renderPlot(
    {plot_cag_details_tax_bars(
      cag_details_df(),
      cag_extended_summary_df()$CAG[input$cag_summary_DT_rows_selected]
    )},
    width = 300,
    height = 600
  )
  
  # Make the taxonomic assignment barplot for a single CAG available for download as a PDF
  output$cag_details_pdf <- downloadHandler(
    filename = paste(input$dataset, "CAG.taxa.pdf", sep="."),
    content = function(file) {
      pdf(file)
      print(plot_cag_details_tax_bars(
        cag_details_df(),
        cag_extended_summary_df()$CAG[input$cag_summary_DT_rows_selected]
      ))
      dev.off()
    }
  )

  # Make the details for a single CAG available for download as a CSV
  output$cag_details_csv <- downloadHandler(
    filename = paste(input$dataset, "CAG.details.csv", sep="."),
    content = function(file) {
      write_csv(cag_details_df(), path = file)
    },
    contentType = "text/csv"
  )
  
  #################
  # CAG ABUNDANCE #
  #################
  
  output$cag_abundance_x <- renderUI({
    selectInput(
      "cag_abundance_x",
      "Horizontal x-axis:",
      manifest_df() %>%
        colnames
    )
  })
  output$cag_abundance_hue <- renderUI({
    selectInput(
      "cag_abundance_hue",
      "Hue:",
      c("None", manifest_df() %>% colnames)
    )
  })
  output$cag_abundance_col <- renderUI({
    selectInput(
      "cag_abundance_col",
      "Column:",
      c("None", manifest_df() %>% colnames)
    )
  })
  output$cag_abundance_group <- renderUI({
    selectInput(
      "cag_abundance_group",
      "Group By:",
      c("None", manifest_df() %>% colnames)
    )
  })
  output$cag_abundance_geom <- renderUI({
    selectInput(
      "cag_abundance_geom",
      "Plot Type:",
      c("point", "boxplot", "jitter", "line", "bar")
    )
  })
  output$cag_abundance_plot <- renderPlot(
    plot_cag_abundance(
      cag_abundance_df(),
      input$cag_abundance_x,
      input$cag_abundance_hue,
      input$cag_abundance_col,
      input$cag_abundance_group,
      input$cag_abundance_geom,
      input$cag_abundance_geom_logy
    )
  )
  
  # Make the single CAG abundance table available for download as a PDF
  output$cag_abundance_pdf <- downloadHandler(
    filename = paste(input$dataset, "CAG.abundance.pdf", sep="."),
    content = function(file) {
      pdf(file)
      print(
        plot_cag_abundance(
          cag_abundance_df(),
          input$cag_abundance_x,
          input$cag_abundance_hue,
          input$cag_abundance_col,
          input$cag_abundance_group,
          input$cag_abundance_geom,
          input$cag_abundance_geom_logy
        )
      )
      dev.off()
    }
  )
  
  # Make the single CAG abundance table available for download as a CSV
  output$cag_details_csv <- downloadHandler(
    filename = paste(input$dataset, "CAG.abundance.csv", sep="."),
    content = function(file) {
      write_csv(cag_abundance_df(), path = file)
    },
    contentType = "text/csv"
  )
  
  ########################
  # MULTIPLE CAG HEATMAP #
  ########################
  output$heatmap_color_by_primary <- renderUI({
    selectInput(
      "heatmap_color_by_primary",
      "Color By (primary):",
      c("NONE", colnames(manifest_df()))
    )
  })
  output$heatmap_color_by_secondary <- renderUI({
    selectInput(
      "heatmap_color_by_secondary",
      "Color By (secondary):",
      c("NONE", colnames(manifest_df()))
    )
  })
  output$heatmap_color_by_tertiary <- renderUI({
    selectInput(
      "heatmap_color_by_tertiary",
      "Color By (tertiary):",
      c("NONE", colnames(manifest_df()))
    )
  })
  output$heatmap_cag_plot <- renderPlot({
    plot_cag_heatmap(
      multiple_cag_abundance_df(),
      input$heatmap_cag_selector,
      input$heatmap_color_by_primary,
      input$heatmap_color_by_secondary,
      input$heatmap_color_by_tertiary,
      input$heatmap_checkbox_sample_name,
      input$heatmap_checkbox_cluster
    )
  })
  
  # Make the multiple CAG abundance table available for download as a PDF
  output$heatmap_cag_pdf <- downloadHandler(
    filename = paste(input$dataset, "multiple.CAG.abundance.pdf", sep="."),
    content = function(file) {
      pdf(file)
      print(
        plot_cag_heatmap(
          multiple_cag_abundance_df(),
          input$heatmap_cag_selector,
          input$heatmap_color_by_primary,
          input$heatmap_color_by_secondary,
          input$heatmap_color_by_tertiary,
          input$heatmap_checkbox_sample_name,
          input$heatmap_checkbox_cluster
        )
      )
      dev.off()
    }
  )
  
  # Make the single CAG abundance table available for download as a CSV
  output$heatmap_cag_csv <- downloadHandler(
    filename = paste(input$dataset, "multiple.CAG.abundance.csv", sep="."),
    content = function(file) {
      write_csv(multiple_cag_abundance_df(), path = file)
    },
    contentType = "text/csv"
  )
  
}

shinyApp(ui, server)

