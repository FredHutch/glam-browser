# This function will make a volcano plot for one of the metadata
# features included by the user in the corncob results
plot_corncob_results <- function(corncob_results_df, parameter){
  
  plot_df <- corncob_results_df
  
  # If the parameter is NULL, return an empty plot
  if(!is.null(parameter) && "p_value" %in% colnames(plot_df)){
    # Filter to this parameter
    plot_df <- plot_df[plot_df$parameter == parameter,]

    # Add the -log10 p-value
    plot_df <- plot_df %>%
      add_column(
        neg_log_pvalue = sapply(
          pull(plot_df, p_value),
          function(v){-log10(v)}
        )
      )

    p <- ggplot(
      data = plot_df,
      aes(
        x = estimate,
        y = neg_log_pvalue
      )
    ) + geom_point(
    ) + geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color="grey"
    ) + geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color="grey"
    ) + ggtitle(
      parameter
    ) + xlab(
      "Estimated Coefficient"
    ) + ylab(
      "p-value (-log10)"
    ) + theme_minimal(
    ) + theme(
      plot.title = element_text(hjust = 0.5)
    )
    
    return(p)
  } else {
    return(ggplot() + geom_blank())
  }
  
}

# Function to plot readcounts
plot_readcounts <- function(readcounts_df, plot_type){
  if(plot_type == "total"){
    p <- ggplot(
      data = readcounts_df, aes(x = n_reads)
    ) + geom_histogram() + xlab(
      "Number of Total Reads"
    )
  }else{
    if(plot_type == "aligned"){
      p <- ggplot(
        data = readcounts_df, aes(x = aligned_reads)
      ) + geom_histogram() + xlab(
        "Number of Aligned Reads"
      )
    } else {
      if(plot_type == "proportion"){
        p <- ggplot(
          data = readcounts_df, aes(x = prop_aligned)
        ) + geom_histogram() + xlab(
          "Proportion of Aligned Reads"
        ) + xlim(
          0, 1
        )
      } else {
        stopifnot(plot_type == "scatter")
        p <- ggplot(
          data = readcounts_df, aes(x = n_reads, y=aligned_reads)
        ) + geom_point() + xlab(
          "Number of Total Reads"
        ) + ylab(
          "Number of Aligned Reads"
        )
      }
    }
  }
  return(
    p + ylab(
      "Number of Samples"
    ) + theme_minimal()
  )
}

# Function to plot PCA for CAGs or samples
plot_dataset_pca <- function(cag_abund_df, manifest_df, dataset_summary_plot_type, color_pca_by){
  
  if(is.null(dataset_summary_plot_type) || is.null(color_pca_by)){
    return(ggplot() + geom_blank())
  }
  
  if( dataset_summary_plot_type == "Samples" ){

    # Run PCA
    pca <- prcomp(
      cag_abund_df
    )$rotation %>%
      as_tibble
    pca["specimen"] <- colnames(cag_abund_df)
    pca <- pca %>% full_join(
        manifest_df
      )
    
    p <- ggplot(
      data = pca,
      aes(
        x = PC1,
        y = PC2,
        color = pca[[color_pca_by]]
      )
    ) + geom_point(
    ) + guides(
      color=guide_legend(title=color_pca_by)
    )
    
    
  } else {
    stopifnot(dataset_summary_plot_type == "CAGs")
    
    # Run PCA
    pca <- prcomp(
      cag_abund_df %>% t
    )$rotation %>%
      as_tibble
    
    p <- ggplot(
      data = pca,
      aes(
        x = PC1,
        y = PC2
      )
    ) + geom_point(
    )
    
  }
  
  p <- p + ggtitle(
    paste("PCA -", dataset_summary_plot_type)
  ) + theme_minimal(
  ) + theme(
    plot.title = element_text(hjust = 0.5)
  )
  
  return(p)
  
}

# Function to make a plot based on various summary metrics for each CAG
plot_cag_summary <- function(cag_summary_df, x_val, y_val, hue_val){
  
  p <- ggplot(
    data = cag_summary_df,
    aes(
      x = cag_summary_df[[x_val]],
      y = cag_summary_df[[y_val]],
      color = cag_summary_df[[hue_val]]
    )
  ) + geom_point(
  ) + ggtitle(
    "CAG Summary"
  ) + guides(
    color=guide_legend(title=hue_val)
  ) + xlab(
    x_val
  ) + ylab(
    y_val
  ) + theme_minimal(
  ) + theme(
    plot.title = element_text(hjust = 0.5)
  )
  
  return(p)
}

# This function will plot the abundance of a single CAG across
# all samples, as a function of different aspects of sample metadata
# The intent here is that the user has quite a bit of flexibility
# in how to display and visualize the data, depending on the
# experimental design and question that they are asking
plot_cag_abundance <- function(cag_id, x, hue, group, col, row){

    # Read the abundance for this CAG, along with the sample metadata
    plot_df <- read_cag_abundance(cag_id)
    
}