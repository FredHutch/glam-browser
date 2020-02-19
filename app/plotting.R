
# Function to make a plot summarizing CAGs on the basis of size or prevalence
plot_cag_hist <- function(cag_summary_df, plot_type){

  if(plot_type == "size"){
    p <- ggplot(
      data = cag_summary_df,
      aes(
        x = size
      )
    ) + ggtitle(
      "CAG Size"
    ) + xlab(
      "Number of Genes"
    )
  } else {
    stopifnot(plot_type == "prevalence")
    p <- ggplot(
      data = cag_summary_df,
      aes(
        x = prevalence
      )
    ) + ggtitle(
      "CAG Prevalence"
    ) + xlab(
      "Proportion of Samples"
    )
  }
  return(
    p + geom_histogram(
    ) + xlab(
      paste("CAG", plot_type)
    ) + ylab(
      "Number of CAGs per bin"
    ) + theme_minimal(
    ) + theme(
      plot.title = element_text(hjust = 0.5)
    )
  )
}

plot_cag_size_prevalence_distribution <- function(cag_summary_df){
  p <- ggplot(
    data = cag_summary_df,
    aes(
      x = size,
      y = prevalence
    )
  ) + ggtitle(
    "CAG Size and Prevalence"
  ) + xlab(
    "Number of Genes"
  ) + ylab(
    "Proportion of Samples"
  ) + geom_hex(
  ) + theme_minimal(
  ) + theme(
    plot.title = element_text(hjust = 0.5)
  )
  return(p)
}


# Function to plot PCA for samples
plot_ordination_scatter <- function(ordination_df, manifest_df, color_ordination_by, title_text){
  
  if(is.null(color_ordination_by)){
    return(ggplot() + geom_blank())
  }

  # Join the annotations to the PCA table
  ordination_df <- ordination_df %>% full_join(
      manifest_df
    )
  
  p <- ggplot(
    data = ordination_df,
    aes(
      x = ordination_df[[names(ordination_df)[2]]],
      y = ordination_df[[names(ordination_df)[3]]],
      color = ordination_df[[color_ordination_by]]
    )
  ) + geom_point(
  ) + guides(
    color=guide_legend(title=color_ordination_by)
  ) + ggtitle(
    title_text
  ) + xlab(
    names(ordination_df)[2]
  ) + ylab(
    names(ordination_df)[3]
  ) + theme_minimal(
  ) + theme(
    plot.title = element_text(hjust = 0.5)
  )
  
  return(p)
  
}

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

plot_cag_abundance <- function(
  cag_abundance_df,
  x,
  hue,
  col,
  geom
){
  if(nrow(cag_abundance_df) == 0){
    return(ggplot() + geom_blank())
  } else {
    
    if(col != "None"){
      cag_abundance_df <- cag_abundance_df %>% rename(
        col_variable = col
      )
    }
      
    if(hue == "None"){
      p <- ggplot(
        data = cag_abundance_df,
        aes_string(
          x = x,
          y = "abundance"
        )
      )
    } else {
      p <- ggplot(
        data = cag_abundance_df,
        aes_string(
          x = x,
          y = "abundance",
          fill = hue
        )
      )
    }
    if(col != "None"){
      p <- p + facet_wrap(
        vars(col_variable), nrow = 1
      )
    }
    p <- p + xlab(
      x
    ) + ylab(
      "Relative Abundance"
    )
    if(geom == "point"){
      p <- p + geom_point()
    }
    if(geom == "boxplot"){
      p <- p + geom_boxplot()
    }
    if(geom == "jitter"){
      p <- p + geom_jitter()
    }
    if(geom == "line"){
      p <- p + geom_line()
    }
    if(geom == "bar"){
      p <- p + geom_bar()
    }
    return(
      p + theme_minimal(
      ) + theme(
        plot.title = element_text(hjust = 0.5)
      )
    )
  }
}