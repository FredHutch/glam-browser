library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# Function to plot the number of reads per sample
plot_sample_summary <- function(df, x_val, plot_type, show_ylabels){
  
  if(plot_type == "Values"){
    g <- ggplot(
      data = df,
      aes(
        x = df$Specimen,
        y = df[[x_val]]
      )
    ) + geom_bar(
      stat = "identity"
    ) + xlab(
      "Specimen Name"
    ) + ylab(
      x_val
    ) + coord_flip(
    )
  } else {
    if(plot_type == "Histogram"){
      g <- ggplot(
        data = df,
        aes(
          x = df[[x_val]]
        )
      ) + geom_histogram(
      ) + xlab(
        x_val
      ) + ylab(
        "Number of specimens per bin"
      )
    } else {
      stopifnot(plot_type == "Scatter vs. Raw Reads")
      g <- ggplot(
        data = df,
        aes(
          x = df[["Raw Reads"]],
          y = df[[x_val]]
        )
      ) + geom_point(
      ) + ylab(
        x_val
      ) + xlab(
        "Raw Reads"
      )
    }
  }
  if(plot_type == "Values" && show_ylabels == FALSE){
    return(
      g + theme_minimal(
        base_size = 20
      ) + theme(
        axis.text.y = element_blank()
      )
    )
  } else {
    return(
      g + theme_minimal(
        base_size = 20
      )
    )
  }
}

# Function to make a plot summarizing CAGs on the basis of size or abundance
plot_cag_hist <- function(cag_summary_df, plot_type){
  
  if(plot_type == "size"){
    p <- ggplot(
      data = cag_summary_df,
      aes(
        x = size_log10,
        weight = size
      )
    ) + ggtitle(
      "CAG Size"
    ) + xlab(
      "Number of Genes per CAG (log10)"
    )
  } else {
    stopifnot(plot_type == "abundance")
    p <- ggplot(
      data = cag_summary_df,
      aes(
        x = mean_abundance_log10,
        weight = size
      )
    ) + ggtitle(
      "CAG Abundance"
    ) + xlab(
      "Mean Abundance per CAG (log10)"
    )
  }
  return(
    p + geom_histogram(
      bins = 60
    ) + ylab(
      "Number of Genes per Bin"
    ) + theme_minimal(
      base_size = 20
    ) + theme(
      plot.title = element_text(hjust = 0.5)
    )
  )
}

plot_cag_size_abundance_distribution <- function(cag_summary_df){
  p <- ggplot(
    data = cag_summary_df %>%
      filter(size >= 5),
    aes(
      x = size_log10,
      y = mean_abundance_log10,
      weight = size
    )
  ) + ggtitle(
    "CAG Size and Abundance"
  ) + xlab(
    "Number of Genes (log10)"
  ) + ylab(
    "Mean Abundance (log10)"
  ) + geom_hex(
    bins = c(100, 30)
  ) + labs(
    fill = "Number of Genes per Bin"
  ) + theme_minimal(
    base_size = 20
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
    base_size = 20
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
    
    # Filter to p-values < 0.99
    plot_df <- plot_df %>% filter(p_value < 0.99)
    
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
      base_size = 20
    ) + theme(
      plot.title = element_text(hjust = 0.5)
    )
    
    return(p)
  } else {
    return(ggplot() + geom_blank())
  }
  
}

plot_cag_details_tax_bars <- function(
  cag_details_df,
  cag_name
){
  if(nrow(cag_details_df) == 0){
    return(ggplot() + geom_blank())
  }
  if(! "tax_name" %in% colnames(cag_details_df)){
    return(ggplot() + geom_blank())
  }
  g <- ggplot(
    data = cag_details_df %>%
      filter(
        tax_id != "0"
      ) %>%
      count(tax_name) %>%
      arrange(n) %>%
      tail(10),
    aes(
      x = reorder(tax_name, -n),
      y = n
    )
  ) + geom_bar(
    stat = "identity"
  ) + ylab(
    "Number of annotated genes"
  ) + xlab(
    ""
  ) + ggtitle(
    paste("CAG", cag_name)
  ) + theme_minimal(
    base_size = 20
  ) + theme(
    aspect.ratio=2,
    plot.title = element_text(hjust = 0.5)
  )
  return(g)
}

plot_cag_abundance <- function(
  cag_abundance_df,
  x,
  hue,
  col,
  group,
  geom,
  logy
){
  if(nrow(cag_abundance_df) == 0){
    return(ggplot() + geom_blank())
  } else {
    
    if(col != "None"){
      cag_abundance_df$col_variable <- cag_abundance_df %>% pull(col)
    }
    
    if(hue == "None" && group == "None"){
      p <- ggplot(
        data = cag_abundance_df,
        aes_string(
          x = x,
          y = "abundance"
        )
      )
    } else {
      if(hue != "None" && group != "None"){
        p <- ggplot(
          data = cag_abundance_df,
          aes_string(
            x = x,
            y = "abundance",
            color = hue,
            fill = hue,
            group = group
          )
        )
      } else {
        if(hue != "None"){
          p <- ggplot(
            data = cag_abundance_df,
            aes_string(
              x = x,
              y = "abundance",
              color = hue,
              fill = hue
            )
          )
        } else {
          p <- ggplot(
            data = cag_abundance_df,
            aes_string(
              x = x,
              y = "abundance",
              group = group
            )
          )
        }
      }
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
      p <- p + geom_bar(
        stat = "identity"
      )
    }
    if(logy){
      p <- p + scale_y_continuous(
        "Relative Abundance (log10)",
        trans='log10'
      )
    }
    return(
      p + theme_minimal(
        base_size = 20
      ) + theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      )
    )
  }
}

plot_cag_heatmap <- function(
  df,
  cag_list,
  color_by_primary,
  color_by_secondary,
  color_by_tertiary,
  show_column_names,
  checkbox_cluster
){
  if(nrow(df) == 0){
    return(ggplot() + geom_blank())
  }
  # If the user set the tertiary label but not the secondary
  # then just treat the tertiary as though it were the secondary
  if(color_by_secondary == "NONE" && color_by_tertiary != "NONE"){
    color_by_secondary <- color_by_tertiary
    color_by_tertiary <- "NONE"
  }

  # If the user set the secondary label but not the primary
  # then just treat the secondary as though it were the primary
  if(color_by_primary == "NONE" && color_by_secondary != "NONE"){
    color_by_primary <- color_by_secondary
    color_by_secondary <- "NONE"
  }
  
  # Sort the table depending on the user's input
  if(checkbox_cluster || color_by_primary == "NONE"){
    # Sort the table by hierarchical clustering of CAG abundances
    df <- df[hclust(
      dist(
        df %>% select(cag_list),
        method = "euclidean"
      ),
      method = "ward.D"
    )$order,]
  }else{
    # First sort by the tertiary key
    if(color_by_tertiary != "NONE"){
      df <- arrange(df, df[[color_by_tertiary]])
    }
    # Then sort by the secondary key
    if(color_by_secondary != "NONE"){
      df <- arrange(df, df[[color_by_secondary]])
    }
    # Then the primary key
    df <- arrange(df, df[[color_by_primary]])
  }
  
  # Duplicate the specimen label as a column name
  df$INDEX <- df$specimen
  df <- df %>% column_to_rownames(var = "INDEX")
  
  # Find the lowest non-zero value for the dataset,
  # which will be used to fill in the zero values
  # for plotting purposes (on the log scale)
  max_val <- df %>% 
    select(cag_list) %>%
    max
  min_val <- df %>%
    select(cag_list) %>%
    na_if(0) %>%
    mutate_if(is.numeric , replace_na, replace = max_val) %>%
    min

  # Format the table for plotting CAG abundance values
  abund_df <- df %>%
    select(cag_list) %>%
    na_if(0) %>%
    t %>%
    log10 %>%
    replace_na(
      min_val
    )
  
  # Annotation with the primary color_by metadata
  if(color_by_primary != "NONE"){
    g1 <- set_up_heatmap_annotation(df, color_by_primary)
  }

  # Annotation with the secondary color_by metadata
  if(color_by_secondary != "NONE"){
    g2 <- set_up_heatmap_annotation(df, color_by_secondary)
  }

  # Annotation with the tertiary color_by metadata
  if(color_by_tertiary != "NONE"){
    g3 <- set_up_heatmap_annotation(df, color_by_tertiary)
  }
  
  # Heatmap with CAG abundances
  g4 <- Heatmap(
    abund_df,
    name = "Prop. Abund. (log10)",
    rect_gp = gpar(col = "white", lwd = 0.5),
    column_title = "Specimens",
    row_title = "CAGs",
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = show_column_names
  )
  if(color_by_tertiary == "NONE"){
    if(color_by_secondary == "NONE"){
      if(color_by_primary == "NONE"){
        return(g4)
      } else {
        return(g1 %v% g4)
      }
    } else {
      return(g1 %v% g2 %v% g4)
    }
  } else {
    return(g1 %v% g2 %v% g3 %v% g4)
  }
}

set_up_heatmap_annotation <- function(df, color_by){
  # Use a dedicated function to set the colors of the annotation rows
  cmap <- color_by_column(
    df[[color_by]]
  )
  # Set up the heatmap
  g <- Heatmap(
    df %>%
      select(
        color_by
      ) %>%
      t,
    name = color_by,
    rect_gp = gpar(col = "white", lwd = 0.5),
    column_title = "Specimens",
    row_title = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  return(g)
}

color_by_column <- function(col_vals){
  if(all(is.character(col_vals)) || length(unique(col_vals)) < 5){
    # Make a list of colors to user
    n <- length(unique(col_vals))
    color_list <- rep(
      brewer.pal(8, "Paired"), 
      n
    )[1:n]
    cmap <- structure(
      color_list,
      names = unique(as.character(col_vals))
    )
  } else {
    col_vals <- as.numeric(col_vals)
    cmap <- colorRamp2(
      c(
        min(col_vals),
        median(col_vals),
        max(col_vals)
      ),
      c("blue", "white", "red")
    )
  }
  return(cmap)
}
