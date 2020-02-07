# Source the Python / Pandas code needed to read HDF5
use_python("/Users/sminot/miniconda3/bin/python3", required = TRUE)
pandas <- import("pandas")

# Walk through the specified directory and find all files ending with HDF5
find_datasets <- function(directory){
    # Find all of the files in the folder
    files <- list.files(directory)
    
    # Subset to those which end with ".hdf5"
    files <- files[sapply(files, function(n){endsWith(n, ".hdf5")})]
    return(sub(".hdf5", "", files))
}

# Function to determine whether the input dataset has any corncob results
read_hdf_has_corncob <- function(dataset_prefix, data_folder){
  return(
    tryCatch(
      {
        df <- pandas$read_hdf(
          file.path(
            data_folder, 
            paste(dataset_prefix, ".hdf5", sep="")
          ),
          "/stats/cag/corncob"
        )
        return(TRUE)
      },
      error = function(c) FALSE
    )
  )
}
  
# Function to read in the results of statistical analysis from the HDF store
read_hdf_corncob_results <- function(dataset_prefix, data_folder){
  df <- pandas$read_hdf(
    file.path(
      data_folder, 
      paste(dataset_prefix, ".hdf5", sep="")
    ),
    "/stats/cag/corncob"
  ) %>% filter(grepl("mu.", parameter))
  
  df$parameter <- sapply(df$parameter, function(n){str_replace(n, "mu.", "")})
  
  df <- df %>% pivot_wider(names_from = type, values_from = value)
  
  return(df)
}

# Function to read the abundance for all CAGs
read_cag_abundances <- function(dataset_prefix, data_folder){

  return(
    pandas$read_hdf(
      file.path(
        data_folder, 
        paste(dataset_prefix, ".hdf5", sep="")
      ),
      "/abund/cag/wide"
    ) %>% 
      as_tibble %>% 
      column_to_rownames(var = "CAG")
  )
  
}


# Function to read the annotations for every unique gene
read_gene_annotations <- function(dataset_prefix, data_folder){
  
  return(
    pandas$read_hdf(
      file.path(
        data_folder, 
        paste(dataset_prefix, ".hdf5", sep="")
      ),
      "/annot/gene/all"
    ) %>% 
      as_tibble
  )
  
}

# Function to read the sample metadata (manifest) for this study
read_hdf_manifest <- function(dataset_prefix, data_folder){
  return(
    pandas$read_hdf(
      file.path(
        data_folder, 
        paste(dataset_prefix, ".hdf5", sep="")
      ),
      "/manifest"
    ) %>% 
      as_tibble  %>% 
      select(-R1) %>% 
      select(-R2) %>% 
      select(-I1) %>% 
      select(-I2)
  )
}

# Function to read the readcounts table for this study
read_hdf_readcounts <- function(dataset_prefix, data_folder){
  df <- pandas$read_hdf(
    file.path(
      data_folder, 
      paste(dataset_prefix, ".hdf5", sep="")
    ),
    "/summary/readcount"
  ) %>% 
    as_tibble
  return(df %>%
     add_column(
       prop_aligned = df$aligned_reads / df$n_reads
     )
  )
}

# Function to assemble a summary of all CAGs
make_cag_summary <- function(
  gene_annotations_df,
  cag_abund_df,
  corncob_results_df
){
  summary_df <- data.frame(
    table(gene_annotations_df$CAG)
  ) %>% 
    as_tibble %>%
    rename(
      "CAG" = Var1,
      "Number of Genes" = Freq
    ) %>% 
    mutate_if(is.factor, ~ as.numeric(as.character(.x)))
  
  summary_df["Prevalence"] <- apply(
    cag_abund_df,
    1,
    function(v){mean(v > 0)}
  )
  summary_df["Mean Abund"] <- apply(cag_abund_df, 1, mean)
  summary_df["Max Abund"] <- apply(cag_abund_df, 1, max)
  summary_df["Min Abund"] <- apply(cag_abund_df, 1, min)
  
  if("CAG" %in% colnames(corncob_results_df)){
    summary_df <- summary_df %>%
      full_join(corncob_results_df, by = "CAG") %>%
      add_column(
        neg_log_pvalue = sapply(
          pull(corncob_results_df, p_value),
          function(v){-log10(v)}
        )
      )
  }
  return(summary_df)
}
