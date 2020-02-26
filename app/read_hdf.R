# Source the Python / Pandas code needed to read HDF5
use_python("/Users/sminot/miniconda3/bin/python3", required = TRUE)
pandas <- import("pandas")

# Function to the summary of the entire expeirment
read_hdf_summary <- function(dataset_prefix, data_folder){
  df <- pandas$read_hdf(
    file.path(
      data_folder, 
      paste(dataset_prefix, ".hdf5", sep="")
    ),
    "/summary/experiment"
  )
  return(df)
}

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
          "/stats/cag/corncob_wide"
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
    "/stats/cag/corncob_wide"
  )
  return(df)
}

# Function to read the abundance for all CAGs
read_hdf_cag_abundances <- function(dataset_prefix, data_folder){

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


# Function to read the annotations for every unique gene within a CAG
read_hdf_cag_details <- function(dataset_prefix, data_folder, cag_id){
  return(
    pandas$read_hdf(
      file.path(
        data_folder, 
        paste(dataset_prefix, ".hdf5", sep="")
      ),
      "/annot/gene/all",
      where = paste("CAG ==", as.character(cag_id))
    ) %>% 
      as_tibble
  )
  
}

# Function to read the sample metadata (manifest) for this study
read_hdf_manifest <- function(dataset_prefix, data_folder){
  df <- pandas$read_hdf(
      file.path(
        data_folder, 
        paste(dataset_prefix, ".hdf5", sep="")
      ),
      "/manifest"
    ) %>% 
      as_tibble  %>% 
      select(-R1) %>% 
      select(-R2)
  if("I1" %in% colnames(df)){
    df <- df %>% 
      select(-I1)
  }
  if("I2" %in% colnames(df)){
    df <- df %>% 
      select(-I2)
  }
  # Convert columns with very few unique values to factors
  if(nrow(df) > 5){
    for(n in colnames(df)){
      if(length(unique(df[[n]])) <= 5){
        if(is.numeric(df[[n]])){
          df[[n]] <- as.factor(df[[n]])
          fct_inseq(df[[n]])
        } else {
          df[[n]] <- as.factor(df[[n]])
        }
        
      }
    }
  }
  return(df)
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

# Read in the CAG summary table
read_hdf_cag_summary <- function(dataset_prefix, data_folder){
  summary_df <- pandas$read_hdf(
    file.path(
      data_folder, 
      paste(dataset_prefix, ".hdf5", sep="")
    ),
    "/annot/cag/all"
  ) %>% 
    as_tibble
  # Add the log10 size
  summary_df$size_log10 <- sapply(
    summary_df$size,
    log10
  )
  summary_df$mean_abundance_log10 <- sapply(
    summary_df$mean_abundance,
    log10
  )
  return(summary_df)
}

# Read in the PCA results
read_hdf_pca <- function(dataset_prefix, data_folder){
  return(pandas$read_hdf(
    file.path(
      data_folder, 
      paste(dataset_prefix, ".hdf5", sep="")
    ),
    "/ordination/pca"
  ) %>% 
    as_tibble)
}

# Read in the t-SNE results
read_hdf_tsne <- function(dataset_prefix, data_folder){
  return(pandas$read_hdf(
    file.path(
      data_folder, 
      paste(dataset_prefix, ".hdf5", sep="")
    ),
    "/ordination/tsne"
  ) %>% 
    as_tibble)
}
