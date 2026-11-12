#' Checking the difference between virtual and artificial mixture (for data sets containing multiple mixtures)
#'
#' A function that calculates and visualizes the percentage deviation between 
#' an experimental artificial mixture and a theoretical virtual mixture calculated from user-defined source proportions.
#' Source proportions should be defined as described in Option A: manual list or Option B: CSV file 
#' 
#'
#' @param data Data frame. Contains source and mixture data in the standard format.
#' @param proportions Data.frame. With one row per mixture and one column per source.
#'   Can be created in two ways:
#'   
#'   Option A - manually:
#'   \code{proportions = data.frame(rbind(c(0.5, 0.25, 0.25, 0), c(0.3, 0.3, 0.2, 0.2)))}
#'   
#'   Option B - from a CSV file (one row per mixture, one column per source):
#'   \code{proportions = read.csv()} 
#'   
#' @param sample_col Character. Defines which column of the data contains the sample names.
#' @param source_col Character. Defines which column of the data contains the source names. Indicate the column number in which the sources are defined independently from the sample number. 
#' 
#' @return No return value, called for side effects (draws one plot per mixture
#'   showing the deviation between the virtual and artificial mixture).
#' 
#' @importFrom utils read.csv capture.output
#' @importFrom dplyr slice bind_rows
#' 
#' @export


checking_diff_mult_mix <- function(data, 
                                   proportions = props_table,
                                   sample_col = 1,
                                   source_col = 2)
{
  
  if (!is.character(source_col)) source_col <- colnames(data)[as.numeric(source_col)]
  if (!is.character(sample_col)) sample_col <- colnames(data)[as.numeric(sample_col)]
  
  # Resolve proportions
  if (is.character(proportions) && file.exists(proportions)) {
    props_table <- read.csv(proportions)
  } else {
    props_table <- as.data.frame(proportions)
  }
  
  # Detect sources and mixtures
  all_unique       <- unique(data[[source_col]])
  all_unique_clean <- gsub("\\d+$", "", all_unique)
  mix_base         <- tail(unique(all_unique_clean), 1)
  mix_ids          <- all_unique[grepl(paste0("^", mix_base, "\\d*$"), all_unique, ignore.case = TRUE)]
  source_names     <- all_unique[!all_unique %in% mix_ids]
  n_mixtures       <- sum(data[[source_col]] %in% mix_ids)
  
  # Loop over each mixture
  for (i in 1:n_mixtures) {
    
    message("--- Processing mixture: ", i, " of ", n_mixtures, " ---")
    
    # Take proportions from row i
    props_i <- as.numeric(props_table[i, sapply(props_table, is.numeric)])
    
    # Subset with sources + row i of the mixture
    mix_row <- data %>%
      filter(.data[[source_col]] %in% mix_ids) %>%
      slice(i)
    
    data_sub <- bind_rows(
      data %>% filter(.data[[source_col]] %in% source_names),
      mix_row
    )
    
    # Build subtitle with iteration info
    props_label <- paste(source_names, props_i, sep = " = ", collapse = " | ")
    mix_name    <- mix_row[[sample_col]]
    
    subtitule <- paste0(
      "Mixture: ", i, " of ", n_mixtures, "\n",
      "Proportions: ", props_label
    )
    
    # Call the original function (it draws the plot itself)
    p <- checking_diff_single_mix(data_sub, 
                                  proportions = props_i,
                                  sample_col  = sample_col,
                                  source_col  = source_col,
                                  subtitle_extra = subtitule)
  }
}