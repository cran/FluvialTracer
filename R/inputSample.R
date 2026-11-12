#' Extract mixture samples from a fingerprinting dataset
#'
#' Identifies and returns the rows corresponding to the sediment mixture(s)
#' in a combined source-mixture data frame. The mixture is assumed to be the
#' last unique group in the second column of \code{data}. Both raw and
#' averaged (\code{Mean_}/\code{SD_}) formats are supported.
#'
#' @param data Data frame in the standard FluvialTracer format, where the second
#'   column encodes the group identity (sources and mixture).
#'
#' @return A data frame containing only the mixture rows, with column names
#'   standardised to the FluvialTracer convention.
#'
#' @seealso \code{\link{inputSource}} for extracting source rows.
#' 
#'
#' @export
inputSample <- function(data) {
  
  # --- Internal Helper: Detect data format ---
  is_averaged_internal <- function(d) {
    # Check for the existence of the 'n' column (case-insensitive)
    if (sum(grepl("^n$", tolower(colnames(d)))) == 1) {
      # Check if column structure matches (ID, Group, [tracers...], n)
      ntracer <- (ncol(d) - 3) / 2
      # Verify if even number of columns for mean/sd exists
      if (ntracer > 0 && ntracer %% 1 == 0) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
  
  # --- Case 1: Averaged Data ---
  if (is_averaged_internal(data)) {
    # Calculate number of tracers based on (Mean + SD + n) structure
    ntracer <- (ncol(data) - 3) / 2
    # Extract only the mixture row (last row) and the mean columns (1 to ntracer+2)
    # Re-calling inputSample on the subset will skip this block and go to Case 2 logic
    mixtures <- inputSample(data[, c(1:(ntracer + 2))])
    return(mixtures)
  }
  
  # --- Case 2: Raw Data ---
  # Ensure groups are factors and maintain original order
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  # Extract the last group (assumed to be the mixture)
  mixtures <- data[data[, 2] == levels(data[, 2])[nlevels(data[, 2])], ]
  
  # Standardize output format
  rownames(mixtures) <- NULL
  # Remove ID and Group columns for processing
  samples <- as.data.frame(mixtures[, -(1:2)])
  
  # Clean names: remove existing 'Mean_' or 'mean_' prefixes and re-apply 'mean_'
  colnames(samples) <- gsub("(?i)^Mean_", "", colnames(samples))
  colnames(samples) <- paste0('Mean_', colnames(samples))
  
  # Re-attach the ID/Group identifiers
  id <- as.vector(tail(data[, 2], n = nrow(samples)))
  samples <- cbind(id, samples)
  
  return(samples)
}


# --- Backwards Compatibility Aliases ---
#' @rdname inputSample
#' @export
inputMixture <- inputSample