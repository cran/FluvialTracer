#' @title Consistent Tracer Selection (CTS) for multiple sources
#'
#' @description Identifies the most discriminant minimal tracer combinations (n-1) 
#' by solving determined systems and assessing solution dispersion. Lower dispersion 
#' indicates higher discriminant capacity.
#'
#' @param data Data frame containing sediment source and mixture data.
#' @param iter Number of Monte Carlo iterations for variability analysis.
#' @param seed Integer for reproducibility.
#'
#' @return A data frame with one row per tracer combination, sorted by
#'   \code{max_dw} (lowest dispersion first), with columns: \code{id}
#'   (tracer combination), \code{w1}, \code{w2}, ... (source apportionments),
#'   \code{dw1}, \code{dw2}, ... (solution dispersion), \code{\% in sources range}
#'   (percentage of physically feasible solutions) and \code{max_dw}
#'   (maximum dispersion).
#'
#' @references
#' Latorre, B., Lizaga, I., Gaspar, L., & Navas, A. (2021). A novel method for analysing
#' consistency and unravelling multiple solutions in sediment fingerprinting.
#'
#' @export
#' 
CTS_n_sources <- function(data, iter = 1000, seed = 123456)
{
  source <- inputSource(data)
  mixture <- inputSample(data)
  n <- nrow(source)
  
  colnames(source) <- gsub("(?i)^mean_", "", colnames(source))
  colnames(mixture) <- gsub("(?i)^mean_", "", colnames(mixture))
  
  if(n == 2) {
    seeds <- CTS_2_sources(source, mixture, iter, seed)
  }
  else if(n == 3) {
    seeds <- CTS_3_sources(source, mixture, iter, seed)
  }
  else if(n == 4) {
    seeds <- CTS_4_sources(source, mixture, iter, seed)
  }
  else if(n == 5) {
    seeds <- CTS_5_sources(source, mixture, iter, seed)
  }
  else {
    stop(paste0("The CTS method currently supports datasets with 2 to 5 sediment sources."))
  }
  
  colnames(seeds)[which(colnames(seeds) == "cons")] <- "% in sources range"
  seeds$`% in sources range` <- seeds$`% in sources range` * 100.0
  
  return(seeds)
}

# --- Compatibility Aliases ---
#' @rdname CTS_n_sources
#' @export
CTS_seeds <- CTS_n_sources

    