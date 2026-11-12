#' @title CTS engine for 2 sources
#'
#' @description Solves determined systems for 1 tracer and 2 sources to assess 
#' solution dispersion and discriminant capacity.
#'
#' @param source Data frame of sediment sources.
#' @param mixture Data frame of mixture samples.
#' @param iter Number of Monte Carlo iterations.
#' @param seed Integer for reproducibility.
#'
#' @return A data frame with one row per tracer, sorted by \code{max_dw}
#'   (most discriminant first), with columns: \code{id} (tracer name),
#'   \code{w1}, \code{w2} (source contributions from mean values),
#'   \code{dw1}, \code{dw2} (dispersion of each contribution, half the
#'   16-84\% quantile range), \code{cons} (proportion of physically
#'   feasible solutions) and \code{max_dw} (maximum dispersion).
#'
#' @export
#' 
CTS_2_sources <- function(source, mixture, iter = 1000, seed = 123456) {
  
  set.seed(seed)
  # Convert to matrix excluding the first column (ID)
  source_mat <- data.matrix(source[-1])
  mixture_mat <- data.matrix(mixture[-1])
  
  ntracers <- (ncol(source_mat) - 1) / 2
  tracer_names <- colnames(source_mat)[1:ntracers]
  n_col_idx <- ncol(source_mat) # Column for 'n'
  
  df_results <- data.frame()
  pb <- txtProgressBar(min = 0, max = ntracers, style = 3, width = 50, char = "=")
  
  for (i in 1:ntracers) {
    w1 <- numeric(iter)
    w2 <- numeric(iter)
    sc <- 0 # Consistency counter
    
    for (l in 1:iter) {
      if (l == 1) {
        # Mean values for the first iteration
        s11 <- source_mat[1, i]
        s21 <- source_mat[2, i]
        m1  <- mixture_mat[i]
      } else {
        # Variability using t-distribution and the 'n' column
        s11 <- source_mat[1, i] + source_mat[1, i + ntracers] * rt(1, source_mat[1, n_col_idx]) / sqrt(source_mat[1, n_col_idx])
        s21 <- source_mat[2, i] + source_mat[2, i + ntracers] * rt(1, source_mat[2, n_col_idx]) / sqrt(source_mat[2, n_col_idx])
        m1  <- mixture_mat[i]
      }
      
      # System determinant (s21 - s11)
      det <- s21 - s11
      
      if (det != 0) {
        x1 <- (s21 - m1) / det
        x2 <- (m1 - s11) / det
        
        # Check physical feasibility (0-1)
        if (min(x1, x2) >= 0 && max(x1, x2) <= 1) sc <- sc + 1
        
        w1[l] <- x1
        w2[l] <- x2
      }
    }
    
    # Calculate dispersion (Dw) using 16-84% quantiles
    # This represents the uncertainty of the source apportionment
    dw1 <- (quantile(w1, 0.84, na.rm = TRUE) - quantile(w1, 0.16, na.rm = TRUE)) / 2
    dw2 <- (quantile(w2, 0.84, na.rm = TRUE) - quantile(w2, 0.16, na.rm = TRUE)) / 2
    
    res <- data.frame(
      id = tracer_names[i],
      w1 = w1[1],
      w2 = w2[1],
      dw1 = dw1,
      dw2 = dw2,
      cons = sc / iter,
      max_dw = max(dw1, dw2)
    )
    
    df_results <- rbind(df_results, res)
    setTxtProgressBar(pb, i)
  }
  
  # Sort by lowest dispersion (most discriminant)
  df_results <- df_results[order(df_results$max_dw), ]
  rownames(df_results) <- NULL
  
  close(pb)
  return(df_results)
}

# --- Compatibility Aliases ---

#' @rdname CTS_2_sources
#' @export
singles <- CTS_2_sources

#' @rdname CTS_2_sources
#' @export
CTS_seeds_singles <- CTS_2_sources