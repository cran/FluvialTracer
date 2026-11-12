#' @title CTS engine for 3 sources
#'
#' @description Solves determined systems for 2 tracers and 3 sources to assess 
#' solution dispersion and discriminant capacity.
#'
#' @param source Data frame of sediment sources.
#' @param mixture Data frame of mixture samples.
#' @param iter Number of Monte Carlo iterations.
#' @param seed Integer for reproducibility.
#'
#' @return A data frame with one row per tracer pair, sorted by \code{max_dw}
#'   (most discriminant first), with columns: \code{id} (tracer pair),
#'   \code{w1}, \code{w2}, \code{w3} (source contributions from mean values),
#'   \code{dw1}, \code{dw2}, \code{dw3} (dispersion of each contribution, half
#'   the 16-84\% quantile range), \code{cons} (proportion of physically
#'   feasible solutions) and \code{max_dw} (maximum dispersion).
#'
#' @export
#' 
CTS_3_sources <- function(source, mixture, iter = 1000, seed = 123456) {
  
  set.seed(seed)
  # Convert to matrix excluding ID
  source_mat <- data.matrix(source[-1])
  mixture_mat <- data.matrix(mixture[-1])
  
  ntracers <- (ncol(source_mat) - 1) / 2
  tracer_names <- colnames(source_mat)[1:ntracers]
  n_col_idx <- ncol(source_mat)
  
  df_results <- data.frame()
  
  # Combinations of 2 tracers for 3 sources
  comb_matrix <- utils::combn(ntracers, 2)
  n_comb <- ncol(comb_matrix)
  
  pb <- txtProgressBar(min = 0, max = n_comb, style = 3, width = 50, char = "=")
  
  for (i in 1:n_comb) {
    t1 <- comb_matrix[1, i]
    t2 <- comb_matrix[2, i]
    
    w1 <- numeric(iter); w2 <- numeric(iter); w3 <- numeric(iter)
    sc <- 0
    
    for (l in 1:iter) {
      if (l == 1) {
        s11 <- source_mat[1, t1]; s12 <- source_mat[1, t2]
        s21 <- source_mat[2, t1]; s22 <- source_mat[2, t2]
        s31 <- source_mat[3, t1]; s32 <- source_mat[3, t2]
        m1  <- mixture_mat[t1];   m2  <- mixture_mat[t2]
      } else {
        # Variability using t-distribution and 'n' column
        deg_f <- source_mat[, n_col_idx]
        s11 <- source_mat[1, t1] + source_mat[1, t1 + ntracers] * rt(1, deg_f[1]) / sqrt(deg_f[1])
        s12 <- source_mat[1, t2] + source_mat[1, t2 + ntracers] * rt(1, deg_f[1]) / sqrt(deg_f[1])
        s21 <- source_mat[2, t1] + source_mat[2, t1 + ntracers] * rt(1, deg_f[2]) / sqrt(deg_f[2])
        s22 <- source_mat[2, t2] + source_mat[2, t2 + ntracers] * rt(1, deg_f[2]) / sqrt(deg_f[2])
        s31 <- source_mat[3, t1] + source_mat[3, t1 + ntracers] * rt(1, deg_f[3]) / sqrt(deg_f[3])
        s32 <- source_mat[3, t2] + source_mat[3, t2 + ntracers] * rt(1, deg_f[3]) / sqrt(deg_f[3])
        m1  <- mixture_mat[t1];   m2  <- mixture_mat[t2]
      }
      
      # Optimized Determinant (Maxima script)
      det <- -s22 * s31 + s12 * s31 + s21 * s32 - s11 * s32 + s11 * s22 - s12 * s21
      
      if (det != 0) {
        x1 <- (s21 * s32 - s22 * s31 + m1 * (s22 - s32) + m2 * (s31 - s21)) / det
        x2 <- (s12 * s31 - s11 * s32 + m1 * (s32 - s12) + m2 * (s11 - s31)) / det
        x3 <- (s11 * s22 - s12 * s21 + m1 * (s12 - s22) + m2 * (s21 - s11)) / det
        
        if (min(x1, x2, x3) >= 0 && max(x1, x2, x3) <= 1) sc <- sc + 1
        
        w1[l] <- x1; w2[l] <- x2; w3[l] <- x3
      }
    }
    
    # Uncertainty (16-84% quantiles)
    dw1 <- (quantile(w1, 0.84, na.rm = TRUE) - quantile(w1, 0.16, na.rm = TRUE)) / 2
    dw2 <- (quantile(w2, 0.84, na.rm = TRUE) - quantile(w2, 0.16, na.rm = TRUE)) / 2
    dw3 <- (quantile(w3, 0.84, na.rm = TRUE) - quantile(w3, 0.16, na.rm = TRUE)) / 2
    
    res <- data.frame(
      id = paste(tracer_names[t1], tracer_names[t2], sep = " "),
      w1 = w1[1], w2 = w2[1], w3 = w3[1],
      dw1 = dw1, dw2 = dw2, dw3 = dw3,
      cons = sc / iter,
      max_dw = max(dw1, dw2, dw3)
    )
    
    df_results <- rbind(df_results, res)
    setTxtProgressBar(pb, i)
  }
  
  df_results <- df_results[order(df_results$max_dw), ]
  rownames(df_results) <- NULL
  
  close(pb)
  return(df_results)
}

# --- Aliases ---
#' @rdname CTS_3_sources
#' @export
pairs <- CTS_3_sources

#' @rdname CTS_3_sources
#' @export
CTS_seeds_pairs <- CTS_3_sources