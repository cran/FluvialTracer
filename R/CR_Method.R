#' @title Consensus Ranking (CR) method for tracer selection
#'
#' @description Computes the Consensus Ranking (CR) method, an ensemble technique to identify non-conservative and dissenting tracers. The method combines predictions from single-tracer models and is based on a scoring function derived from a series of random "debates" between tracers, utilizing advanced thresholding logic for extreme datasets.
#'
#' @param data Data frame. Contains source and mixture data in the standard FluvialTracer format.
#' @param debates Integer. Target number of random debates per tracer. Default 1000.
#' @param completion_quantile Numeric. Value between 0-1. Quantile threshold for debate completion. For example, 0.5 stops the simulation once half the tracers hit the debates target. Default 1.
#' @param apportionment_threshold Numeric. Tolerance boundary for valid debate apportionments (e.g., defining valid as -threshold to 1 + threshold). Default 0.
#' @param seed Integer. Random seed for reproducibility. Default 123456.
#'
#' @return A data frame containing the CR score, debates won, and debates lost per tracer.
#'
#' @references
#' Lizaga, I., Latorre, B., Gaspar, L., & Navas, A. (2020). Consensus ranking as a method to identify non-conservative and dissenting tracers in fingerprinting studies. *Science of The Total Environment*, *720*, 137537. https://doi.org/10.1016/j.scitotenv.2020.137537
#' 
#' @export
#' 
CR_Method <- function(data, debates = 1000, completion_quantile = 1, apportionment_threshold = 0, seed = 123456) {
  
  # Data extraction using internal library functions
  source_df <- inputSource(data)
  mixture_df <- inputSample(data)
  
  source <- data.matrix(source_df[-1])
  mixture <- data.matrix(mixture_df[-1])
  
  set.seed(seed)
  
  # Dimensions and parameter setup
  ntracer <- (ncol(source) - 1) / 2
  nsource <- nrow(source)
  
  # Technical normalization
  for (col in 1:ntracer) {
    mx <- max(source[, col] + source[, ntracer + col])
    mn <- min(source[, col] - source[, ntracer + col])
    if (mx != mn) {
      source[, col] <- (source[, col] - mn) / (mx - mn)
      source[, ntracer + col] <- source[, ntracer + col] / (mx - mn)
      mixture[, col] <- (mixture[, col] - mn) / (mx - mn)
    }
  }
  
  tracer_names <- colnames(source)[1:ntracer]
  tracer_names <- gsub("^Mean_", "", tracer_names)
  
  # Counters for participations (cr1) and losses/dissent (cr2)
  cr1 <- rep(0, ntracer)
  cr2 <- rep(0, ntracer)
  
  pb <- txtProgressBar(min = 0, max = debates, style = 3, width = 50, char = "=")
  
  # Main simulation loop using improved quantile logic
  while (stats::quantile(cr1, completion_quantile) < debates) {
    var_idx <- sample(1:ntracer, nsource + 1)
    
    if (min(cr1[var_idx]) < debates) {
      gof_scores <- c()
      max_w <- c()
      min_w <- c()
      
      # Variability matrix using Student's t distribution
      rnd_matrix <- matrix(0, nrow = nsource, ncol = nsource + 1)
      n_col_idx <- ncol(source)
      for (i in 1:(nsource + 1)) {
        for (k in 1:nsource) {
          rnd_matrix[k, i] <- stats::rt(1, source[k, n_col_idx]) / sqrt(source[k, n_col_idx])
        }
      }
      
      # Debate sub-loop (Leave-one-out)
      for (i in 1:(nsource + 1)) {
        active_debate_idx <- setdiff(1:(nsource + 1), i)
        
        y_vec <- c()
        x_mat <- matrix(NA, nrow = nsource, ncol = nsource - 1)
        csource_sim <- matrix(NA, nrow = nsource, ncol = nsource)
        ls_val <- c()
        
        for (jj in 1:nsource) {
          k_idx <- active_debate_idx[jj]
          l_idx <- var_idx[k_idx]
          
          csource_sim[nsource, jj] <- source[nsource, l_idx] + source[nsource, ntracer + l_idx] * rnd_matrix[nsource, k_idx]
          ls_val <- c(ls_val, csource_sim[nsource, jj])
          
          y_vec <- c(y_vec, mixture[l_idx] - ls_val[jj])
          for (kk in 1:(nsource - 1)) {
            csource_sim[kk, jj] <- source[kk, l_idx] + source[kk, ntracer + l_idx] * rnd_matrix[kk, k_idx]
            x_mat[jj, kk] <- csource_sim[kk, jj] - ls_val[jj]
          }
        }
        
        # Solve via linear least squares
        model <- lm.fit(x_mat, y_vec)
        coeffs <- as.vector(stats::coef(model))
        weights <- c(coeffs, 1 - sum(coeffs))
        
        # Goodness of Fit (GOF) calculation
        err <- 0
        for (jj in 1:nsource) {
          l_orig_idx <- var_idx[active_debate_idx[jj]]
          pred_val <- sum(csource_sim[, jj] * weights)
          err <- err + (pred_val - mixture[l_orig_idx])^2
        }
        
        gof_scores <- c(gof_scores, err)
        max_w <- c(max_w, max(weights))
        min_w <- c(min_w, min(weights))
      }
      
      # Consensus vs. Dissent logic with thresholds
      sorted_gof <- sort(gof_scores)
      if (sorted_gof[nsource + 1] - sorted_gof[nsource] > 10.0 * (sorted_gof[nsource] - sorted_gof[1])) {
        dissenter_idx <- which.max(gof_scores)
        if (cr1[var_idx[dissenter_idx]] < debates) {
          cr1[var_idx[dissenter_idx]] <- cr1[var_idx[dissenter_idx]] + 1
          cr2[var_idx[dissenter_idx]] <- cr2[var_idx[dissenter_idx]] + 1
        }
      } else {
        best_fit_idx <- which.min(gof_scores)
        if (max_w[best_fit_idx] <= (1.0 + apportionment_threshold) && 
            min_w[best_fit_idx] >= (0.0 - apportionment_threshold)) {
          for (idx in 1:(nsource + 1)) {
            if (cr1[var_idx[idx]] < debates) {
              cr1[var_idx[idx]] <- cr1[var_idx[idx]] + 1
              if (idx == best_fit_idx) cr2[var_idx[idx]] <- cr2[var_idx[idx]] + 1
            }
          }
        }
      }
    }
    setTxtProgressBar(pb, stats::quantile(cr1, completion_quantile))
  }
  
  # Full Pro Output
  CR_score <- 100 - ((cr2 * 100) / cr1)
  cr_results <- data.frame(tracer = tracer_names, 
                           CR_score = round(CR_score, 2), 
                           debates = cr1, 
                           debates_lost = cr2)
  
  cr_results <- cr_results[order(-cr_results$CR_score), ]
  rownames(cr_results) <- NULL
  
  close(pb)
  return(cr_results)
}

# --- Backwards Compatibility Aliases ---
#' @rdname CR_Method
#' @export
CR <- CR_Method

#' @rdname CR_Method
#' @export
cr_ns <- CR_Method