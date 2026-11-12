#' Unmix sediment or marine particle mixtures
#'
#' Estimates the fractional contribution of each potential source to one or
#' more mixtures. Users can seamlessly switch between a lightning-fast Frequentist 
#' least-squares approach (with optional Linear Variability Propagation), a 
#' standard Bayesian MCMC approach, or a Machine Learning approach.
#'
#' Both raw (individual observations) and averaged (\code{Mean_}/\code{SD_}
#' columns) input formats are accepted. When multiple mixture samples are
#' present, the model processes each independently.
#'
#' @param data Data frame in the standard FluvialTracer format.
#' @param model Character. \code{"frequentist"} (default), \code{"bayesian"} or \code{"machine_learning"}.
#' @param iter Integer. Number of iterations for frequentist Monte Carlo/LVP.
#'   Default \code{1000}.
#' @param variability Character. Uncertainty metric: \code{"SEM"} (Standard
#'   Error of the Mean, default) or \code{"SD"} (Standard Deviation).
#' @param lvp Logical. \code{TRUE} (default) uses Linear Variability
#'   Propagation; \code{FALSE} uses classical Monte Carlo resampling.
#' @param seed Integer. Random seed for reproducibility. Default \code{123456}.
#' @param ml_type Character. Only used if \code{model = "machine_learning"}. 
#'   Options are \code{"RandomForest"} (default) or \code{"XGBoost"}.
#' @param ... Additional arguments passed to the Bayesian or ML engines (e.g., \code{run = "test"}, \code{n_trees = 500}, \code{stretch_bounds = 1.2}).
#'
#' @return For the \code{"frequentist"} model, a data frame with one row per iteration 
#'   and mixture. For the \code{"bayesian"} model, a named list containing the 
#'   formatted apportionment data frame, the MixSIAR summary statistics, and the 
#'   raw JAGS model object. For the \code{"machine_learning"} model, a named list
#'   containing the apportionment predictions and a statistical summary (\code{rf_summary} or \code{xgb_summary}).
#'
#' @references
#' Latorre, B., Lizaga, I., Gaspar, L., & Navas, A. (2025). Evaluating the
#' impact of high source variability and extreme contributing sources on
#' sediment fingerprinting models. \emph{Water Resources Management}.
#' \doi{10.1007/s11269-025-04169-8}
#'
#' @export
#' 
unmixing <- function(data, 
                     model = "frequentist", 
                     iter = 1000L, 
                     variability = "SEM",
                     lvp = TRUE, 
                     seed = 123456L, 
                     ml_type = "RandomForest",
                     ...) {
  
  # ── Internal helper: detect averaged format (Mean_/SD_ or mean_/sd_) ──────
  is_averaged_internal <- function(d) {
    if (sum(grepl("^n$", tolower(colnames(d)))) != 1) return(FALSE)
    if ((ncol(d) - 3) %% 2 != 0) return(FALSE)
    ntracer <- (ncol(d) - 3) / 2
    nms <- colnames(d)[-c(1:2)][1:(2 * ntracer)]
    v1  <- nms[grepl("^[Mm]ean_", nms)]
    v2  <- nms[grepl("^[Ss][Dd]_", nms)]
    length(v1) == ntracer && length(v2) == ntracer &&
      all(gsub("^[Mm]ean_", "", v1) == gsub("^[Ss][Dd]_", "", v2))
  }
  
  # ── Multiple mixtures: process each individually ──────────────────────────
  mixture_n <- nrow(inputSample(data))
  if (mixture_n > 1) {
    if (tolower(model) == "bayesian") {
      message("----------------------------------------------------------------")
      message("WARNING: Multiple mixture samples detected for Bayesian unmixing.")
      message("MixSIAR MCMC chains will run sequentially for each mixture.")
      message("Depending on your 'run' length, this may take a considerable amount of time.")
      message("----------------------------------------------------------------")
    } else if (tolower(model) == "machine_learning") {
      message("Multiple mixture samples detected. Training and predicting sequentially.")
    } else {
      warning("Multiple mixture samples detected. Each sample will be processed individually.")
    }
    
    results <- NULL
    for (i in seq_len(mixture_n)) {
      exclude      <- setdiff(seq_len(mixture_n), i)
      exclude_rows <- nrow(data) - mixture_n + exclude
      data_one     <- data[-exclude_rows, ]
      
      # Recursive call using the function name (passing the new ml_type parameter)
      res_i        <- unmixing(data_one, model, iter, variability, lvp, seed, ml_type, ...)
      
      if (tolower(model) %in% c("bayesian", "machine_learning")) {
        # For Bayesian and ML, we must carefully append lists
        if (is.null(results)) {
          results <- res_i
        } else {
          results$apportionment <- rbind(results$apportionment, res_i$apportionment)
          
          if (tolower(model) == "bayesian") {
            results$bayesian_summary <- rbind(results$bayesian_summary, res_i$bayesian_summary)
            # We only keep the last jags_model in the list to save RAM
            results$jags_model <- res_i$jags_model 
          } else if (tolower(model) == "machine_learning") {
            # Check which summary to bind based on ml_type
            if (tolower(ml_type) %in% c("xgboost", "xgb")) {
              results$xgb_summary <- rbind(results$xgb_summary, res_i$xgb_summary)
            } else {
              results$rf_summary <- rbind(results$rf_summary, res_i$rf_summary)
            }
          }
        }
      } else {
        # For Frequentist, just rbind the dataframes
        results <- if (is.null(results)) res_i else rbind(results, res_i)
      }
    }
    return(results)
  }
  
  # ── Single mixture: run unmixing ──────────────────────────────────────────
  if (tolower(model) == "frequentist") {
    
    means <- is_averaged_internal(data)
    
    if (!lvp) {
      results <- frequentist_unconstrained(data, variability = variability,
                                           iter = iter, means = means, seed = seed)
    } else {
      results <- frequentist_unconstrained_lvp(data, variability = variability,
                                               iter = iter, means = means, seed = seed)
    }
    
    colnames(results)[1:2] <- c("ID", "GOF")
    
    # --- CORRECCIÓN ---
    # Identificamos cuál es el nombre del grupo de la mezcla
    groups       <- factor(data[, 2], levels = unique(data[, 2]))
    mixture_name <- levels(groups)[nlevels(groups)]
    
    # Extraemos el ID real de la mezcla directamente de la primera columna de los datos
    actual_id <- as.character(data[data[, 2] == mixture_name, 1][1])
    
    # Asignamos el ID limpio
    results$ID <- actual_id
    # ------------------
    
    return(results)
    
  } else if (tolower(model) == "bayesian") {
    
    message("----------------------------------------------------------------")
    message("Running Bayesian apportionment using the MixSIAR framework.")
    message("Powered by JAGS (Markov Chain Monte Carlo).")
    message("----------------------------------------------------------------")
    
    results <- unmix_mixsiar(data = data, ...)
    return(results)
    
  } else if (tolower(model) == "machine_learning") {
    
    message("----------------------------------------------------------------")
    if (tolower(ml_type) %in% c("xgboost", "xgb")) {
      message("Running Machine Learning apportionment using XGBoost.")
      message("----------------------------------------------------------------")
      results <- unmix_xgboost(data = data, ...)
    } else {
      message("Running Machine Learning apportionment using Random Forest.")
      message("----------------------------------------------------------------")
      results <- unmix_rf(data = data, ...)
    }
    return(results)
    
  } else {
    stop("The 'model' argument must be either 'frequentist', 'bayesian', or 'machine_learning'.")
  }
}
