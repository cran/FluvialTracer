#' Machine Learning source apportionment using XGBoost (Gradient Boosting)
#'
#' Trains an ensemble of XGBoost models on simulated virtual mixtures. 
#' Uses bootstrapping to generate uncertainty distributions (density curves).
#'
#' @param data Data frame. Contains source and mixture data.
#' @param n_sim Integer. Number of virtual mixtures for the base training set. Default 5000.
#' @param n_boot Integer. Number of XGBoost models to train (iterations for uncertainty). Default 100.
#' @param stretch_bounds Numeric. Expansion factor for training boundaries. Default 1.2.
#' @param seed Integer. Random seed. Default 123456..
#'
#' @return Invisibly returns a list with two elements: \code{apportionment},
#'   a data frame with one row per bootstrapped model (\code{ID}, \code{GOF}
#'   and one column per source), compatible with
#'   \code{\link{plot_apportionment}}; and \code{xgb_summary}, a data frame
#'   with the mean, SD, median and 95\% interval of each source contribution.
#'
#' @importFrom utils read.csv capture.output
#' @importFrom dplyr slice bind_rows
#'
#' @export
#' 
unmix_xgboost <- function(data, n_sim = 5000L, n_boot = 100L, stretch_bounds = 1.2, seed = 123456L) {
  
  if (!requireNamespace("xgboost", quietly = TRUE)) stop("Please install 'xgboost'.", call. = FALSE)
  set.seed(seed)
  
  id_col <- names(data)[1]
  group_col <- names(data)[2]
  tracers <- names(data)[3:ncol(data)]
  
  data[[group_col]] <- factor(data[[group_col]], levels = unique(data[[group_col]]))
  mix_label <- levels(data[[group_col]])[nlevels(data[[group_col]])]
  
  sources_raw <- data[data[[group_col]] != mix_label, ]
  mix_df <- data[data[[group_col]] == mix_label, ]
  source_names <- unique(as.character(sources_raw[[group_col]]))
  n_sources <- length(source_names)
  n_tracers <- length(tracers)
  
  message(crayon::cyan(crayon::bold("Starting XGBoost Unmixing:")))
  message("Training ", n_boot, " bootstrapped XGBoost models...")
  
  # Calculate source distributions
  source_stats <- list()
  for (src in source_names) {
    src_data <- sources_raw[sources_raw[[group_col]] == src, tracers]
    source_stats[[src]] <- list(mean = colMeans(src_data, na.rm = TRUE), sd = apply(src_data, 2, sd, na.rm = TRUE))
    source_stats[[src]]$sd[is.na(source_stats[[src]]$sd) | source_stats[[src]]$sd == 0] <- 1e-6
  }
  
  # Generate base virtual mixtures
  w <- matrix(stats::rgamma(n_sim * n_sources, shape = 1, scale = 1), nrow = n_sim, ncol = n_sources)
  w <- w / rowSums(w)
  if (stretch_bounds > 1.0) {
    center <- 1.0 / n_sources
    w <- center + stretch_bounds * (w - center)
  }
  colnames(w) <- source_names
  
  synth_mix <- matrix(0, nrow = n_sim, ncol = n_tracers)
  colnames(synth_mix) <- tracers
  for (i in 1:n_sources) {
    src_weight <- w[, i]
    for (j in 1:n_tracers) {
      tr_name <- tracers[j]
      r_vals <- stats::rnorm(n_sim, mean = source_stats[[source_names[i]]]$mean[tr_name], sd = source_stats[[source_names[i]]]$sd[tr_name])
      synth_mix[, j] <- synth_mix[, j] + (src_weight * r_vals)
    }
  }
  
  mix_tracers <- as.matrix(mix_df[, tracers, drop = FALSE])
  raw_predictions <- matrix(0, nrow = n_boot, ncol = n_sources)
  colnames(raw_predictions) <- source_names
  
  pb <- txtProgressBar(min = 0, max = n_boot, style = 3, width = 50, char = "=")
  
  # Train bootstrapped XGBoost models
  for (boot_idx in 1:n_boot) {
    # Resample training data for this iteration
    boot_indices <- sample(1:n_sim, replace = TRUE)
    train_x <- synth_mix[boot_indices, ]
    
    for (i in 1:n_sources) {
      train_y <- w[boot_indices, source_names[i]]
      # Suppress verbose output from xgboost
      capture.output({
        xgb_model <- xgboost::xgboost(data = train_x, label = train_y, nrounds = 50, 
                                      objective = "reg:squarederror", max_depth = 4, eta = 0.1, verbose = 0)
      })
      raw_predictions[boot_idx, i] <- stats::predict(xgb_model, newdata = mix_tracers)
    }
    setTxtProgressBar(pb, boot_idx)
  }
  close(pb)
  
  
  raw_predictions <- raw_predictions / rowSums(raw_predictions)
  
  apportionment_df <- data.frame(ID = rep(as.character(mix_df[[id_col]]), n_boot), GOF = NA)
  apportionment_df <- cbind(apportionment_df, as.data.frame(raw_predictions))
  
  summary_df <- data.frame(
    mean = colMeans(raw_predictions), sd = apply(raw_predictions, 2, sd),
    `2.5%` = apply(raw_predictions, 2, quantile, probs = 0.025), `50%` = apply(raw_predictions, 2, median),
    `97.5%` = apply(raw_predictions, 2, quantile, probs = 0.975), check.names = FALSE
  )
  
  message(crayon::green(crayon::bold("Summary of the model outputs:")))
  message(paste(utils::capture.output(print(summary_df)), collapse = "\n"))
  invisible(list(apportionment = apportionment_df, xgb_summary = summary_df))
}