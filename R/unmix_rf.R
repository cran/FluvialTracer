#' Machine Learning source apportionment using Random Forest (with expanded bounds)
#'
#' Trains a Random Forest model on simulated virtual mixtures. Features a 
#' `stretch_bounds` parameter to expand training data beyond the 0-1 limits,
#' effectively eliminating the boundary bias typically seen in pure Machine Learning
#' and Bayesian unmixing of extreme sources.
#'
#' @param data Data frame. Contains source and mixture data in the standard FluvialTracer format.
#' @param n_sim Integer. Number of virtual mixtures to simulate for training. Default 10000.
#' @param n_trees Integer. Number of trees for the Random Forest model. Default 500.
#' @param stretch_bounds Numeric. Factor to expand training boundaries to avoid 0/1 bias. Default 1.2.
#' @param seed Integer. Random seed for reproducibility. Default 123456.
#'
#' @return A list containing the apportionment data frame (one row per tree) and a summary.
#'
#' @export
#' 
unmix_rf <- function(data, n_sim = 10000L, n_trees = 500L, stretch_bounds = 1.2, seed = 123456L) {
  
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("Please install 'randomForest'.", call. = FALSE)
  set.seed(seed)
  
  # 1. IDENTIFY COLUMNS AND GROUPS
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
  
  message(crayon::cyan(crayon::bold("Starting Random Forest Unmixing:")))
  message("Simulating ", n_sim, " virtual mixtures (Stretch factor: ", stretch_bounds, ")...")
  
  # 2. CALCULATE SOURCE DISTRIBUTIONS
  source_stats <- list()
  for (src in source_names) {
    src_data <- sources_raw[sources_raw[[group_col]] == src, tracers]
    source_stats[[src]] <- list(mean = colMeans(src_data, na.rm = TRUE), sd = apply(src_data, 2, sd, na.rm = TRUE))
    source_stats[[src]]$sd[is.na(source_stats[[src]]$sd) | source_stats[[src]]$sd == 0] <- 1e-6
  }
  
  # 3. GENERATE SYNTHETIC TRAINING DATA
  w <- matrix(stats::rgamma(n_sim * n_sources, shape = 1, scale = 1), nrow = n_sim, ncol = n_sources)
  w <- w / rowSums(w)
  
  # --- THE FIX: STRETCH BOUNDARIES ---
  if (stretch_bounds > 1.0) {
    center <- 1.0 / n_sources
    w <- center + stretch_bounds * (w - center)
  }
  colnames(w) <- source_names
  
  synth_mixtures <- matrix(0, nrow = n_sim, ncol = n_tracers)
  colnames(synth_mixtures) <- tracers
  
  for (i in 1:n_sources) {
    src_weight <- w[, i]
    for (j in 1:n_tracers) {
      tr_name <- tracers[j]
      random_src_values <- stats::rnorm(n_sim, mean = source_stats[[source_names[i]]]$mean[tr_name], 
                                        sd = source_stats[[source_names[i]]]$sd[tr_name])
      synth_mixtures[, j] <- synth_mixtures[, j] + (src_weight * random_src_values)
    }
  }
  synth_data <- as.data.frame(synth_mixtures)
  
  # 4. TRAIN THE RANDOM FOREST MODELS
  rf_models <- list()
  message("Training Random Forest models (", n_trees, " trees per source)...")
  pb <- txtProgressBar(min = 0, max = n_sources, style = 3, width = 50, char = "=")
  for (i in 1:n_sources) {
    rf_models[[source_names[i]]] <- randomForest::randomForest(x = synth_data, y = w[, source_names[i]], ntree = n_trees)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 5. PREDICT ACTUAL MIXTURES
  mix_tracers <- mix_df[, tracers, drop = FALSE]
  raw_predictions <- matrix(0, nrow = n_trees, ncol = n_sources)
  colnames(raw_predictions) <- source_names
  
  for (i in 1:n_sources) {
    pred_obj <- stats::predict(rf_models[[source_names[i]]], newdata = mix_tracers, predict.all = TRUE)
    raw_predictions[, i] <- pred_obj$individual[1, ]
  }
  
  raw_predictions <- raw_predictions / rowSums(raw_predictions)
  
  # 6. FORMAT OUTPUT
  apportionment_df <- data.frame(ID = rep(as.character(mix_df[[id_col]]), n_trees), GOF = NA)
  apportionment_df <- cbind(apportionment_df, as.data.frame(raw_predictions))
  
  summary_df <- data.frame(
    mean = colMeans(raw_predictions), sd = apply(raw_predictions, 2, sd),
    `2.5%` = apply(raw_predictions, 2, quantile, probs = 0.025), `50%` = apply(raw_predictions, 2, median),
    `97.5%` = apply(raw_predictions, 2, quantile, probs = 0.975), check.names = FALSE
  )
  
  message(crayon::green(crayon::bold("Summary of the model outputs:")))
  message(paste(utils::capture.output(print(summary_df)), collapse = "\n"))
  invisible(list(apportionment = apportionment_df, rf_summary = summary_df))
}