#' Bayesian source apportionment using MixSIAR
#'
#' A wrapper function that seamlessly formats an `FluvialTracer` dataset (raw or averaged) 
#' and passes it to the `MixSIAR` Bayesian unmixing engine. 
#'
#' @param data Data frame. Contains source and mixture data in the standard FluvialTracer format.
#' @param tracers Character vector. Names of the tracers to include in the model.
#' @param prior_type Character. "uninformative" (default, alpha=1) or "informative".
#' @param prior_values Numeric vector. Required if `prior_type = "informative"`. Must match the number of sources.
#' @param run Character. MCMC run length. Options: "test", "very short", "short", "normal" (default), "long", "very long", "extreme".
#' @param process_err Logical. Include process error in the JAGS model? Default FALSE.
#' @param resid_err Logical. Include residual error in the JAGS model? Default TRUE.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{apportionment}}{Data frame with the posterior draws of the
#'     global source proportions, in the FluvialTracer format (\code{ID},
#'     \code{GOF} and one column per source), compatible with
#'     \code{\link{plot_apportionment}}.}
#'   \item{\code{bayesian_summary}}{Data frame with the summary statistics
#'     of all model parameters returned by JAGS.}
#'   \item{\code{jags_model}}{The full JAGS model object returned by
#'     \code{MixSIAR::run_model}.}
#' }
#'
#' @export
#' 
unmix_mixsiar <- function(data, 
                          tracers, 
                          prior_type = "uninformative", 
                          prior_values = NULL, 
                          run = "normal", 
                          process_err = FALSE, 
                          resid_err = TRUE) {
  
  # Ensure required packages are installed
  if (!requireNamespace("MixSIAR", quietly = TRUE)) {
    stop("The 'MixSIAR' package is required. Please install it using install.packages('MixSIAR').", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required.", call. = FALSE)
  }
  
  # ── JAGS DEPENDENCY CHECK ─────────────────────────────────────────────────
  # If the standalone JAGS software is missing, rjags will fail to load.
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop(paste0("\n======================================================================\n",
                "CRITICAL ERROR: 'rjags' package or standalone JAGS software is missing.\n",
                "======================================================================\n",
                "The Bayesian model requires the standalone JAGS software to be \n",
                "installed on your operating system (not just the R package).\n\n",
                "Please download and install JAGS from the following link:\n",
                "-> https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/\n\n",
                "Once installed, restart R/RStudio and try again.\n",
                "======================================================================\n"), 
         call. = FALSE)
  }
  
  # ── Auto-detect tracers if not provided ───────────────────────────────────
  if (missing(tracers) || is.null(tracers)) {
    is_avg <- any(grepl("^[Mm]ean_", names(data)))
    if (is_avg) {
      # If averaged, grab all columns that start with Mean_ and strip the prefix
      mean_cols <- grep("^[Mm]ean_", names(data), value = TRUE)
      tracers <- gsub("^[Mm]ean_", "", mean_cols)
    } else {
      # If raw, tracers are simply everything from column 3 onwards
      tracers <- names(data)[3:ncol(data)]
    }
  }
  # Dynamically grab the name of the group column (e.g., "Source", "samples", "Group")
  group_col <- names(data)[2]
  
  # Check if data is already averaged or if it is raw
  is_avg <- any(grepl("^[Mm]ean_", names(data)))
  
  if (!is_avg) {
    # ── RAW DATA: Calculate Means, SDs, and n automatically ──────────────
    data[[group_col]] <- factor(data[[group_col]], levels = unique(data[[group_col]]))
    mix_label <- levels(data[[group_col]])[nlevels(data[[group_col]])]
    
    sources_raw <- data[data[[group_col]] != mix_label, ]
    mix_df <- data[data[[group_col]] == mix_label, c(group_col, tracers), drop = FALSE]
    
    # Aggregate sources
    source_df <- sources_raw %>%
      dplyr::group_by(!!dplyr::sym(group_col)) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(tracers), list(Mean = ~mean(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE))),
        n = dplyr::n(),
        .groups = 'drop'
      )
    
    # Rename columns to strictly match MixSIAR format (MeanMg, SDMg)
    names(source_df) <- gsub("^(.*)_Mean$", "Mean\\1", names(source_df))
    names(source_df) <- gsub("^(.*)_SD$", "SD\\1", names(source_df))
    source_names <- source_df[[1]]
    
  } else {
    # ── AVERAGED DATA: Just rename the columns ───────────────────────────
    mean_cols <- paste0("Mean_", tracers)
    sd_cols <- paste0("SD_", tracers)
    
    # Fallback if lowercase was used
    if (!any(mean_cols %in% names(data))) {
      mean_cols <- paste0("mean_", tracers)
      sd_cols <- paste0("sd_", tracers)
    }
    
    data_unmix <- data[, c("ID", group_col, mean_cols, sd_cols, "n")]
    data_unmix[[group_col]] <- factor(data_unmix[[group_col]], levels = unique(data_unmix[[group_col]]))
    mix_label <- levels(data_unmix[[group_col]])[nlevels(data_unmix[[group_col]])]
    
    source_df <- data_unmix[data_unmix[[group_col]] != mix_label, -1] # Drop ID
    mix_df <- data_unmix[data_unmix[[group_col]] == mix_label, -1, drop = FALSE] # Drop ID
    
    names(source_df)[names(source_df) %in% mean_cols] <- paste0("Mean", tracers)
    names(source_df)[names(source_df) %in% sd_cols]   <- paste0("SD", tracers)
    source_names <- source_df[[1]]
    
    mix_df <- mix_df[, !names(mix_df) %in% c(sd_cols, "n"), drop = FALSE]
    names(mix_df)[names(mix_df) %in% mean_cols] <- tracers
  }
  
  # ── Build TEF (Discrimination) Table ──────────────────────────────────────
  TEF_df <- source_df[, c(names(source_df)[1], paste0("Mean", tracers))]
  for (t in tracers) TEF_df[[paste0("Mean", t)]] <- 0
  for (t in tracers) TEF_df[[paste0("SD", t)]] <- 0
  
  # ── Temporary CRAN-Safe Files ─────────────────────────────────────────────
  temp_dir <- tempdir()
  src_file <- file.path(temp_dir, "temp_source.csv")
  mix_file <- file.path(temp_dir, "temp_mix.csv")
  tef_file <- file.path(temp_dir, "temp_tef.csv")
  model_file <- file.path(temp_dir, "temp_model.txt")
  
  utils::write.table(source_df, src_file, sep = ",", row.names = FALSE)
  utils::write.table(mix_df, mix_file, sep = ",", row.names = FALSE)
  utils::write.table(TEF_df, tef_file, sep = ",", row.names = FALSE)
  
  # ── Load into MixSIAR ─────────────────────────────────────────────────────
  message("Formatting data for MixSIAR...")
  suppressMessages({
    mix_obj <- MixSIAR::load_mix_data(filename = mix_file, 
                                      iso_names = tracers, 
                                      factors = NULL, fac_random = NULL, fac_nested = NULL, cont_effects = NULL)
    
    source_obj <- MixSIAR::load_source_data(filename = src_file, 
                                            source_factors = NULL, conc_dep = FALSE, data_type = "means", mix_obj)
    
    discr_obj <- MixSIAR::load_discr_data(filename = tef_file, mix_obj)
  })
  
  # ── Priors ────────────────────────────────────────────────────────────────
  if (prior_type == "uninformative") {
    alpha <- rep(1, source_obj$n.sources)
  } else if (prior_type == "informative") {
    if (is.null(prior_values)) stop("Provide 'prior_values' for informative prior.", call. = FALSE)
    alpha <- prior_values
  }
  alpha <- alpha * length(alpha) / sum(alpha)
  alpha[alpha == 0] <- 0.01
  
  # ── Run JAGS ──────────────────────────────────────────────────────────────
  message(paste0("Running MixSIAR MCMC chain (Length: ", run, ")... This may take a while."))
  MixSIAR::write_JAGS_model(model_file, process_err, resid_err, mix_obj, source_obj)
  
  jags_out <- MixSIAR::run_model(run = run, mix_obj, source_obj, discr_obj, model_file, 
                                 alpha.prior = alpha, process_err, resid_err)
  
  # Clean up temp files
  file.remove(src_file, mix_file, tef_file, model_file)
  
  # ── Extract Proportions ───────────────────────────────────────────────────
  p_global <- as.data.frame(jags_out$BUGSoutput$sims.list$p.global)
  
  # 1. Assign names using the alphabetical order that MixSIAR forced internally
  colnames(p_global) <- source_obj$source_names
  
  # 2. Reorder the columns back to the original order of your dataset 
  #    (so plot colors will match the Frequentist and ML models)
  p_global <- p_global[, as.character(source_names), drop = FALSE]
  
  # Grab the actual name of the mixture from the dataset
  mix_id <- as.character(mix_df[[group_col]][1])
  
  # Build the FluvialTracer-compatible dataframe
  apportionment_df <- data.frame(
    ID = mix_id,
    GOF = NA # No exact frequentist equivalent for individual MCMC iterations
  )
  apportionment_df <- cbind(apportionment_df, p_global)
  
  message("MixSIAR run complete.")
  
  return(list(
    apportionment = apportionment_df,
    bayesian_summary = as.data.frame(jags_out$BUGSoutput$summary),
    jags_model = jags_out
  ))
}
