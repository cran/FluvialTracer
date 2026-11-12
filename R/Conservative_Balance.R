# Internal helper - not exported
.has_means <- function(data) {
  nms <- colnames(data)
  (any(grepl("^[Mm]ean_", nms)) && any(grepl("^[Ss][Dd]_", nms))) ||
    (any(grepl("^mean_",    nms)) && any(grepl("^sd_",        nms)))
}

# Internal helper - not exported
.normalise_iso_colnames <- function(data) {
  nms <- colnames(data)
  if (any(grepl("^mean_", nms)) || any(grepl("^sd_", nms))) {
    nms <- sub("^mean_", "Mean_", nms)
    nms <- sub("^sd_",   "SD_",   nms)
    colnames(data) <- nms
  }
  data
}

#' Conservative Balance (CB) transformation for isotopic tracers
#'
#' Transforms isotopic ratio and elemental content data into virtual scalar tracers 
#' that can be directly used in classical unmixing models.
#'
#' @param data Data frame. Contains the isotopic tracer characteristics of sediment sources and one mixture. Both averaged and raw formats are supported.
#' @param error_propagation Character. Error propagation method for averaged data: "ratio" (default) or "full".
#' @param impute_nas Logical. If TRUE, automatically imputes missing values (NAs) in source samples using their respective group means. Default FALSE.
#' 
#' @return A data frame with isotopic tracers replaced by their CB-transformed scalar equivalents. 
#' 
#' @export
#' 
Conservative_Balance <- function(data, error_propagation = "ratio", impute_nas = FALSE) {
  
  # ── 1. Input checks & Setup ───────────────────────────────────────────────
  if (!error_propagation %in% c("ratio", "full"))
    stop("'error_propagation' must be either \"ratio\" or \"full\".")
  
  mix_data <- inputSample(data)
  if (nrow(mix_data) == 0)
    stop("The dataset does not contain a mixture. The mixture's isotopic ratios are required.")
  
  if (nrow(mix_data) > 1) {
    warning("Multiple mixtures detected - averaging into a single representative mixture.")
    mix_avg        <- colMeans(mix_data[, -1L, drop = FALSE], na.rm = TRUE)
    mix_row        <- mix_data[1L, ]
    mix_row[, -1L] <- mix_avg
    data           <- rbind(data[seq_len(nrow(data) - nrow(mix_data)), ], mix_row)
  }
  
  data <- .normalise_iso_colnames(data)
  col_names <- colnames(data)
  is_avg <- .has_means(data)
  mix_row_idx <- nrow(data)
  
  # ── 2. Dynamic Column Detection (Order Independent) ───────────────────────
  if (is_avg) {
    c_mean_cols <- sort(grep("^[Mm]ean_.*([Cc]onc|[Cc]ont)", col_names, value = TRUE))
    c_sd_cols   <- sort(grep("^[Ss][Dd]_.*([Cc]onc|[Cc]ont)", col_names, value = TRUE))
    
    r_mean_cols <- sort(setdiff(grep("^[Mm]ean_", col_names, value = TRUE), c_mean_cols))
    r_sd_cols   <- sort(setdiff(grep("^[Ss][Dd]_", col_names, value = TRUE), c_sd_cols))
    
    if (length(r_mean_cols) == 0 || length(c_mean_cols) == 0) {
      stop("Could not automatically identify Ratio and Content columns. Please ensure content columns contain 'Conc' or 'Cont'.")
    }
    
    if (error_propagation == "full" && length(c_sd_cols) == 0) {
      warning("error_propagation = 'full' requires SD of contents, but none were found. Falling back to 'ratio' propagation.")
      error_propagation <- "ratio"
    }
  } else {
    c_cols <- sort(grep("[Cc]onc|[Cc]ont", col_names, value = TRUE))
    r_cols <- sort(setdiff(col_names[-c(1, 2)], c(c_cols, "n")))
    
    if (length(r_cols) == 0 || length(c_cols) == 0) {
      stop("Could not automatically identify Ratio and Content columns. Please ensure content columns contain 'Conc' or 'Cont'.")
    }
  }
  
  # ── 3. Handle Missing Mixture Ratios ──────────────────────────────────────
  if (is_avg) {
    mix_ratios <- data[mix_row_idx, r_mean_cols]
  } else {
    mix_ratios <- data[mix_row_idx, r_cols]
  }
  
  na_tracers <- is.na(mix_ratios)
  if (any(na_tracers)) {
    warning("The mixture is missing ratios for some tracers. These tracers have been removed from the CB calculation.")
    if (is_avg) {
      r_mean_cols <- r_mean_cols[!na_tracers]
      r_sd_cols   <- r_sd_cols[!na_tracers]
      c_mean_cols <- c_mean_cols[!na_tracers]
      if (length(c_sd_cols) > 0) c_sd_cols <- c_sd_cols[!na_tracers]
    } else {
      r_cols <- r_cols[!na_tracers]
      c_cols <- c_cols[!na_tracers]
    }
  }
  
  # ── 4. Impute NAs in Sources ──────────────────────────────────────────────
  if (impute_nas) {
    source_idx <- seq_len(mix_row_idx - 1L)
    groups <- unique(data[[2]][source_idx])
    imputed_count <- 0
    for (g in groups) {
      g_idx <- which(data[[2]][source_idx] == g)
      for (col in 3L:ncol(data)) {
        na_mask <- is.na(data[g_idx, col])
        if (any(na_mask)) {
          g_mean <- mean(data[g_idx, col], na.rm = TRUE)
          if (!is.nan(g_mean)) {
            data[g_idx[na_mask], col] <- g_mean
            imputed_count <- imputed_count + sum(na_mask)
          }
        }
      }
    }
    if (imputed_count > 0) message(sprintf("Imputed %d missing value(s) in sources using group means.", imputed_count))
  }
  
  # ── 5. Perform the CB Transformation ──────────────────────────────────────
  cb_df <- data.frame(ID = data[[1]], Sources = data[[2]])
  
  if (is_avg) {
    tr_names <- gsub("^[Mm]ean_R_|^[Mm]ean_|[Rr]atio_", "", r_mean_cols)
    
    # Use lists to accumulate Means and SDs separately to ensure proper column order
    mean_list <- list()
    sd_list <- list()
    
    for (i in seq_along(r_mean_cols)) {
      r_m <- data[[r_mean_cols[i]]]
      r_s <- data[[r_sd_cols[i]]]
      c_m <- data[[c_mean_cols[i]]]
      mix_ref <- data[mix_row_idx, r_mean_cols[i]]
      
      # Mean CB
      cb_mean <- c_m * (r_m - mix_ref)
      
      # SD CB
      if (error_propagation == "full") {
        c_s <- data[[c_sd_cols[i]]]
        cb_sd <- sqrt((c_s * (r_m - mix_ref))^2 + (c_m * r_s)^2)
      } else {
        cb_sd <- r_s * c_m
      }
      
      mean_list[[paste0("Mean_CB_", tr_names[i])]] <- cb_mean
      sd_list[[paste0("SD_CB_", tr_names[i])]]     <- cb_sd
    }
    
    # Bind them together in the exact order FluvialTracer expects
    cb_df <- cbind(cb_df, as.data.frame(mean_list), as.data.frame(sd_list))
    
    # Format mixture row perfectly
    cb_df[nrow(cb_df), grep("^Mean_", names(cb_df))] <- 0
    cb_df[nrow(cb_df), grep("^SD_", names(cb_df))] <- NA
    
  } else {
    tr_names <- gsub("^R_|^[Rr]atio_", "", r_cols)
    
    for (i in seq_along(r_cols)) {
      r_v <- data[[r_cols[i]]]
      c_v <- data[[c_cols[i]]]
      mix_ref <- data[mix_row_idx, r_cols[i]]
      
      cb_df[[paste0("CB_", tr_names[i])]] <- c_v * (r_v - mix_ref)
    }
    
    cb_df[nrow(cb_df), grep("^CB_", names(cb_df))] <- 0
  }
  
  if ("n" %in% names(data)) cb_df$n <- data$n
  
  return(cb_df)
}

#' @rdname Conservative_Balance
#' @export
CB_Method <- function(data, ...) Conservative_Balance(data, ...)

#' @rdname Conservative_Balance
#' @param ... Additional arguments passed to other methods.
#' @export
cb_method <- function(...) Conservative_Balance(...)