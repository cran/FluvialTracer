#' Correlation matrix chart
#'
#' Generates a correlation matrix scatterplot of dataset properties, grouped by source, 
#' to aid in collinearity diagnosis and optimal tracer selection.
#'
#' @param data Data frame. Contains source and mixture data in the standard FluvialTracer format.
#' @param tracers Numeric vector. Indices of the tracers to include in the matrix. Default 1:(ncol(data) - 2).
#' @param mixtures Logical. If TRUE, includes mixture samples in the correlation plots. Default TRUE.
#' @param nmixtures Integer. Number of mixture samples present in the dataset. Default 2.
#' @param colors Character vector. Custom colors for the scatterplot groups. If NULL (default),
#'   colors are pulled automatically from fluvial_palette: the first N-1 Source colors, plus one
#'   Mixture color if mixtures = TRUE.
#' @param seed Integer. Random seed for reproducibility. Default NULL.
#'
#' @return A GGally::ggpairs object.
#' 
#' @export
#' 
Correlation_Chart <- function(data, tracers = 1:(ncol(data) - 2),
                              mixtures = TRUE, nmixtures = 2, colors = NULL, seed = NULL) {
  
  # ── Internal helpers ────────────────────────────────────────────────────────
  
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
  
  raw_dataset_internal <- function(d) {
    # Identify structure
    ntracer <- (ncol(d) - 3) / 2
    n_col   <- which(tolower(colnames(d)) == "n")
    
    # Identify mixture rows (last unique group) vs source rows
    d[, 2] <- factor(d[, 2], levels = unique(d[, 2]))
    mixture_label <- levels(d[, 2])[nlevels(d[, 2])]
    is_mixture    <- d[, 2] == mixture_label
    
    raw_data <- data.frame()
    
    # ── Sources: expand using n, mean, sd ──────────────────────────────────
    source_rows <- d[!is_mixture, ]
    for (i in seq_len(nrow(source_rows))) {
      n_samples <- source_rows[i, n_col]
      if (n_samples < 1) n_samples <- 1   # safety fallback
      temp_data <- data.frame(ID    = rep(source_rows[i, 1], n_samples),
                              Group = rep(as.character(source_rows[i, 2]), n_samples))
      for (j in seq_len(ntracer)) {
        mean_val <- source_rows[i, 2 + j]
        sd_val   <- source_rows[i, 2 + ntracer + j]
        temp_data[colnames(source_rows)[2 + j]] <-
          if (n_samples > 1) stats::rnorm(n_samples, mean = mean_val, sd = sd_val)
        else               mean_val
      }
      raw_data <- rbind(raw_data, temp_data)
    }
    
    # ── Mixture: append directly (no expansion, n is irrelevant) ───────────
    mix_rows <- d[is_mixture, ]
    for (i in seq_len(nrow(mix_rows))) {
      temp_data <- data.frame(ID    = mix_rows[i, 1],
                              Group = as.character(mix_rows[i, 2]))
      for (j in seq_len(ntracer)) {
        temp_data[colnames(mix_rows)[2 + j]] <- mix_rows[i, 2 + j]
      }
      raw_data <- rbind(raw_data, temp_data)
    }
    
    # Strip Mean_/mean_ prefix from tracer column names
    colnames(raw_data) <- gsub("^[Mm]ean_", "", colnames(raw_data))
    raw_data
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # ── Pre-processing ──────────────────────────────────────────────────────────
  if (is_averaged_internal(data)) {
    data <- raw_dataset_internal(data)
  }
  
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  data1 <- data
  colnames(data1)[2] <- "Sources"
  
  if (mixtures) {
    if (nmixtures < 3) {
      data1[nrow(data1) + 1, ] <- data1[nrow(data1), ]
      data1[nrow(data1) + 1, ] <- data1[nrow(data1), ]
    }
  } else {
    data1 <- data1[which(!data1[, 2] == data1[nrow(data1), 2]), ]
  }
  
  # Default colors: N-1 Source colors + 1 Mixture color when mixtures = TRUE,
  # or N Source colors only when mixtures = FALSE (mixture row excluded above).
  if (is.null(colors)) {
    grp_levels <- levels(data1$Sources)
    n_groups   <- length(grp_levels)
    if (mixtures) {
      colors <- c(.default_source_colors(n_groups - 1),
                  .default_source_colors(1, mixture = TRUE))
    } else {
      colors <- .default_source_colors(n_groups)
    }
    names(colors) <- grp_levels
  }
  
  # ── Chart ───────────────────────────────────────────────────────────────────
  p <- GGally::ggpairs(data = data1, columns = tracers + 2,
                       upper  = list(continuous = "cor"),
                       lower  = list(combo = "facetdensity"),
                       title  = "Correlation Matrix",
                       mapping = ggplot2::aes(color = Sources, alpha = 0.25),
                       legend = c(1, 1), 
                       progress = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  p <- p +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors)
  
  return(p)
}

# ── Aliases ──────────────────────────────────────────────────────────────────
#' @rdname Correlation_Chart
#' @export
correlationPlot <- Correlation_Chart

#' @rdname Correlation_Chart
#' @export
correlation_plot <- Correlation_Chart
