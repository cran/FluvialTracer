#' Principal Component Analysis biplot
#'
#' Performs a PCA on the tracer data and displays a biplot with source group
#' ellipses. Averaged datasets (with \code{Mean_}/\code{SD_} or 
#' \code{mean_}/\code{sd_} columns) are automatically expanded to individual
#' observations by sampling from a normal distribution before the PCA is
#' computed, so all formats are accepted.
#'
#' @param data A data frame containing source and mixture data in the standard
#'   fingerprinting format (raw or averaged).
#' @param components A numeric vector of length 2 indicating which principal
#'   components to display. Default \code{c(1, 2)}.
#' @param colors An optional character vector of colors, one per group. If
#'   \code{NULL} (default), colors are pulled automatically from fluvial_palette:
#'   the first N-1 Source colors, plus one Mixture color for the last group.
#' @param seed An integer value used to initialize the random number generator.
#'
#' @return A ggplot object.
#'
#' @export
#' 
plot_pca <- function(data, components = c(1L, 2L), colors = NULL, seed = NULL) {
  
  # ── Internal helpers ────────────────────────────────────────────────────────
  
  # Detects averaged format regardless of prefix case (Mean_/SD_ or mean_/sd_)
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
  
  # ── If averaged, expand to synthetic raw observations ──────────────────────
  if (is_averaged_internal(data)) {
    src      <- inputSource(data)   # always returns Mean_/SD_/n layout
    n_col    <- ncol(src)
    tracer_n <- as.integer((n_col - 2L) / 2L)
    
    rows <- list()
    for (i in seq_len(nrow(src))) {
      n_obs <- as.integer(src[i, n_col])
      means <- as.numeric(src[i, 2L:(1L + tracer_n)])
      sds   <- as.numeric(src[i, (2L + tracer_n):(1L + 2L * tracer_n)])
      mat   <- mapply(function(m, s) stats::rnorm(n_obs, m, s), means, sds)
      if (n_obs == 1L) mat <- matrix(mat, nrow = 1L)
      df <- as.data.frame(mat)
      colnames(df) <- gsub("^[Mm]ean_", "", colnames(src)[2L:(1L + tracer_n)])
      df$samples <- as.character(src[i, 1L])
      rows[[i]] <- df
    }
    
    # Append mixture row (single observation, no SD needed)
    mix    <- inputSample(data)
    mix_df <- as.data.frame(mix[, -1L, drop = FALSE])
    colnames(mix_df) <- gsub("^[Mm]ean_", "", colnames(mix_df))
    mix_df$samples <- as.character(mix[1L, 1L])
    rows[[nrow(src) + 1L]] <- mix_df
    
    raw      <- do.call(rbind, rows)
    data_pca <- raw[, !colnames(raw) %in% "samples"]
    groups   <- factor(raw$samples, levels = unique(raw$samples))
    
  } else {
    # ── Raw format ─────────────────────────────────────────────────────────
    data[, 2L] <- factor(data[, 2L], levels = unique(data[, 2L]))
    data_pca   <- data[, 3L:ncol(data)]
    groups     <- data[, 2L]
  }
  
  # ── PCA ───────────────────────────────────────────────────────────────────
  pca_fit <- prcomp(data_pca, center = TRUE, scale. = TRUE)
  
  pca_scores        <- as.data.frame(pca_fit$x[, components])
  pca_scores$groups <- groups
  
  if (any(is.na(pca_scores))) {
    message("NA values in PCA scores were replaced with the group minimum.")
    pca_scores <- pca_scores |>
      dplyr::group_by(groups) |>
      dplyr::mutate(dplyr::across(
        dplyr::where(is.numeric),
        ~ ifelse(is.na(.), min(., na.rm = TRUE), .)
      )) |>
      as.data.frame()
  }
  
  # Default colors: N-1 Source colors + 1 Mixture color, matching the group
  # levels (mixture assumed to be the last level, per package convention).
  num_groups <- length(levels(groups))
  if (is.null(colors)) {
    colors <- c(.default_source_colors(num_groups - 1),
                .default_source_colors(1, mixture = TRUE))
    names(colors) <- levels(groups)
  } else {
    if (length(colors) < num_groups)
      stop("Not enough colors: need ", num_groups,
           " but only ", length(colors), " supplied.")
  }
  
  # ── Plot ──────────────────────────────────────────────────────────────────
  p <- ggbiplot(pca_fit,
                choices        = components,
                obs.scale      = 1,
                var.scale      = 1,
                ellipse        = FALSE,
                groups         = pca_scores$groups,
                ellipse.prob   = 0.95,
                circle         = FALSE,
                varname.size   = 5,
                varname.adjust = 3,
                var.axes       = TRUE) +
    geom_point(aes(shape = groups, colour = groups), size = 3) +
    geom_hline(yintercept = 0, colour = "black",
               linetype = "longdash", linewidth = 1) +
    geom_vline(xintercept = 0, colour = "black",
               linetype = "longdash", linewidth = 1) +
    stat_ellipse(type = "t", linewidth = 1,
                 geom = "polygon", alpha = 0.2,
                 aes(fill = groups), level = 0.85) +
    stat_ellipse(type = "t", linewidth = 1,
                 aes(colour = groups), level = 0.85) +
    guides(groups = "none") +
    ggtitle("PCA") +
    theme(legend.direction = "horizontal",
          legend.position  = "top",
          text             = element_text(size = 15),
          plot.title       = element_text(hjust = 0.5)) +
    scale_colour_manual(values = colors) +
    scale_fill_manual(values   = colors)
  
  return(p)
}

# ── Aliases ──────────────────────────────────────────────────────────────────
#' @rdname plot_pca
#' @export
PCAPlot <- plot_pca

#' @rdname plot_pca
#' @export
PCA_plot <- plot_pca