#' Linear discriminant analysis chart
#'
#' Performs a linear discriminant analysis (LDA) to project data into maximum-separation
#' dimensions. Includes a comprehensive discrimination summary detailing variance explained,
#' a leave-one-out cross-validation confusion matrix, and misclassification rates to
#' rigorously assess source separability.
#'
#' @param data Data frame. Contains source and mixture data in the standard FluvialTracer format.
#' @param Plot3D Logical. If TRUE, renders a 3D scatter plot instead of a 2D projection. Default FALSE.
#' @param text Logical. If TRUE, overlays the sample identification numbers on the plot. Default FALSE.
#' @param colors Character vector. Custom colors for the source groups. If NULL (default),
#'   the first N Source colors from fluvial_palette are used automatically, where N is the
#'   number of source groups detected in data.
#' @param interactive Logical. If TRUE, renders an interactive HTML widget for 3D plots. Default FALSE.
#' @param seed Integer. Random seed for reproducibility. Default NULL.
#'
#' @return Invisibly returns a list with: \code{lda} (the fitted
#'   \code{MASS::lda} object), \code{confusion} (leave-one-out confusion
#'   matrix), \code{prop_var} (percentage of variance explained by each
#'   discriminant axis) and \code{plot} (the plot object).
#'
#' @export
#'
LDA_Chart <- function(data, Plot3D = FALSE, text = FALSE, colors = NULL, interactive = FALSE, seed = NULL) {
  
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
    ntracer <- (ncol(d) - 3) / 2
    n_col   <- which(tolower(colnames(d)) == "n")
    d[, 2] <- factor(d[, 2], levels = unique(d[, 2]))
    mixture_label <- levels(d[, 2])[nlevels(d[, 2])]
    is_mixture    <- d[, 2] == mixture_label
    raw_data <- data.frame()
    source_rows <- d[!is_mixture, ]
    for (i in seq_len(nrow(source_rows))) {
      n_samples <- source_rows[i, n_col]
      if (n_samples < 1) n_samples <- 1
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
    mix_rows <- d[is_mixture, ]
    for (i in seq_len(nrow(mix_rows))) {
      temp_data <- data.frame(ID    = mix_rows[i, 1],
                              Group = as.character(mix_rows[i, 2]))
      for (j in seq_len(ntracer)) {
        temp_data[colnames(mix_rows)[2 + j]] <- mix_rows[i, 2 + j]
      }
      raw_data <- rbind(raw_data, temp_data)
    }
    colnames(raw_data) <- gsub("^[Mm]ean_", "", colnames(raw_data))
    raw_data
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # ── Pre-processing ──────────────────────────────────────────────────────────
  if (is_averaged_internal(data)) {
    data <- raw_dataset_internal(data)
  }
  
  # ── Plotting ────────────────────────────────────────────────────────────────
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  groups  <- data[, 2]
  mixture <- levels(groups)[nlevels(groups)]
  sources <- data[!groups == mixture, ]
  s_groups <- droplevels(sources[, 2])
  
  data.lda1 <- sources[3:ncol(sources)]
  data.lda1$groups <- as.factor(s_groups)
  
  if (any(is.na(data.lda1))) {
    message("NA values were filled with the minimum value of the same group.")
  }
  
  data.lda1 <- data.lda1 |>
    dplyr::group_by(groups) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.na(.), min(., na.rm = TRUE), .)))
  
  data.lda1$groups <- as.factor(s_groups)
  data.lda      <- MASS::lda(groups ~ ., data = data.lda1)
  data.lda.pred <- predict(data.lda)
  data.lda.temp <- data.frame(data.lda.pred$x, Sources = data.lda1$groups)
  
  num_groups <- length(unique(data.lda.temp$Sources))
  # Default colors: first N Source colors from fluvial_palette, named to
  # match the actual source group levels.
  if (is.null(colors)) {
    default_colors <- .default_source_colors(num_groups)
    names(default_colors) <- levels(data.lda.temp$Sources)
  } else {
    if (length(colors) < num_groups) stop("Insufficient number of colors provided.")
    default_colors <- colors
  }
  
  if (Plot3D) {
    if (interactive) {
      plot <- plotly::plot_ly(data.lda.temp, x = ~LD1, y = ~LD2, z = ~LD3,
                              type = "scatter3d", mode = "markers",
                              color = ~Sources, colors = default_colors,
                              text = ~Sources, marker = list(size = 8, opacity = 0.7)) |>
        plotly::add_markers()
    } else {
      plot <- with(data.lda.temp,
                   car::scatter3d(LD1, LD2, LD3,
                                  groups = Sources,
                                  surface.col = unname(default_colors)[seq_len(num_groups)],
                                  speed = 8, bg.col = "white",
                                  model.summary = TRUE, surface.alpha = 0.2,
                                  ellipsoid = TRUE, ellipsoid.alpha = 0.3, level = 0.8))
    }
  } else {
    base_plot <- ggplot(data.lda.temp, aes(LD1, LD2, colour = Sources, fill = Sources)) +
      geom_point(aes(shape = Sources), size = 3, alpha = 0.7) +
      geom_hline(yintercept = 0, colour = "black", linetype = "longdash") +
      geom_vline(xintercept = 0, colour = "black", linetype = "longdash") +
      stat_ellipse(type = "t", linewidth = 1, geom = "polygon", alpha = 0.2, level = 0.9) +
      stat_ellipse(type = "t", linewidth = 1, alpha = 0.7, level = 0.9) +
      ggtitle("LDA") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = default_colors) +
      scale_fill_manual(values = default_colors)
    
    if (text) {
      data.lda.temp$row_num <- seq_len(nrow(data.lda.temp))
      plot <- base_plot +
        geom_text(data = data.lda.temp,
                  aes(label = row_num, hjust = 2, vjust = 0.5), size = 4)
    } else {
      plot <- base_plot
    }
  }
  
  print(plot)
  
  # ── Discrimination summary ──────────────────────────────────────────────────
  
  # 1. Variance explained by each discriminant axis
  prop_var <- data.lda$svd^2 / sum(data.lda$svd^2) * 100
  ld_labels <- paste0("LD", seq_along(prop_var))
  message("\n--- Variance explained by each discriminant axis -----------------")
  for (k in seq_along(prop_var)) {
    message(sprintf("  %-6s  %.1f%%", ld_labels[k], prop_var[k]))
  }
  
  # 2. Leave-one-out cross-validation confusion matrix
  data.lda.cv <- MASS::lda(groups ~ ., data = data.lda1, CV = TRUE)
  conf <- table(Actual = data.lda1$groups, Predicted = data.lda.cv$class)
  
  message("\n--- Leave-one-out classification accuracy per source -------------")
  acc_per_group <- diag(conf) / rowSums(conf) * 100
  for (grp in names(acc_per_group)) {
    n_correct <- diag(conf)[grp]
    n_total   <- rowSums(conf)[grp]
    message(sprintf("  %-20s  %d/%d  (%.0f%%)", grp, n_correct, n_total, acc_per_group[grp]))
  }
  overall <- sum(diag(conf)) / sum(conf) * 100
  message(sprintf("  %-20s  %d/%d  (%.0f%%)", "Overall", sum(diag(conf)), sum(conf), overall))
    
  # 3. Confusion matrix — shows which sources are confused with which
  # Only print if at least one off-diagonal cell is non-zero
  if (any(conf[row(conf) != col(conf)] > 0)) {
    message("\n--- Confusion matrix (rows = actual, cols = predicted) -----------")
    message(paste(utils::capture.output(print(conf)), collapse = "\n"))
    # Highlight the most confused pairs
    conf_offdiag <- conf
    diag(conf_offdiag) <- 0L
    if (sum(conf_offdiag) > 0) {
      worst_idx <- which(conf_offdiag == max(conf_offdiag), arr.ind = TRUE)[1, ]
      message(sprintf("\n  Most confused pair: '%s' misclassified as '%s' (%d time%s)",
                      rownames(conf)[worst_idx["Actual"]],
                      colnames(conf)[worst_idx["Predicted"]],
                      max(conf_offdiag),
                      if (max(conf_offdiag) == 1) "" else "s"))
    }
  } else {
    message("\n  Perfect classification - all sources fully separated.")
  }
  message("------------------------------------------------------------------\n")
  
  invisible(list(lda = data.lda, confusion = conf, prop_var = prop_var))
}

# ── Aliases ──────────────────────────────────────────────────────────────────
#' @rdname LDA_Chart
#' @export
LDA_plot <- LDA_Chart

#' @rdname LDA_Chart
#' @export
LDAPlot <- LDA_Chart