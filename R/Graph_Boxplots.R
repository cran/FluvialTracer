#' Box and whiskers plot
#'
#' Visualizes the distribution of continuous tracer variables across sources and 
#' mixtures. It compactly displays five summary statistics (the median, two hinges, 
#' and two whiskers), alongside individual outliers.
#'
#' @param data Data frame. Contains source and mixture data in the standard FluvialTracer format.
#' @param tracers Numeric vector. Indices of the tracers to plot (column 1 refers to the first tracer variable). Default 1:(ncol(data) - 2).
#' @param ncol Integer. Number of plots per row in the grid layout. Default 3.
#' @param colors Character vector. Optional hex codes or names for custom boxplot colors, one
#'   per group (sources + mixture). If NULL (default), colors are pulled automatically from
#'   fluvial_palette: the first N-1 Source colors, plus one Mixture color for the last group.
#' @param seed Integer. Random seed for reproducible jittering or sampling. Default NULL.
#' @param max_per_page Integer. Maximum number of tracer panels to draw per page/device call.
#'   If the number of requested tracers exceeds this, the plots are split across multiple
#'   pages, drawn sequentially (each replacing the previous on-screen, or as successive pages
#'   if called inside a pdf()/dev device). Default 16.
#' @param label_size Numeric. Font size for the n= and mean annotation labels inside each
#'   panel. Lower this if labels overlap when many panels are shown. Default 3.
#'
#' @return Invisibly returns a list of \code{gtable} objects, one per page.
#'   The plots are also drawn on the current graphics device.
#'
#' @export
#' 
Graph_boxplots <- function(data, tracers = 1:(ncol(data) - 2), ncol = 3, colors = NULL,
                           seed = NULL, max_per_page = 16, label_size = 3) {
  
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
  
  # ── Plotting ────────────────────────────────────────────────────────────────
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  group_col <- colnames(data)[2]
  
  # Default colors: N-1 Source colors + 1 Mixture color, matching the group
  # levels of column 2 (sources first, mixture assumed to be the last level).
  if (is.null(colors)) {
    grp_levels <- levels(data[, 2])
    n_groups   <- length(grp_levels)
    colors <- c(.default_source_colors(n_groups - 1),
                .default_source_colors(1, mixture = TRUE))
    names(colors) <- grp_levels
  }
  
  dat_plot  <- reshape2::melt(data, id = c(1, 2))
  
  g_legend <- function(a.gplot) {
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) > 0) tmp$grobs[[leg]] else NULL
  }
  
  blank <- grid::grid.rect(gp = grid::gpar(col = "white"))
  
  # ── Build a single page of panels for a given subset of tracer indices ─────
  .build_page <- function(tracer_subset) {
    
    funclist <- list()
    for (i in seq_along(tracer_subset)) {
      idx <- tracer_subset[i]
      funclist[[i]] <- local({
        idx_local <- idx
        function(x) {
          y <- dat_plot[dat_plot$variable == unique(dat_plot$variable)[idx_local], 4]
          data.frame(y = min(y) - 0.05 * (max(y) - min(y)), label = paste0('n=', length(x)))
        }
      })
    }
    
    funclist2 <- list()
    for (i in seq_along(tracer_subset)) {
      idx <- tracer_subset[i]
      funclist2[[i]] <- local({
        idx_local <- idx
        function(x) {
          y <- dat_plot[dat_plot$variable == unique(dat_plot$variable)[idx_local], 4]
          data.frame(y = median(x) + 0.03 * (max(y) - min(y)), label = round(mean(x), 1))
        }
      })
    }
    
    glist <- list()
    j <- 1
    
    for (i in seq_along(tracer_subset)) {
      idx <- tracer_subset[i]
      p <- ggplot2::ggplot(
        dat_plot[dat_plot$variable == unique(dat_plot$variable)[idx], ],
        ggplot2::aes(x = .data[[group_col]], y = value, color = .data[[group_col]])
      ) +
        ggplot2::geom_boxplot() +
        ggplot2::ggtitle(unique(dat_plot$variable)[idx]) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::xlab("") + ggplot2::ylab("") +
        ggplot2::stat_summary(fun.data = funclist[[i]],  geom = "text",
                              fontface = "bold", color = "black", alpha = 0.5, size = label_size) +
        ggplot2::stat_summary(fun.data = funclist2[[i]], geom = "text",
                              fontface = "bold", color = "black", alpha = 0.5, size = label_size) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                       legend.position = "none")
      
      p <- p + ggplot2::scale_color_manual(values = colors)
      
      glist[[j]] <- p
      
      if (i == 1) {
        mylegend <- g_legend(p)
        if (!is.null(mylegend)) { j <- j + 1; glist[[j]] <- mylegend }
      } else if (i %% ncol == 0) {
        j <- j + 1; glist[[j]] <- blank
      }
      j <- j + 1
    }
    
    gridExtra::arrangeGrob(grobs = glist, ncol = ncol + 1,
                           widths = c(rep(10, each = ncol), 3))
  }
  
  # ── Split tracers into pages and draw each ──────────────────────────────────
  pages <- split(tracers, ceiling(seq_along(tracers) / max_per_page))
  
  plots_out <- list()
  for (pg in seq_along(pages)) {
    p <- .build_page(pages[[pg]])
    grid::grid.newpage()
    grid::grid.draw(p)
    plots_out[[pg]] <- p
  }
  
  invisible(plots_out)
}

# ── Aliases ──────────────────────────────────────────────────────────────────
#' @rdname Graph_boxplots
#' @export
boxPlot <- Graph_boxplots

#' @rdname Graph_boxplots
#' @export
box_plot <- Graph_boxplots