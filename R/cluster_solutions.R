#' Cluster CTS solutions by source contribution similarity
#'
#' Takes the output of \code{\link{CTS_n_sources}} and groups the tracer
#' combinations into clusters of similar source contributions (w1, w2, ...).
#' This helps identify the dominant apportionment solutions and how many
#' individual tracers support each one.
#'
#' Hierarchical clustering with Ward.D2 linkage is used because contribution
#' clusters are rarely spherical or equal-sized. The optimal number of clusters
#' is selected automatically via the average silhouette width criterion, but
#' can be overridden with the \code{k} argument.
#'
#' Visualisation adapts to the number of sources:
#' \itemize{
#'   \item 2 sources: single 2D panel (w1 vs w2).
#'   \item 3 sources: two 2D panels (w1~w2, w1~w3) plus optional interactive
#'     3D plotly scatter.
#'   \item 4+ sources: all pairwise 2D panels as a grid. Interactive 3D is not
#'     available (dimensionality reduction would obscure interpretation).
#' }
#'
#' The physically feasible region (0-100% contribution on both axes) is
#' shaded white; everything outside it (negative or >100% contributions,
#' i.e. not physically possible) is shaded grey.
#'
#' @param cts_result Data frame returned by \code{CTS_n_sources} or any of the
#'   \code{CTS_*_sources} functions.
#' @param top_n Integer. Number of top rows to use for clustering. Default 50.
#' @param k Integer or \code{NULL}. Number of clusters. If \code{NULL}
#'   (default), chosen automatically via silhouette width.
#' @param max_k Integer. Maximum k to evaluate when \code{k = NULL}. Default 6.
#' @param interactive Boolean. If \code{TRUE} and n_sources == 3, produces an
#'   interactive 3D plotly scatter. Ignored for 2 or 4+ sources. Default
#'   \code{FALSE}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{summary}}{Data frame: one row per cluster with mean/SD of each
#'     contribution, n_pairs, n_unique_tracers, avg_silhouette.}
#'   \item{\code{tracers}}{Named list. Each element has \code{$pairs} (character
#'     vector) and \code{$individual_tracers} (data frame ranked by frequency).}
#'   \item{\code{plot}}{ggplot2 or plotly visualisation.}
#' }
#'
#' @references
#' Wynants, M., Verboom, C., Lizaga, I., de Boer, H., van der Perk, M.,
#' Bennett, W., Boeckx, P., & Brooks, A. Sediment and phosphorus source
#' hysteresis and tracing during an extreme flood in a subtropical
#' Australian river.
#'
#' @export
#' 
cluster_solutions <- function(cts_result,
                              top_n       = 50L,
                              k           = NULL,
                              max_k       = 6L,
                              interactive = FALSE,
                              seed        = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # ── 1. Identify w columns ─────────────────────────────────────────────────
  w_cols <- grep("^w[0-9]+$", colnames(cts_result), value = TRUE)
  n_src  <- length(w_cols)
  if (n_src < 2L) stop("No 'w1, w2, ...' columns found in cts_result.")
  
  # ── 2. Subset top_n ───────────────────────────────────────────────────────
  top_n <- min(top_n, nrow(cts_result))
  dat   <- cts_result[seq_len(top_n), ]
  w_mat <- as.matrix(dat[, w_cols])
  if (nrow(w_mat) < 3L)
    stop("Need at least 3 rows. Increase top_n or check cts_result.")
  
  # ── 3. Hierarchical clustering (Ward.D2) ──────────────────────────────────
  dist_w <- dist(w_mat, method = "euclidean")
  hc     <- hclust(dist_w, method = "ward.D2")
  
  # ── 4. Choose k ───────────────────────────────────────────────────────────
  if (is.null(k)) {
    max_k   <- min(max_k, nrow(w_mat) - 1L)
    sil_avg <- vapply(2L:max_k, function(ki) {
      cl  <- cutree(hc, k = ki)
      sil <- cluster::silhouette(cl, dist_w)
      mean(sil[, "sil_width"])
    }, numeric(1L))
    k_opt <- which.max(sil_avg) + 1L
    message(sprintf(
      "Optimal k = %d (avg silhouette = %.3f). Use k= to override.",
      k_opt, max(sil_avg)))
  } else {
    k_opt <- as.integer(k)
  }
  
  clusters  <- cutree(hc, k = k_opt)
  sil_obj   <- cluster::silhouette(clusters, dist_w)
  sil_width <- sil_obj[, "sil_width"]
  
  # ── 5. Helper: split "Ba Sr" → c("Ba","Sr") ───────────────────────────────
  .split_tracers <- function(id_vec) unlist(strsplit(trimws(id_vec), "\\s+"))
  
  # ── 6. Summary table ──────────────────────────────────────────────────────
  summary_rows <- lapply(seq_len(k_opt), function(ci) {
    idx    <- which(clusters == ci)
    w_sub  <- w_mat[idx, , drop = FALSE]
    w_mean <- colMeans(w_sub)
    w_sd   <- if (nrow(w_sub) > 1L) apply(w_sub, 2L, sd) else rep(0, n_src)
    indiv  <- .split_tracers(dat$id[idx])
    
    row <- data.frame(cluster          = ci,
                      n_pairs          = length(idx),
                      n_unique_tracers = length(unique(indiv)),
                      avg_silhouette   = round(mean(sil_width[idx]), 3L),
                      stringsAsFactors = FALSE)
    for (i in seq_len(n_src)) {
      row[[paste0("mean_w", i)]] <- round(w_mean[i], 4L)
      row[[paste0("sd_w",   i)]] <- round(w_sd[i],   4L)
    }
    row
  })
  summary_df <- do.call(rbind, summary_rows)
  
  # ── 7. Tracers per cluster ─────────────────────────────────────────────────
  tracer_list <- lapply(seq_len(k_opt), function(ci) {
    idx   <- which(clusters == ci)
    pairs <- dat$id[idx]
    indiv <- .split_tracers(pairs)
    freq  <- sort(table(indiv), decreasing = TRUE)
    list(
      pairs              = pairs,
      individual_tracers = data.frame(tracer           = names(freq),
                                      supporting_pairs = as.integer(freq),
                                      stringsAsFactors = FALSE)
    )
  })
  names(tracer_list) <- paste0("cluster_", seq_len(k_opt))
  
  # ── 8. Colour palette ─────────────────────────────────────────────────────
  pal <- scales::hue_pal()(k_opt)
  
  # ── 9. Build plot_df and per-cluster centroid label_df ────────────────────
  plot_df <- as.data.frame(w_mat)
  colnames(plot_df) <- w_cols
  plot_df$cluster <- factor(clusters)
  plot_df$tracer  <- dat$id
  
  # Generic centroid builder for any pair of w columns
  .make_label_df <- function(xvar, yvar) {
    do.call(rbind, lapply(seq_len(k_opt), function(ci) {
      sub <- plot_df[plot_df$cluster == ci, ]
      data.frame(cluster = factor(ci),
                 x       = mean(sub[[xvar]]),
                 y       = mean(sub[[yvar]]),
                 label   = sprintf("C%d (n=%d)", ci, nrow(sub)),
                 stringsAsFactors = FALSE)
    }))
  }
  
  # Generic 2D panel builder
  .panel <- function(xvar, yvar, show_legend = TRUE) {
    ldf <- .make_label_df(xvar, yvar)
    
    # Zoom the view to the actual data range (with a margin), so the
    # background shading rectangles (which span -Inf/Inf and 0/100) don't
    # force the whole 0-100 square into view and shrink the real points.
    x_rng  <- range(plot_df[[xvar]], na.rm = TRUE)
    y_rng  <- range(plot_df[[yvar]], na.rm = TRUE)
    x_pad  <- diff(x_rng) * 0.15; if (x_pad == 0) x_pad <- 1
    y_pad  <- diff(y_rng) * 0.15; if (y_pad == 0) y_pad <- 1
    
    # Gridline positions, computed over the zoomed view — drawn as explicit
    # layers (below) so they sit ON TOP of the shading rectangles instead of
    # being hidden underneath them (theme_bw()'s own grid always renders
    # first, behind every geom layer, so it can't show through an opaque
    # rect layer added afterwards).
    x_breaks <- scales::extended_breaks()(c(x_rng[1] - x_pad, x_rng[2] + x_pad))
    y_breaks <- scales::extended_breaks()(c(y_rng[1] - y_pad, y_rng[2] + y_pad))
    
    ggplot2::ggplot(plot_df,
                    ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]],
                                 colour = cluster, fill = cluster)) +
      # ── Shade the physically infeasible region grey, and mark the valid
      # 0-100% contribution square with a white fill AND a clear dashed
      # border, so the boundary between feasible and infeasible is visible.
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                        fill = "grey70") +
      ggplot2::annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                        fill = "white", colour = "black",
                        linetype = "dashed", linewidth = 0.6) +
      # ── Gridlines redrawn on top of the shading. grey60 was chosen
      # because it stays visible against both the white (feasible) and
      # grey70 (infeasible) backgrounds.
      ggplot2::geom_hline(yintercept = y_breaks, colour = "grey60",
                          linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = x_breaks, colour = "grey60",
                          linewidth = 0.3) +
      ggplot2::geom_point(size = 2.5, alpha = 0.75) +
      ggplot2::geom_text(ggplot2::aes(label = tracer),
                         size = 2.6, hjust = -0.1, vjust = 0.4,
                         show.legend = FALSE) +
      stat_ellipse(type = "t", linewidth = 0.8,
                   geom = "polygon", alpha = 0.12, level = 0.85) +
      stat_ellipse(type = "t", linewidth = 0.8, alpha = 0.8, level = 0.85) +
      ggplot2::geom_label(data = ldf,
                          ggplot2::aes(x = x, y = y, label = label,
                                       colour = cluster),
                          fill = "white", size = 3,
                          label.padding = ggplot2::unit(0.12, "lines"),
                          show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = pal) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = xvar, y = yvar,
                    colour = "Cluster", fill = "Cluster") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = if (show_legend) "right" else "none",
                     panel.grid = ggplot2::element_blank()) +
      ggplot2::coord_cartesian(xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
                               ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
                               expand = FALSE)
  }
  
  main_title <- sprintf(
    "CTS solution clusters  |  top %d tracers  |  k = %d", top_n, k_opt)
  
  # ── 10a. 3D interactive (only for exactly 3 sources) ──────────────────────
  if (interactive && n_src == 3L) {
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Package 'plotly' is required for interactive = TRUE.")
    
    p <- plotly::plot_ly(
      data   = plot_df,
      x      = ~w1, y = ~w2, z = ~w3,
      color  = ~cluster,
      colors = pal,
      text   = ~paste0(tracer, "<br>Cluster ", cluster),
      type   = "scatter3d",
      mode   = "markers",
      marker = list(size = 5, opacity = 0.8)
    ) |>
      plotly::layout(
        title = list(text = main_title),
        scene = list(xaxis = list(title = "w1 (source 1)"),
                     yaxis = list(title = "w2 (source 2)"),
                     zaxis = list(title = "w3 (source 3)"))
      )
    
    if (n_src > 3L) {
      message("interactive = TRUE is only available for 3 sources. ",
              "Switching to static 2D panels.")
      interactive <- FALSE
    }
  }
  
  # ── 10b. Static 2D panels ─────────────────────────────────────────────────
  if (!interactive || n_src != 3L) {
    
    if (n_src == 2L) {
      # Single panel
      p <- .panel("w1", "w2") +
        ggplot2::ggtitle(main_title) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
    } else if (n_src == 3L) {
      # Two panels: w1~w2 and w1~w3
      p <- gridExtra::grid.arrange(
        .panel("w1", "w2"),
        .panel("w1", "w3", show_legend = FALSE),
        ncol = 2L,
        top  = grid::textGrob(main_title,
                              gp = grid::gpar(fontsize = 12L,
                                              fontface = "bold")))
      
    } else {
      # 4+ sources: all pairwise combinations
      pairs     <- utils::combn(w_cols, 2L, simplify = FALSE)
      n_pairs   <- length(pairs)
      n_col_plt <- min(3L, n_pairs)
      panels    <- lapply(seq_along(pairs), function(i) {
        show_leg <- (i == 1L)
        .panel(pairs[[i]][1L], pairs[[i]][2L], show_legend = show_leg)
      })
      p <- gridExtra::grid.arrange(
        grobs = panels,
        ncol  = n_col_plt,
        top   = grid::textGrob(main_title,
                               gp = grid::gpar(fontsize = 12L,
                                               fontface = "bold")))
    }
  }
  
  invisible(list(summary = summary_df,
                 tracers = tracer_list,
                 plot    = p))
}