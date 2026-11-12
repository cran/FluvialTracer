#' Checking the difference between virtual and artificial mixture (for data sets containing only one mixture)
#'
#' A function that calculates and visualizes the percentage deviation between 
#' an experimental artificial mixture and a theoretical virtual mixture calculated from user-defined source proportions.
#' 
#'
#' @param data Data frame. Contains source and mixture data in the standard format.
#' @param proportions Numeric vector. Defines the proportions of each source in the mixture. Required to calculate virtual mixture. 
#' @param sample_col Character. Defines which column of the data contains the sample names.
#' @param source_col Character. Defines which column of the data contains the source names. Indicate the column number in which the sources are defined independently from the sample number. 
#' @param subtitle_extra Character. Custom subtitle for the plot; overrides the automatic one if provided. Default NULL.
#'
#' @return Invisibly returns a \code{gtable} object with the combined plot
#'   (both panels, title, legend and caption). The plot is also drawn on the
#'   current graphics device.
#'
#' @export
checking_diff_single_mix <- function(data, 
                                     proportions = c(1, 0, 0, 0), 
                                     sample_col = 1,   
                                     source_col = 2,   
                                     subtitle_extra = NULL) {
  
  # -- Internal helpers: detect and expand averaged (Mean_/SD_) format --------
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
  
  # -- Pre-processing: expand to raw observations if data is in Mean_/SD_ format
  if (is_averaged_internal(data)) {
    data <- raw_dataset_internal(data)
  }
  
  if (!is.character(source_col)) source_col <- colnames(data)[as.numeric(source_col)]
  if (!is.character(sample_col)) sample_col <- colnames(data)[as.numeric(sample_col)]
  
  source_names   <- head(unique(data[[source_col]]), -1)
  artificial_mix <- tail(unique(data[[source_col]]), 1)   # la última = la mezcla real
  
  # 1. Calculate the raw source means
  source_means <- data |>       
    dplyr::select(-all_of(sample_col)) |> 
    dplyr::filter(.data[[source_col]] %in% source_names) |>  
    dplyr::group_by(.data[[source_col]]) |> 
    dplyr::summarise(dplyr::across(tidyselect::where(is.numeric), \(x) mean(x, na.rm = TRUE)))
  
  # 2 & 3. Define proportions and calculate virtual mixture
  virtual_mix <- tibble::tibble(
    !!source_col := source_names,
    Proportion   = proportions
  ) |>
    dplyr::left_join(source_means, by = source_col) |>
    dplyr::mutate(dplyr::across(-c(dplyr::all_of(source_col), Proportion), as.numeric)) |>
    dplyr::mutate(dplyr::across(-c(dplyr::all_of(source_col), Proportion), ~ .x * Proportion)) |>
    dplyr::summarise(dplyr::across(tidyselect::where(is.numeric), sum)) |>
    dplyr::mutate(!!source_col := "virtual_mix")
  
  # 4. Generate a comparison dataset with the virtual mix and the lab mix
  comparison <- dplyr::bind_rows(
    data |> 
      dplyr::filter(.data[[source_col]] == artificial_mix) |>
      dplyr::select(-all_of(sample_col)),
    
    virtual_mix |> 
      dplyr::mutate(!!source_col := "virtual_mix")
  )
  
  # 5. Calculate the Percentage Difference
  artificial_mix_data <- comparison |>
    dplyr::filter(.data[[source_col]] == artificial_mix) |>
    dplyr::select(tidyselect::where(is.numeric))
  
  virtual_mix_data <- comparison |>
    dplyr::filter(.data[[source_col]] == "virtual_mix") |>
    dplyr::select(tidyselect::where(is.numeric))
  
  perc_diff <- ((artificial_mix_data - virtual_mix_data) / virtual_mix_data) * 100
  
  # Identificar elementos negativos
  neg_tracers <- colnames(artificial_mix_data)[
    as.logical(artificial_mix_data[1, ] < 0 | virtual_mix_data[1, ] < 0)
  ]
  
  # 6. Reshape the data for plotting
  perc_diff_long <- perc_diff |>
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Tracer", values_to = "Difference") |>
    dplyr::mutate(
      Difference_Plot = ifelse(Tracer %in% neg_tracers, NA, Difference),
      Direction = ifelse(Difference_Plot > 0, "Overestimated", "Underestimated"),
      
      Label_HTML = dplyr::case_when(
        Tracer %in% neg_tracers ~ paste0("<span style='color:red'>", Tracer, "*</span>"),
        abs(Difference) < 10    ~ paste0("**", Tracer, "**"), 
        TRUE                    ~ Tracer
      )
    )
  
  #----PLOTS-----
  
  # 7. Create Plot 1: Sorted by Magnitude (The "Nice" one)
  perc_diff_long_p1 <- perc_diff_long |>
    dplyr::mutate(Sort_Val = ifelse(is.na(Difference_Plot), -Inf, Difference_Plot)) |>
    dplyr::mutate(Label_HTML = stats::reorder(Label_HTML, Sort_Val))
  
  p1 <- ggplot2::ggplot(perc_diff_long_p1, ggplot2::aes(x = Label_HTML, y = Difference_Plot, fill = Direction)) +
    ggplot2::geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Sorted by Magnitude",
                  x = "Tracer", y = "Difference (%)") +
    ggplot2::scale_fill_manual(
      values = c("Overestimated" = "indianred", "Underestimated" = "steelblue"),
      na.translate = FALSE 
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggtext::element_markdown() 
    )
  
  # 8. Create Plot 2: Sorted by Dataset Column Order
  original_order <- colnames(data)[3:ncol(data)]
  
  html_mapping <- stats::setNames(perc_diff_long$Label_HTML, perc_diff_long$Tracer)
  ordered_html_levels <- html_mapping[original_order]
  
  perc_diff_long_p2 <- perc_diff_long |>
    dplyr::mutate(Label_Orig = factor(Label_HTML, levels = rev(ordered_html_levels)))
  
  p2 <- ggplot2::ggplot(perc_diff_long_p2, ggplot2::aes(x = Label_Orig, y = Difference_Plot, fill = Direction)) +
    ggplot2::geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Dataset Column Order",
                  x = "", y = "Difference (%)") +
    ggplot2::scale_fill_manual(
      values = c("Overestimated" = "indianred", "Underestimated" = "steelblue"),
      na.translate = FALSE
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggtext::element_markdown() 
    )
  
  # Extract the legend from p2 (built with legend before stripping) to place
  # it once, below both panels, so neither panel loses vertical space to it.
  p2_with_legend <- p2 + ggplot2::theme(legend.position = "bottom")
  g_legend <- function(a.gplot) {
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) > 0) tmp$grobs[[leg]] else NULL
  }
  legend_grob <- g_legend(p2_with_legend)
  
  # 9. Combine them side-by-side using gridExtra (avoids ggplot2 4.x / patchwork
  #    S7-dispatch incompatibility seen when combining plots with `+`)
  props_label <- paste(source_names, proportions, sep = " = ", collapse = " | ")
  mix_name    <- data[nrow(data), 1] 
  
  subtitle_txt <- if (!is.null(subtitle_extra)) {
    subtitle_extra
  } else {
    paste0("Percentage Difference (%) between probabilistic and deterministic virtual mixtures  |  ",
           "Artificial mixture: ", mix_name,
           "  |  Proportions: ", props_label)
  }
  
  title_grob <- grid::textGrob(
    "Deviation of Experimental Artificial Mixture from Theoretical Virtual Mixture",
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )
  subtitle_grob <- grid::textGrob(subtitle_txt, gp = grid::gpar(fontsize = 10))
  
  caption_txt <- "* Red elements contained negative source values (instrumental noise) and were excluded from calculation."
  caption_grob <- grid::textGrob(caption_txt, gp = grid::gpar(fontsize = 9, col = "red"))
  
  if (!is.null(legend_grob)) {
    combined_plot <- gridExtra::arrangeGrob(
      title_grob, subtitle_grob,
      gridExtra::arrangeGrob(p1, p2, ncol = 2),
      legend_grob,
      caption_grob,
      ncol = 1,
      heights = c(0.6, 0.5, 10, 0.5, 0.6)
    )
  } else {
    combined_plot <- gridExtra::arrangeGrob(
      title_grob, subtitle_grob,
      gridExtra::arrangeGrob(p1, p2, ncol = 2),
      caption_grob,
      ncol = 1,
      heights = c(0.6, 0.5, 10, 0.6)
    )
  }
  
  grid::grid.newpage()
  grid::grid.draw(combined_plot)
  
  return(invisible(combined_plot))
}