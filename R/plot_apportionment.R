#' Plot apportionment results
#'
#' Displays the relative contribution of sediment sources to each mixture as
#' density or violin plots. Designed to take the output of the unmixing
#' functions directly.
#'
#' @param data A data frame, typically the output from \code{unmix()}, containing
#'   the relative source contributions. The first two columns must be the mixture
#'   identifier (e.g. \code{id}) and the goodness-of-fit (\code{gof}), followed
#'   by one column per source. Only the source columns are plotted.
#' @param violin Logical. If \code{TRUE}, violin + boxplot charts are used
#'   instead of density plots. Default \code{FALSE}.
#' @param bounds A numeric vector of length 2 specifying the lower and upper
#'   display bounds (the camera window). Default \code{c(0.0, 1.0)}.
#' @param scaled Logical. If \code{TRUE}, density curves are area-scaled to a
#'   common maximum. Only used when \code{violin = FALSE}. Default \code{TRUE}.
#' @param y_high Maximum value for the y-axis. Only used when
#'   \code{violin = FALSE} and \code{scaled = FALSE}. Default \code{1}.
#' @param colors An optional named character vector of colors, one per source.
#'   If \code{NULL} (default), the first N colors of \code{fluvial_palette}'s
#'   Source series are used automatically, where N is the number of sources
#'   detected in \code{data}.
#' @param vlines An optional named numeric vector specifying true/expected values 
#'   for specific sources. Draws a dashed vertical line to allow visual validation 
#'   of the model's accuracy (e.g., \code{c("SourceA" = 0.45)}).
#' @param title Optional character string to add a main title to the plot.
#' @param bw The smoothing bandwidth to be used. Can be a numeric value (e.g., 0.06) 
#'   or a character string like "nrd" (the default automatic method).
#'   
#' @return A \code{ggplot2} object.
#'
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_density geom_violin geom_boxplot geom_vline stat_summary coord_cartesian ylim scale_fill_manual scale_color_manual labs ggtitle after_stat
#'
#' @export
#' 
plot_apportionment <- function(data, violin = FALSE, bounds = c(0.0, 1.0), 
                               scaled = TRUE, y_high = 1, colors = NULL, 
                               vlines = NULL, title = NULL, bw = "nrd") {
  
  # Melt the data (assuming Col 1 is ID and Col 2 is GOF)
  data_plots <- reshape2::melt(data, id.vars = 1:2)
  
  # Default colors: pull the first N "Source" colors from fluvial_palette,
  # where N is the number of source columns detected in the data.
  if (is.null(colors)) {
    n_sources <- ncol(data) - 2
    colors <- .default_source_colors(n_sources)
    # Name them to match the actual source column names so scale_*_manual
    # matches them correctly regardless of column order/naming.
    names(colors) <- colnames(data)[3:ncol(data)]
  }
  
  if (!violin) {
    if (scaled) {
      p <- ggplot(data_plots, aes(x = value)) +
        # Removed the 'bounds' argument from geom_density so it smooths naturally
        geom_density(aes(group = variable, colour = variable, fill = variable,
                         y = after_stat(scaled)),
                     n = 2000L, bw = bw, alpha = 0.35)
    } else {
      p <- ggplot(data_plots, aes(x = value)) +
        geom_density(aes(group = variable, colour = variable, fill = variable),
                     n = 2000L, bw = bw, alpha = 0.35) +
        ylim(0, y_high)
    }
    
    # coord_cartesian acts as the camera. It crops the view without distorting the density.
    p <- p + coord_cartesian(xlim = bounds, expand = TRUE) +
      labs(x = "Proportional Contribution", y = ifelse(scaled, "Scaled Posterior Density", "Density"))
    
    # Add student's true-value vertical lines
    if (!is.null(vlines)) {
      for (varname in names(vlines)) {
        val <- vlines[[varname]]
        # Match the line color to the density fill color if colors are provided
        col <- if (!is.null(colors) && varname %in% names(colors)) colors[[varname]] else "black"
        
        if (!is.na(val)) {
          p <- p + geom_vline(xintercept = val, color = col, linetype = "dashed", linewidth = 0.8)
        }
      }
    }
    
  } else {
    p <- ggplot(data_plots, aes(x = variable, y = value, colour = variable)) +
      geom_violin(bw = bw, trim = FALSE, alpha = 0.65) +
      geom_boxplot(width = 0.2, alpha = 0.65, colour = "black",
                   aes(fill = variable)) +
      stat_summary(fun = mean, geom = "point", size = 2.5,
                   color = "red", shape = 8) +
      labs(x = "Source", y = "Apportionment") +
      coord_cartesian(ylim = bounds, expand = TRUE)
  }
  
  # Apply Colors
  if (!is.null(colors)) {
    p <- p +
      scale_fill_manual(values  = colors) +
      scale_color_manual(values = colors)
  }
  
  # Add Title
  if (!is.null(title)) {
    p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(p)
}


# --- Backwards Compatibility Aliases ---
#' @rdname plot_apportionment
#' @export
plotResults <- plot_apportionment

#' @rdname plot_apportionment
#' @export
plot_results <- plot_apportionment