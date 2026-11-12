#' @title Analyze the Marginal Effect (Pull) of Individual Tracers
#'
#' @description Evaluates the directional bias introduced by individual tracers within 
#' Consistent Tracer Selection (CTS) subset models. This function unpacks the unconstrained 
#' source apportionments to reveal if specific tracers are systematically pulling 
#' source contributions into physically impossible spaces (e.g., >1 or <0).
#'
#' @details 
#' Global unmixing models (using all tracers simultaneously) can suffer from "tracer compensation," 
#' where extreme mathematical contradictions between non-conservative tracers cancel each other out, 
#' resulting in a misleadingly plausible central solution. 
#' 
#' This function exposes that illusion. By taking the output of a CTS run (e.g., \code{\link{CTS_n_sources}}), 
#' it isolates every debate a specific tracer participated in. It then calculates the marginal effect 
#' (median, mean, and SD) that the tracer has on the fractional weight of each source.
#' 
#' Visualizations are automatically generated to aid interpretation:
#' \itemize{
#'   \item \strong{Boxplots:} Display the distribution of calculated weights for each source, grouped by tracer. 
#'   Red dashed lines denote the physical boundaries of 0 and 1.
#'   \item \strong{Heatmap:} Displays the median weight pull. Deep red indicates a strong positive pull 
#'   (overestimation), while deep blue indicates a negative pull (pushing the source into negative mass space).
#' }
#' 
#' @param cts_result A data frame containing the output from a CTS run. It must contain an 
#' \code{id} column with space-separated tracer names, and weight columns named \code{w1}, \code{w2}, etc.
#'   
#' @return A named list containing four elements:
#' \describe{
#'   \item{\code{summary}}{A \code{tibble} containing the Mean, Median, and SD of the weights forced by each tracer on each source.}
#'   \item{\code{boxplot}}{A \code{ggplot} object showing the distribution of source weights grouped by tracer.}
#'   \item{\code{heatmap}}{A \code{ggplot} object visualizing the median tracer bias across all sources.}
#'   \item{\code{raw_expanded}}{A \code{tibble} of the raw, pivoted data used for calculations and plotting.}
#' }
#'
#' @importFrom dplyr select group_by summarise rename all_of
#' @importFrom tidyr separate_rows pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline facet_wrap theme_bw theme element_text element_rect labs geom_tile scale_fill_gradient2 theme_minimal
#' 
#' @export
#' 
tracer_pull_analysis <- function(cts_result) {
  
  # 1. Identify the weight columns dynamically (w1, w2, w3, etc.)
  w_cols <- grep("^w[0-9]+$", colnames(cts_result), value = TRUE)
  
  # 2. Reshape the data: Unpack the tracer combinations
  # We separate the space-delimited 'id' column so each tracer gets its own row
  # Then we pivot the weights into a long format for plotting
  df_long <- cts_result %>%
    select(id, all_of(w_cols)) %>%
    separate_rows(id, sep = " ") %>%
    rename(Tracer = id) %>%
    pivot_longer(cols = all_of(w_cols), names_to = "Source", values_to = "Weight")
  
  # 3. Calculate the central tendency of the "pull"
  summary_pull <- df_long %>%
    group_by(Tracer, Source) %>%
    summarise(
      Mean_Weight = mean(Weight, na.rm = TRUE),
      Median_Weight = median(Weight, na.rm = TRUE),
      SD_Weight = sd(Weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 4. Plot 1: Boxplots of Tracer Pull (Shows extreme compensations)
  p_box <- ggplot(df_long, aes(x = Tracer, y = Weight, fill = Tracer)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.3) +
    # Add red dashed lines to show the physical bounds (0 to 1)
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ Source, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.background = element_rect(fill = "grey90")) +
    labs(title = "Tracer Pull: Distribution of Source Contributions",
         subtitle = "Values outside 0-1 (red lines) indicate the tracer is forcing impossible math.",
         x = "Tracer Included in the Model",
         y = "Calculated Source Weight")
  
  # 5. Plot 2: Heatmap of Median Pull (Shows the specific source bias)
  p_heat <- ggplot(summary_pull, aes(x = Source, y = Tracer, fill = Median_Weight)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", 
                         midpoint = 0.5, # Midpoint set to 0.5 (equal mixing)
                         name = "Median\nWeight") +
    theme_minimal() +
    labs(title = "Heatmap of Tracer Bias",
         subtitle = "Red = Pulls strongly towards source. Blue = Pushes into negative space.",
         x = "Sediment Source",
         y = "Tracer")
  
  # Print the plots
  print(p_heat)
  print(p_box)
  
  return(list(summary = summary_pull, boxplot = p_box, heatmap = p_heat, raw_expanded = df_long))
}
