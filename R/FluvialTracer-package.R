#' An integrative tool for the sediment fingerprinting community
#' 
#' @description
#' FluvialTracer is designed to support sediment fingerprinting workflows in river and catchment systems, from tracer screening through source unmixing models. It is built for researchers working on soil erosion, land degradation, and sediment connectivity in fluvial environments, with particular attention to catchments where source signatures are affected by land-use change. All statistical analysis and models are integrated within a single framework, streamlining the workflow for researchers and offering multiple model applications for comparing possible results.
#' @author
#' **Authors:**
#' * Margarita Villalonga (mvillalonga@ipe.csic.es)
#' * Esther Cepeda
#' * Maarten Wynants
#' * Ivan Lizaga - Project and thesis supervisor (ivan.lizaga@ipe.csic.es)
#' 
#' **Contributors:**
#' * Samuel Bodé, Pascal Boeckx, William Blake, Jan Vanaverbeke, Ulrike Braeckman, Aiming Liang, Borja Latorre, Yuri Jacques Agra Bezerra da Silva, Leticia Gaspar, Jorge Crespo, Ana Navas, Mohammadi Maziar, Vince Q Vu
#'
#' @md
#'
#' @section Tool development:
#'
#' FluvialTracer is being developed as a community-driven tool, created by 
#' and for the sediment fingerprinting research community. Its main goal is 
#' to bring together the different methods currently used for sediment 
#' fingerprinting within a single, integrated unmixing framework, enabling 
#' researchers to easily apply, explore, and compare alternative approaches.
#'
#' @keywords internal
#' 
#' @import ggplot2
#' @importFrom stats lm.fit coef rt rnorm runif prcomp predict cor hclust as.dist cutree as.formula kruskal.test aggregate sd quantile qchisq var dist complete.cases median
#' @importFrom utils txtProgressBar setTxtProgressBar head tail
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom graphics par title
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select filter group_by ungroup across everything rename_with select_if rowwise
#' @importFrom tidyr unnest
#' @importFrom tidyselect starts_with ends_with
#' @importFrom reshape2 melt
#' @importFrom Ternary TernaryPlot TernaryPoints
#' @importFrom klaR greedy.wilks
#' @importFrom plotly plot_ly add_markers
#'
#' @examples
#' \donttest{
#' # Load the 'FluvialTracer' package to access its functions.
#' library('FluvialTracer')
#' 
#' # Load the example dataset for a 3-source mixing problem.
#' data <- fluvialdata_example
#' 
#' # 1. Visualize the distribution of each tracer across sources and mixture.
#' Graph_boxplots(data, max_per_page = 6)
#' 
#' # 2. Check how well sources separate using Linear Discriminant Analysis.
#' LDA_Chart(data)
#' 
#' # 3. Evaluate tracer conservativeness with the Conservativeness Index.
#' ci_results <- CI_Method(data, points = 2000)
#' 
#' # 4. Estimate source contributions to the mixture (frequentist unmixing).
#' rs <- unmixing(data, iter = "short")
#' plot_apportionment(rs, violin = TRUE)
#' }
"_PACKAGE"