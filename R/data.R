#' Default color palette for FluvialTracer plots
#'
#' A named vector of hex color codes used as the default palette across
#' FluvialTracer's plotting functions (e.g. \code{\link{plot_apportionment}},
#' \code{\link{LDA_Chart}}). Designed to remain visually distinct for up to
#' five sources plus the mixture.
#'
#' @format A named character vector of length 20.
#'
#' @source Custom palette defined by the FluvialTracer authors.
"fluvial_palette"


#' Example dataset with 3 sources and one mixture
#'
#' Virtual dataset in the standard FluvialTracer format, used in the
#' package examples.
#'
#' @format A data frame with one row per sample: sample ID, source/mixture
#'   group, and one column per tracer.
#'
#' @examples
#' head(fluvialdata_example)
"fluvialdata_example"