# Internal helper — not exported
# Computes the normalized CTS error for each tracer given pre-processed
# source and mixture matrices and a solution weight vector.
# Works for any number of sources via a dot product.
.cts_core <- function(source, mixture, solution) {

  source  <- data.matrix(source[-1L])
  mixture <- data.matrix(mixture[-1L])

  cols <- (ncol(source) - 1L) / 2L

  # Normalize each tracer by the range of its source confidence interval
  for (col in seq_len(cols)) {
    mx <- max(source[, col] + source[, cols + col])
    mn <- min(source[, col] - source[, cols + col])
    rng <- mx - mn
    if (rng == 0) rng <- 1   # avoid division by zero for depleted tracers
    source[, col]        <- (source[, col] - mn) / rng
    source[, cols + col] <- source[, cols + col] / rng
    mixture[, col]       <- (mixture[, col] - mn) / rng
  }

  tracer <- colnames(source)[seq_len(cols)]

  # Dot product: predicted = sum_k(source[k, col] * solution[k])
  err <- vapply(seq_len(cols), function(col) {
    predicted <- sum(source[seq_len(length(solution)), col] * solution)
    abs(predicted - mixture[1L, col])
  }, numeric(1L))

  data.frame(tracer = tracer, consistency_error = err,
             stringsAsFactors = FALSE)
}


#' Consistency error of a tracer selection for an apportionment solution
#'
#' Assesses the mathematical consistency of a tracer selection by computing
#' the normalized error between the predicted and observed tracer concentrations
#' in the mixture. A low error (typically below 0.05) for all tracers indicates
#' a consistent selection. High errors point to non-conservative tracers or
#' multiple partial solutions in the dataset.
#'
#' The error for each tracer is the absolute difference between the observed
#' mixture concentration and the virtual mixture reconstructed from the
#' apportionment solution, after normalizing by the range of the source
#' confidence intervals.
#'
#' @param data A data frame containing the characteristics of sediment sources
#'   and one mixture, in the standard fingerprinting format.
#' @param solution Either a numeric vector of apportionment weights (one per
#'   source, in the same order as they appear in \code{data}), or a single-row
#'   data frame with columns named \code{w1}, \code{w2}, ... as returned by
#'   the unmixing functions.
#'
#' @return A data frame with columns \code{tracer} and \code{consistency_error},
#'   one row per tracer.
#'
#' @references
#' Latorre, B., Lizaga, I., Gaspar, L., & Navas, A. (2021). A novel method for
#' analysing consistency and unravelling multiple solutions in sediment
#' fingerprinting. \emph{Science of The Total Environment}, \emph{789}, 147804.
#' 
#' @export
#' 
consistency_error <- function(data, solution) {

  source   <- inputSource(data)
  mixture  <- inputSample(data)
  source_n <- nrow(source)

  # Strip Mean_/mean_ prefix so column indexing inside .cts_core works
  # regardless of whether the user passed averaged or raw data
  colnames(source)  <- gsub("^[Mm]ean_", "", colnames(source))
  colnames(mixture) <- gsub("^[Mm]ean_", "", colnames(mixture))

  # Accept either a named data frame (unmix output) or a plain numeric vector
  if (is.data.frame(solution)) {
    w_cols <- grep("^w[0-9]+$", colnames(solution), value = TRUE)
    if (length(w_cols) == 0L)
      stop("'solution' data frame must contain columns named w1, w2, ...")
    sol <- as.numeric(solution[1L, w_cols])
  } else if (is.numeric(solution)) {
    sol <- solution
  } else {
    stop("'solution' must be a numeric vector or a data frame with w1, w2, ... columns.")
  }

  if (length(sol) != source_n)
    stop("Length of solution (", length(sol), ") does not match number of sources (", source_n, ").")

  .cts_core(source, mixture, sol)
}


# --- Backwards Compatibility Aliases ---
#' @rdname consistency_error
#' @export
CTS_error <- consistency_error

#' @rdname consistency_error
#' @export
cts_error <- consistency_error