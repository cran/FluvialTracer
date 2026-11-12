#' Input sediment sources
#'
#' Selects and extracts the source samples from the dataset, computing mean,
#' SD and n when the input is in raw format. When a different format dataset 
#' is detected, column names are normalised to \code{Mean_}/\code{SD_} in the 
#' output so that all downstream functions work correctly.
#'
#' @param data Data frame containing source and mixture data.
#' @param na.omit Boolean. If \code{TRUE} (default), \code{NA} values are
#'   omitted when computing the mean and SD.
#'
#' @return A data frame with one row per source, containing the source
#'   identifier (\code{id}), the mean of each tracer (\code{Mean_} columns),
#'   the standard deviation of each tracer (\code{SD_} columns) and the
#'   number of samples (\code{n}).
#'
#' @export
#' 
inputSource <- function(data, na.omit = TRUE) {

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

  # Returns TRUE when prefixes are lowercase
  is_legacy_format <- function(d) {
    any(grepl("^mean_", colnames(d))) || any(grepl("^sd_", colnames(d)))
  }

  # Normalises legacy mean_/sd_ prefixes to OceanTracer Mean_/SD_
  normalise_prefixes <- function(d) {
    colnames(d) <- gsub("^mean_", "Mean_", colnames(d))
    colnames(d) <- gsub("^sd_",   "SD_",   colnames(d))
    d
  }

  # ── Averaged format ─────────────────────────────────────────────────────────
  if (is_averaged_internal(data)) {

    if (is_legacy_format(data)) {
      # message("Note: lowercase 'mean_'/'sd_' column prefixes detected. ",
      #         "Column names have been normalised to 'Mean_'/'SD_' for compatibility with OceanTracer.")
      data <- normalise_prefixes(data)
    }

    id      <- head(data[, 2], n = -1)
    sources <- cbind(id, data[-nrow(data), 3:ncol(data)])
    return(sources)
  }

  # ── Raw format: aggregate to Mean_ / SD_ / n per source group ───────────────
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))

  # Exclude the last group (assumed to be the mixture)
  sources2 <- data[data[, 2] != levels(data[, 2])[nlevels(data[, 2])], ]
  sources2 <- sources2[order(sources2[, 2]), ]
  s_groups <- droplevels(sources2[, 2])

  data_mean <- stats::aggregate(sources2[, 3:ncol(sources2)],
                                list(s_groups), mean, na.rm = na.omit)
  data_mean[[1]] <- NULL
  colnames(data_mean) <- paste0("Mean_", colnames(data_mean))

  data_sd <- stats::aggregate(sources2[, 3:ncol(sources2)],
                              list(s_groups), stats::sd, na.rm = na.omit)
  data_sd <- as.data.frame(data_sd[, -1])
  colnames(data_sd) <- paste0("SD_", colnames(data_sd))

  n_counts      <- as.data.frame(table(sources2[, 2]))
  n_counts[[1]] <- NULL
  n_counts      <- as.data.frame(n_counts[-nrow(n_counts), ])
  colnames(n_counts) <- "n"

  id      <- head(levels(data[, 2]), n = -1)
  sources <- cbind(id, data_mean, data_sd, n_counts)

  return(sources)
}
