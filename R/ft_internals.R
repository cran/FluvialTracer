# Internal helper - not exported
#
# Returns the first `n` colors from fluvial_palette's Source (or Mixture)
# series, unnamed. Used as the automatic default in plotting functions when
# the user does not supply their own `colors` argument.
#
# @param n Integer. Number of colors needed.
# @param mixture Logical. If TRUE, pulls from the Mixture series instead of
#   the Source series. Default FALSE.
.default_source_colors <- function(n, mixture = FALSE) {
  prefix <- if (mixture) "Mixture" else "Source"
  unname(fluvial_palette[paste0(prefix, seq_len(n))])
}