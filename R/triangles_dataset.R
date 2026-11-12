#' Conservative Index evaluation using randomly sampled companion tracers
#'
#' For each tracer of interest, solves a series of minimal determined linear
#' systems to assess how consistently that tracer constrains the source
#' apportionment. In each iteration the tracer under evaluation is kept fixed
#' while the remaining tracers needed to complete the system are drawn at
#' random from the other measured tracers in the dataset (dataset approach).
#' Because the companion tracers are real geochemical measurements, this
#' strategy reflects realistic co-variation patterns among tracers and
#' produces a less conservative — but more ecologically meaningful — estimate
#' of tracer informativeness compared to the virtual approach.
#'
#' This function is called internally by \code{CI_Method()} and is not usually
#' invoked directly by the user.
#'
#' @param sources Data frame of source statistics: ID, one mean column per
#'   tracer, one SD column per tracer, and sample size \code{n}.
#' @param mixtures Data frame of mixture samples with the same tracer columns.
#' @param tracer Integer (0-based). Index of the tracer under evaluation.
#' @param iter Integer. Number of Monte Carlo iterations per mixture sample.
#' @param seed Integer. Random seed for reproducibility. Default \code{69512}.
#' @param pb An optional \code{txtProgressBar} object created externally.
#'   When provided, the function updates this shared bar instead of creating
#'   its own, allowing \code{CI_Method()} to display a single continuous
#'   progress bar across all tracers. Default \code{NULL}.
#' @param pb_offset Integer. Starting offset for the external progress bar.
#'   Default \code{0L}.
#' @param tracer_current Integer. Position of this tracer in the outer loop,
#'   used to print a progress message. Default \code{NULL}.
#' @param tracer_total Integer. Total number of tracers in the outer loop.
#'   Default \code{NULL}.
#'
#' @return A data frame with \code{nrow(mixtures) * iter} rows. Columns are
#'   \code{ID} (mixture index), \code{GOF} (set to 1 for all rows), and one
#'   weight column per source (\code{w.S1}, \code{w.S2}, ...).
#'
#' @seealso \code{\link{CI_Method}}, \code{\link{triangles_virtual}}
#'
#' @export
#' 
triangles_dataset <- function(sources, mixtures, tracer = 0L, iter = 100L, seed = 69512L,
                               pb = NULL, pb_offset = 0L,
                               tracer_current = NULL, tracer_total = NULL) {

  set.seed(seed)

  # ── Input parsing ──────────────────────────────────────────────────────────
  col_names  <- colnames(sources)
  vars       <- ncol(mixtures) - 1L    # number of tracers
  nsource    <- nrow(sources)
  npoints    <- nrow(mixtures)

  # Basic checks (mirror C++ stops)
  if (ncol(sources) != 2L + (ncol(mixtures) - 1L) * 2L)
    stop("Wrong number of columns in source data frame")
  if (!tolower(col_names[1]) %in% c("id"))
    stop("First data frame column must be ID")
  if (tracer < 0L || tracer >= vars)
    stop("Requested tracer is out of range")
  if (!nsource %in% c(3L, 4L))
    stop("triangles_dataset only supports 3 or 4 sources")

  # ── Read source means, SDs, and n ─────────────────────────────────────────
  src_mean <- as.matrix(sources[, 2L:(vars + 1L),               drop = FALSE])
  src_sd   <- as.matrix(sources[, (vars + 2L):(vars * 2L + 1L), drop = FALSE])
  src_n    <- as.integer(sources[, ncol(sources)])

  # ── Read mixture values ────────────────────────────────────────────────────
  mix_vals <- as.matrix(mixtures[, 2L:(vars + 1L), drop = FALSE])

  # ── Source names for output columns ───────────────────────────────────────
  source_ids <- as.character(sources[, 1L])

  # ── Output storage ─────────────────────────────────────────────────────────
  nrows_out <- npoints * iter
  out_id  <- numeric(nrows_out)
  out_gof <- numeric(nrows_out)
  out_w   <- matrix(0.0, nrow = nrows_out, ncol = nsource)

  # ── Helper: perturb a source value with t-distribution (SEM) ──────────────
  correct_sem <- function(avg, dev, n) {
    if (n > 1L)
      return(avg + rt(1L, df = n - 1L) * dev / sqrt(n))
    return(avg)
  }

  # Available tracer indices (0-based, to match C++ logic)
  all_tracers <- seq_len(vars) - 1L

  # ── Progress bar ───────────────────────────────────────────────────────────
  local_pb <- is.null(pb)
  if (local_pb)
    pb <- txtProgressBar(min = 0, max = npoints * iter, style = 3, width = 50, char = "=")

  # ── Tracer counter message ─────────────────────────────────────────────────
  if (!is.null(tracer_current) && !is.null(tracer_total))
    message(sprintf("Tracer %d out of %d", tracer_current, tracer_total))

  row_idx <- 1L
  for (i in seq_len(npoints)) {

    mix_i <- mix_vals[i, ]
    t1 <- tracer + 1L   # convert to 1-based for R matrix indexing

    for (ii in seq_len(iter)) {

      # Perturb sources
      sc <- matrix(0.0, nrow = nsource, ncol = vars)
      for (s in seq_len(nsource)) {
        for (v in seq_len(vars)) {
          sc[s, v] <- correct_sem(src_mean[s, v], src_sd[s, v], src_n[s])
        }
      }

      setTxtProgressBar(pb, pb_offset + (i - 1L) * iter + ii)

      # ── 3 sources ────────────────────────────────────────────────────────
      if (nsource == 3L) {
        s11 <- sc[1L, t1];  s21 <- sc[2L, t1];  s31 <- sc[3L, t1]
        m1  <- mix_i[t1]

        # Pick a random second tracer (≠ tracer, 0-based, then convert)
        other_tracers <- all_tracers[all_tracers != tracer]
        k  <- sample(other_tracers, 1L) + 1L   # 1-based
        s12 <- sc[1L, k];  s22 <- sc[2L, k];  s32 <- sc[3L, k]
        m2  <- mix_i[k]

        det <- -s22*s31 + s12*s31 + s21*s32 - s11*s32 + s11*s22 - s12*s21
        if (det != 0.0) {
          w1 <- (s21*s32 - s22*s31 + m1*(s22 - s32) + m2*(s31 - s21)) / det
          w2 <- (s12*s31 - s11*s32 + m1*(s32 - s12) + m2*(s11 - s31)) / det
          w3 <- (s11*s22 - s12*s21 + m1*(s12 - s22) + m2*(s21 - s11)) / det
          w  <- c(w1, w2, w3)
        } else {
          w <- c(0.0, 0.0, 0.0)
        }
      }

      # ── 4 sources ────────────────────────────────────────────────────────
      else {  # nsource == 4
        s11 <- sc[1L, t1];  s21 <- sc[2L, t1]
        s31 <- sc[3L, t1];  s41 <- sc[4L, t1];  m1 <- mix_i[t1]

        # Pick second random tracer (≠ tracer)
        other_tracers <- all_tracers[all_tracers != tracer]
        k  <- sample(other_tracers, 1L) + 1L
        s12 <- sc[1L, k];  s22 <- sc[2L, k];  s32 <- sc[3L, k];  s42 <- sc[4L, k]
        m2  <- mix_i[k]

        # Pick third random tracer (≠ tracer, ≠ k)
        remaining <- all_tracers[all_tracers != tracer & all_tracers != (k - 1L)]
        l  <- sample(remaining, 1L) + 1L
        s13 <- sc[1L, l];  s23 <- sc[2L, l];  s33 <- sc[3L, l];  s43 <- sc[4L, l]
        m3  <- mix_i[l]

        det <- s21*(s32*s43 - s33*s42) - s11*(s32*s43 - s33*s42) -
               s31*(s22*s43 - s23*s42) + s11*(s22*s43 - s23*s42) +
               s31*(s12*s43 - s13*s42) - s21*(s12*s43 - s13*s42) +
               (s22*s33 - s23*s32)*s41 - (s12*s33 - s13*s32)*s41 +
               (s12*s23 - s13*s22)*s41 - s11*(s22*s33 - s23*s32) +
               s21*(s12*s33 - s13*s32) - (s12*s23 - s13*s22)*s31

        if (det != 0.0) {
          w1 <- (s21*(s32*s43-s33*s42) - s31*(s22*s43-s23*s42) + (s22*s33-s23*s32)*s41
                 + m1*(-s32*s43+s22*s43+s33*s42-s23*s42-s22*s33+s23*s32)
                 + m2*(s31*s43-s21*s43-s33*s41+s23*s41+s21*s33-s23*s31)
                 + m3*(-s31*s42+s21*s42+s32*s41-s22*s41-s21*s32+s22*s31)) / det

          w2 <- (-s11*(s32*s43-s33*s42) + s31*(s12*s43-s13*s42) - (s12*s33-s13*s32)*s41
                 + m1*(s32*s43-s12*s43-s33*s42+s13*s42+s12*s33-s13*s32)
                 + m2*(-s31*s43+s11*s43+s33*s41-s13*s41-s11*s33+s13*s31)
                 + m3*(s31*s42-s11*s42-s32*s41+s12*s41+s11*s32-s12*s31)) / det

          w3 <- (s11*(s22*s43-s23*s42) - s21*(s12*s43-s13*s42) + (s12*s23-s13*s22)*s41
                 + m1*(-s22*s43+s12*s43+s23*s42-s13*s42-s12*s23+s13*s22)
                 + m2*(s21*s43-s11*s43-s23*s41+s13*s41+s11*s23-s13*s21)
                 + m3*(-s21*s42+s11*s42+s22*s41-s12*s41-s11*s22+s12*s21)) / det

          w4 <- (-s11*(s22*s33-s23*s32) + s21*(s12*s33-s13*s32) - (s12*s23-s13*s22)*s31
                 + m1*(s22*s33-s12*s33-s23*s32+s13*s32+s12*s23-s13*s22)
                 + m2*(-s21*s33+s11*s33+s23*s31-s13*s31-s11*s23+s13*s21)
                 + m3*(s21*s32-s11*s32-s22*s31+s12*s31+s11*s22-s12*s21)) / det

          w <- c(w1, w2, w3, w4)
        } else {
          w <- rep(0.0, 4L)
        }
      }

      # ── Store result ────────────────────────────────────────────────────
      out_id[row_idx]  <- i
      out_gof[row_idx] <- 1.0   # GOF fixed at 1.0 as in C++ original
      out_w[row_idx, ] <- w
      row_idx <- row_idx + 1L
    }
  }

  # ── Build output data frame ────────────────────────────────────────────────
  result <- data.frame(
    ID  = out_id,
    GOF = out_gof,
    out_w,
    stringsAsFactors = FALSE
  )
  colnames(result)[3L:(2L + nsource)] <- make.names(paste0("w ", source_ids))

  if (local_pb) close(pb)
  return(result)
}
