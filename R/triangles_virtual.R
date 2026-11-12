#' Conservative Index evaluation using virtually generated companion tracers
#'
#' For each tracer of interest, solves a series of minimal determined linear
#' systems to assess how consistently that tracer constrains the source
#' apportionment. In each iteration the tracer under evaluation is kept fixed
#' while the remaining tracers needed to complete the system are replaced by
#' uniformly random values in \eqn{[0, 1]} (virtual approach). The mixture
#' value for each virtual tracer is drawn inside the range spanned by the
#' virtual source values, ensuring the system is geometrically solvable.
#' Because the companion tracers carry no real geochemical information, any
#' consistency in the resulting source weights is attributable solely to the
#' tracer under evaluation, making this a conservative (worst-case) estimator
#' of tracer informativeness.
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
#' @seealso \code{\link{CI_Method}}, \code{\link{triangles_dataset}}
#'
#'
#' @export
#' 
triangles_virtual <- function(sources, mixtures, tracer = 0L, iter = 100L, seed = 69512L,
                                pb = NULL, pb_offset = 0L,
                                tracer_current = NULL, tracer_total = NULL) {

  set.seed(seed)

  # ── Input parsing ──────────────────────────────────────────────────────────
  col_names  <- colnames(sources)
  col_names2 <- colnames(mixtures)
  vars       <- ncol(mixtures) - 1L      # number of tracers
  nsource    <- nrow(sources)
  npoints    <- nrow(mixtures)

  # Basic checks (mirror C++ stops)
  if (ncol(sources) != 2L + (ncol(mixtures) - 1L) * 2L)
    stop("Wrong number of columns in source data frame")
  if (!tolower(col_names[1]) %in% c("id"))
    stop("First data frame column must be ID")
  if (tracer < 0L || tracer >= vars)
    stop("Requested tracer is out of range")
  if (!nsource %in% c(2L, 3L, 4L, 5L))
    stop("Number of sources not implemented")

  # ── Read source means, SDs, and n ─────────────────────────────────────────
  # sources layout: ID | mean_t1 ... mean_tV | sd_t1 ... sd_tV | n
  src_mean <- as.matrix(sources[, 2L:(vars + 1L),          drop = FALSE])  # nsource x vars
  src_sd   <- as.matrix(sources[, (vars + 2L):(vars * 2L + 1L), drop = FALSE])
  src_n    <- as.integer(sources[, ncol(sources)])

  # ── Read mixture values ────────────────────────────────────────────────────
  mix_vals <- as.matrix(mixtures[, 2L:(vars + 1L), drop = FALSE])  # npoints x vars

  # ── Source names for output columns--───────────────────────────────────────
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
  
  # ── Progress bar ───────────────────────────────────────────────────────────
  # If no external pb is passed, create a local one for standalone use
  local_pb <- is.null(pb)
  if (local_pb)
    pb <- txtProgressBar(min = 0, max = npoints * iter, style = 3, width = 50, char = "=")

  # ── Tracer counter message ─────────────────────────────────────────────────
  if (!is.null(tracer_current) && !is.null(tracer_total))
    message(sprintf("Tracer %d out of %d", tracer_current, tracer_total))
  
  # ── Main loop ──────────────────────────────────────────────────────────────
  row_idx <- 1L

  for (i in seq_len(npoints)) {

    mix_i <- mix_vals[i, ]   # mixture values for point i

    for (ii in seq_len(iter)) {

      # Perturb sources with SEM-based t-distribution
      sc <- matrix(0.0, nrow = nsource, ncol = vars)
      for (s in seq_len(nsource)) {
        for (v in seq_len(vars)) {
          sc[s, v] <- correct_sem(src_mean[s, v], src_sd[s, v], src_n[s])
        }
      }
      
      setTxtProgressBar(pb, pb_offset + (i - 1L) * iter + ii)
      
      # Tracer of interest (1-based index internally)
      t1 <- tracer + 1L

      # ── 2 sources ────────────────────────────────────────────────────────
      if (nsource == 2L) {
        s11 <- sc[1L, t1];  s21 <- sc[2L, t1];  m1 <- mix_i[t1]
        det <- s11 - s21
        if (det != 0.0) {
          w1 <- (m1 - s21) / det
          w  <- c(w1, 1.0 - w1)
        } else {
          w <- c(0.0, 0.0)
        }
      }

      # ── 3 sources ────────────────────────────────────────────────────────
      else if (nsource == 3L) {
        s11 <- sc[1L, t1];  s21 <- sc[2L, t1];  s31 <- sc[3L, t1]
        m1  <- mix_i[t1]

        # Build virtual second tracer: 3 uniform values, mixture inside their range
        v2  <- sort(runif(3L))
        s12 <- v2[1L];  s22 <- v2[2L];  s32 <- v2[3L]  # sorted so they span a range
        # Re-draw independently (matching C++ logic: each drawn separately, min/max tracked)
        rnd <- runif(3L)
        s12 <- rnd[1L];  s22 <- rnd[2L];  s32 <- rnd[3L]
        mn2 <- min(rnd);  mx2 <- max(rnd)
        m2  <- mn2 + runif(1L) * (mx2 - mn2)

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
      else if (nsource == 4L) {
        s11 <- sc[1L, t1];  s21 <- sc[2L, t1]
        s31 <- sc[3L, t1];  s41 <- sc[4L, t1];  m1 <- mix_i[t1]

        # Virtual tracer 2
        rnd <- runif(4L)
        s12 <- rnd[1L];  s22 <- rnd[2L];  s32 <- rnd[3L];  s42 <- rnd[4L]
        mn2 <- min(rnd);  mx2 <- max(rnd)
        m2  <- mn2 + runif(1L) * (mx2 - mn2)

        # Virtual tracer 3
        rnd <- runif(4L)
        s13 <- rnd[1L];  s23 <- rnd[2L];  s33 <- rnd[3L];  s43 <- rnd[4L]
        mn3 <- min(rnd);  mx3 <- max(rnd)
        m3  <- mn3 + runif(1L) * (mx3 - mn3)

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

      # ── 5 sources ────────────────────────────────────────────────────────
      else {  # nsource == 5
        s11 <- sc[1L, t1];  s21 <- sc[2L, t1];  s31 <- sc[3L, t1]
        s41 <- sc[4L, t1];  s51 <- sc[5L, t1];  m1  <- mix_i[t1]

        # Virtual tracer 2
        rnd <- runif(5L)
        s12 <- rnd[1L]; s22 <- rnd[2L]; s32 <- rnd[3L]; s42 <- rnd[4L]; s52 <- rnd[5L]
        mn2 <- min(rnd); mx2 <- max(rnd); m2 <- mn2 + runif(1L)*(mx2-mn2)

        # Virtual tracer 3
        rnd <- runif(5L)
        s13 <- rnd[1L]; s23 <- rnd[2L]; s33 <- rnd[3L]; s43 <- rnd[4L]; s53 <- rnd[5L]
        mn3 <- min(rnd); mx3 <- max(rnd); m3 <- mn3 + runif(1L)*(mx3-mn3)

        # Virtual tracer 4
        rnd <- runif(5L)
        s14 <- rnd[1L]; s24 <- rnd[2L]; s34 <- rnd[3L]; s44 <- rnd[4L]; s54 <- rnd[5L]
        mn4 <- min(rnd); mx4 <- max(rnd); m4 <- mn4 + runif(1L)*(mx4-mn4)

        det <- s21*(s32*(s43*s54-s44*s53)-s42*(s33*s54-s34*s53)+(s33*s44-s34*s43)*s52)-
               s11*(s32*(s43*s54-s44*s53)-s42*(s33*s54-s34*s53)+(s33*s44-s34*s43)*s52)-
               s31*(s22*(s43*s54-s44*s53)-s42*(s23*s54-s24*s53)+(s23*s44-s24*s43)*s52)+
               s11*(s22*(s43*s54-s44*s53)-s42*(s23*s54-s24*s53)+(s23*s44-s24*s43)*s52)+
               s31*(s12*(s43*s54-s44*s53)-s42*(s13*s54-s14*s53)+(s13*s44-s14*s43)*s52)-
               s21*(s12*(s43*s54-s44*s53)-s42*(s13*s54-s14*s53)+(s13*s44-s14*s43)*s52)+
               s41*(s22*(s33*s54-s34*s53)-s32*(s23*s54-s24*s53)+(s23*s34-s24*s33)*s52)-
               s11*(s22*(s33*s54-s34*s53)-s32*(s23*s54-s24*s53)+(s23*s34-s24*s33)*s52)-
               s41*(s12*(s33*s54-s34*s53)-s32*(s13*s54-s14*s53)+(s13*s34-s14*s33)*s52)+
               s21*(s12*(s33*s54-s34*s53)-s32*(s13*s54-s14*s53)+(s13*s34-s14*s33)*s52)+
               s41*(s12*(s23*s54-s24*s53)-s22*(s13*s54-s14*s53)+(s13*s24-s14*s23)*s52)-
               s31*(s12*(s23*s54-s24*s53)-s22*(s13*s54-s14*s53)+(s13*s24-s14*s23)*s52)-
               (s22*(s33*s44-s34*s43)-s32*(s23*s44-s24*s43)+(s23*s34-s24*s33)*s42)*s51+
               (s12*(s33*s44-s34*s43)-s32*(s13*s44-s14*s43)+(s13*s34-s14*s33)*s42)*s51-
               (s12*(s23*s44-s24*s43)-s22*(s13*s44-s14*s43)+(s13*s24-s14*s23)*s42)*s51+
               (s12*(s23*s34-s24*s33)-s22*(s13*s34-s14*s33)+(s13*s24-s14*s23)*s32)*s51+
               s11*(s22*(s33*s44-s34*s43)-s32*(s23*s44-s24*s43)+(s23*s34-s24*s33)*s42)-
               s21*(s12*(s33*s44-s34*s43)-s32*(s13*s44-s14*s43)+(s13*s34-s14*s33)*s42)+
               s31*(s12*(s23*s44-s24*s43)-s22*(s13*s44-s14*s43)+(s13*s24-s14*s23)*s42)-
               (s12*(s23*s34-s24*s33)-s22*(s13*s34-s14*s33)+(s13*s24-s14*s23)*s32)*s41

        if (det != 0.0) {
          w1 <- (m1*(((s22-s32)*s43+(s33-s23)*s42-s22*s33+s23*s32)*s54+((s32-s22)*s44+(s24-s34)*s42+s22*s34-s24*s32)*s53+((s23-s33)*s44+(s34-s24)*s43-s23*s34+s24*s33)*s52+(s22*s33-s23*s32)*s44+(s24*s32-s22*s34)*s43+(s23*s34-s24*s33)*s42)+m2*(((s31-s21)*s43+(s23-s33)*s41+s21*s33-s23*s31)*s54+((s21-s31)*s44+(s34-s24)*s41-s21*s34+s24*s31)*s53+((s33-s23)*s44+(s24-s34)*s43+s23*s34-s24*s33)*s51+(s23*s31-s21*s33)*s44+(s21*s34-s24*s31)*s43+(s24*s33-s23*s34)*s41)+m3*(((s21-s31)*s42+(s32-s22)*s41-s21*s32+s22*s31)*s54+((s31-s21)*s44+(s24-s34)*s41+s21*s34-s24*s31)*s52+((s22-s32)*s44+(s34-s24)*s42-s22*s34+s24*s32)*s51+(s21*s32-s22*s31)*s44+(s24*s31-s21*s34)*s42+(s22*s34-s24*s32)*s41)+((s21*s32-s22*s31)*s43+(s23*s31-s21*s33)*s42+(s22*s33-s23*s32)*s41)*s54+m4*(((s31-s21)*s42+(s22-s32)*s41+s21*s32-s22*s31)*s53+((s21-s31)*s43+(s33-s23)*s41-s21*s33+s23*s31)*s52+((s32-s22)*s43+(s23-s33)*s42+s22*s33-s23*s32)*s51+(s22*s31-s21*s32)*s43+(s21*s33-s23*s31)*s42+(s23*s32-s22*s33)*s41)+((s22*s31-s21*s32)*s44+(s21*s34-s24*s31)*s42+(s24*s32-s22*s34)*s41)*s53+((s21*s33-s23*s31)*s44+(s24*s31-s21*s34)*s43+(s23*s34-s24*s33)*s41)*s52+((s23*s32-s22*s33)*s44+(s22*s34-s24*s32)*s43+(s24*s33-s23*s34)*s42)*s51) / det

          w2 <- (m1*(((s32-s12)*s43+(s13-s33)*s42+s12*s33-s13*s32)*s54+((s12-s32)*s44+(s34-s14)*s42-s12*s34+s14*s32)*s53+((s33-s13)*s44+(s14-s34)*s43+s13*s34-s14*s33)*s52+(s13*s32-s12*s33)*s44+(s12*s34-s14*s32)*s43+(s14*s33-s13*s34)*s42)+m2*(((s11-s31)*s43+(s33-s13)*s41-s11*s33+s13*s31)*s54+((s31-s11)*s44+(s14-s34)*s41+s11*s34-s14*s31)*s53+((s13-s33)*s44+(s34-s14)*s43-s13*s34+s14*s33)*s51+(s11*s33-s13*s31)*s44+(s14*s31-s11*s34)*s43+(s13*s34-s14*s33)*s41)+m3*(((s31-s11)*s42+(s12-s32)*s41+s11*s32-s12*s31)*s54+((s11-s31)*s44+(s34-s14)*s41-s11*s34+s14*s31)*s52+((s32-s12)*s44+(s14-s34)*s42+s12*s34-s14*s32)*s51+(s12*s31-s11*s32)*s44+(s11*s34-s14*s31)*s42+(s14*s32-s12*s34)*s41)+((s12*s31-s11*s32)*s43+(s11*s33-s13*s31)*s42+(s13*s32-s12*s33)*s41)*s54+m4*(((s11-s31)*s42+(s32-s12)*s41-s11*s32+s12*s31)*s53+((s31-s11)*s43+(s13-s33)*s41+s11*s33-s13*s31)*s52+((s12-s32)*s43+(s33-s13)*s42-s12*s33+s13*s32)*s51+(s11*s32-s12*s31)*s43+(s13*s31-s11*s33)*s42+(s12*s33-s13*s32)*s41)+((s11*s32-s12*s31)*s44+(s14*s31-s11*s34)*s42+(s12*s34-s14*s32)*s41)*s53+((s13*s31-s11*s33)*s44+(s11*s34-s14*s31)*s43+(s14*s33-s13*s34)*s41)*s52+((s12*s33-s13*s32)*s44+(s14*s32-s12*s34)*s43+(s13*s34-s14*s33)*s42)*s51) / det

          w3 <- (m1*(((s12-s22)*s43+(s23-s13)*s42-s12*s23+s13*s22)*s54+((s22-s12)*s44+(s14-s24)*s42+s12*s24-s14*s22)*s53+((s13-s23)*s44+(s24-s14)*s43-s13*s24+s14*s23)*s52+(s12*s23-s13*s22)*s44+(s14*s22-s12*s24)*s43+(s13*s24-s14*s23)*s42)+m2*(((s21-s11)*s43+(s13-s23)*s41+s11*s23-s13*s21)*s54+((s11-s21)*s44+(s24-s14)*s41-s11*s24+s14*s21)*s53+((s23-s13)*s44+(s14-s24)*s43+s13*s24-s14*s23)*s51+(s13*s21-s11*s23)*s44+(s11*s24-s14*s21)*s43+(s14*s23-s13*s24)*s41)+m3*(((s11-s21)*s42+(s22-s12)*s41-s11*s22+s12*s21)*s54+((s21-s11)*s44+(s14-s24)*s41+s11*s24-s14*s21)*s52+((s12-s22)*s44+(s24-s14)*s42-s12*s24+s14*s22)*s51+(s11*s22-s12*s21)*s44+(s14*s21-s11*s24)*s42+(s12*s24-s14*s22)*s41)+((s11*s22-s12*s21)*s43+(s13*s21-s11*s23)*s42+(s12*s23-s13*s22)*s41)*s54+m4*(((s21-s11)*s42+(s12-s22)*s41+s11*s22-s12*s21)*s53+((s11-s21)*s43+(s23-s13)*s41-s11*s23+s13*s21)*s52+((s22-s12)*s43+(s13-s23)*s42+s12*s23-s13*s22)*s51+(s12*s21-s11*s22)*s43+(s11*s23-s13*s21)*s42+(s13*s22-s12*s23)*s41)+((s12*s21-s11*s22)*s44+(s11*s24-s14*s21)*s42+(s14*s22-s12*s24)*s41)*s53+((s11*s23-s13*s21)*s44+(s14*s21-s11*s24)*s43+(s13*s24-s14*s23)*s41)*s52+((s13*s22-s12*s23)*s44+(s12*s24-s14*s22)*s43+(s14*s23-s13*s24)*s42)*s51) / det

          w4 <- (m1*(((s22-s12)*s33+(s13-s23)*s32+s12*s23-s13*s22)*s54+((s12-s22)*s34+(s24-s14)*s32-s12*s24+s14*s22)*s53+((s23-s13)*s34+(s14-s24)*s33+s13*s24-s14*s23)*s52+(s13*s22-s12*s23)*s34+(s12*s24-s14*s22)*s33+(s14*s23-s13*s24)*s32)+m2*(((s11-s21)*s33+(s23-s13)*s31-s11*s23+s13*s21)*s54+((s21-s11)*s34+(s14-s24)*s31+s11*s24-s14*s21)*s53+((s13-s23)*s34+(s24-s14)*s33-s13*s24+s14*s23)*s51+(s11*s23-s13*s21)*s34+(s14*s21-s11*s24)*s33+(s13*s24-s14*s23)*s31)+m3*(((s21-s11)*s32+(s12-s22)*s31+s11*s22-s12*s21)*s54+((s11-s21)*s34+(s24-s14)*s31-s11*s24+s14*s21)*s52+((s22-s12)*s34+(s14-s24)*s32+s12*s24-s14*s22)*s51+(s12*s21-s11*s22)*s34+(s11*s24-s14*s21)*s32+(s14*s22-s12*s24)*s31)+((s12*s21-s11*s22)*s33+(s11*s23-s13*s21)*s32+(s13*s22-s12*s23)*s31)*s54+m4*(((s11-s21)*s32+(s22-s12)*s31-s11*s22+s12*s21)*s53+((s21-s11)*s33+(s13-s23)*s31+s11*s23-s13*s21)*s52+((s12-s22)*s33+(s23-s13)*s32-s12*s23+s13*s22)*s51+(s11*s22-s12*s21)*s33+(s13*s21-s11*s23)*s32+(s12*s23-s13*s22)*s31)+((s11*s22-s12*s21)*s34+(s14*s21-s11*s24)*s32+(s12*s24-s14*s22)*s31)*s53+((s13*s21-s11*s23)*s34+(s11*s24-s14*s21)*s33+(s14*s23-s13*s24)*s31)*s52+((s12*s23-s13*s22)*s34+(s14*s22-s12*s24)*s33+(s13*s24-s14*s23)*s32)*s51) / det

          w5 <- (m1*(((s12-s22)*s33+(s23-s13)*s32-s12*s23+s13*s22)*s44+((s22-s12)*s34+(s14-s24)*s32+s12*s24-s14*s22)*s43+((s13-s23)*s34+(s24-s14)*s33-s13*s24+s14*s23)*s42+(s12*s23-s13*s22)*s34+(s14*s22-s12*s24)*s33+(s13*s24-s14*s23)*s32)+m2*(((s21-s11)*s33+(s13-s23)*s31+s11*s23-s13*s21)*s44+((s11-s21)*s34+(s24-s14)*s31-s11*s24+s14*s21)*s43+((s23-s13)*s34+(s14-s24)*s33+s13*s24-s14*s23)*s41+(s13*s21-s11*s23)*s34+(s11*s24-s14*s21)*s33+(s14*s23-s13*s24)*s31)+m3*(((s11-s21)*s32+(s22-s12)*s31-s11*s22+s12*s21)*s44+((s21-s11)*s34+(s14-s24)*s31+s11*s24-s14*s21)*s42+((s12-s22)*s34+(s24-s14)*s32-s12*s24+s14*s22)*s41+(s11*s22-s12*s21)*s34+(s14*s21-s11*s24)*s32+(s12*s24-s14*s22)*s31)+((s11*s22-s12*s21)*s33+(s13*s21-s11*s23)*s32+(s12*s23-s13*s22)*s31)*s44+m4*(((s21-s11)*s32+(s12-s22)*s31+s11*s22-s12*s21)*s43+((s11-s21)*s33+(s23-s13)*s31-s11*s23+s13*s21)*s42+((s22-s12)*s33+(s13-s23)*s32+s12*s23-s13*s22)*s41+(s12*s21-s11*s22)*s33+(s11*s23-s13*s21)*s32+(s13*s22-s12*s23)*s31)+((s12*s21-s11*s22)*s34+(s11*s24-s14*s21)*s32+(s14*s22-s12*s24)*s31)*s43+((s11*s23-s13*s21)*s34+(s14*s21-s11*s24)*s33+(s13*s24-s14*s23)*s31)*s42+((s13*s22-s12*s23)*s34+(s12*s24-s14*s22)*s33+(s14*s23-s13*s24)*s32)*s41) / det

          w <- c(w1, w2, w3, w4, w5)
        } else {
          w <- rep(0.0, 5L)
        }
      }

      # ── Store result ────────────────────────────────────────────────────
      out_id[row_idx]   <- i
      out_gof[row_idx]  <- 1.0   # GOF fixed at 1.0 as in C++ original
      out_w[row_idx, ]  <- w
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
