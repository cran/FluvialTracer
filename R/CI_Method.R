#' Conservativeness Index (CI) for individual tracer evaluation
#'
#' Quantifies how consistently each tracer constrains the mixture apportionment
#' when paired with randomly chosen complementary tracers. For each tracer of
#' interest, the function solves a large number of minimal determined systems
#' and records the spread of the resulting source contributions. A narrow
#' spread indicates that the tracer is conservative and highly informative,
#' regardless of which other tracers accompany it.
#'
#' Two strategies are available for selecting the accompanying tracers:
#' \itemize{
#'   \item \strong{Virtual method}: remaining tracers are replaced by values
#'     drawn uniformly from \eqn{[0, 1]}, isolating the tracer of interest
#'     from correlations with other measured properties.
#'   \item \strong{Dataset method}: remaining tracers are sampled at random
#'     from the measured dataset, reflecting realistic geochemical combinations.
#' }
#'
#' @param data Data frame in the standard FluvialTracer format (sources +
#'   mixture). Raw and averaged formats are both accepted.
#' @param points Integer. Number of Monte Carlo draws per tracer.
#' @param Means Logical or \code{NULL} (default). When \code{NULL}, the
#'   format is detected automatically from column names. Set explicitly to
#'   \code{TRUE} or \code{FALSE} to override.
#' @param seed Integer. Random seed for reproducibility.
#' @param triangles Logical. If \code{TRUE} (default), both virtual and dataset
#'   strategies are computed.
#'
#' @return A list with one data frame per tracer containing the predicted
#'   source contributions across all iterations, followed by two named
#'   elements: \code{DB}, the input data, and \code{CI}, a data frame with
#'   the CI score, normalized score, category and status for each tracer.
#'
#' @references
#' Lizaga, I., Latorre, B., Bodé, S., Gaspar, L., Boeckx, P., & Navas, A. 
#' (2024). Combining isotopic and elemental tracers for enhanced sediment source
#' partitioning in complex catchments. \emph{Journal of Hydrology*}, \emph{631}
#' , 130768.\doi{10.1016/j.jhydrol.2024.130768}
#' 
#' Lizaga, I., Latorre, B., Gaspar, L., & Navas, A. (2020). Consensus ranking
#' as a method to identify non-conservative and dissenting tracers in
#' fingerprinting studies. \emph{Science of The Total Environment}, \emph{720},
#' 137537. \doi{10.1016/j.scitotenv.2020.137537}
#' 
#' @export
#' 
CI_Method <- function(data, points = 2000, Means = NULL, triangles = "virtual", seed = 123456L)
 {
  set.seed(seed)

  # ── Internal helper: detect averaged format (Mean_/SD_ or mean_/sd_) ──────
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
  # Auto-detect if not supplied; legacy TRUE/FALSE still accepted
  if (is.null(Means)) Means <- is_averaged_internal(data)
  
  source_n <- nrow(inputSource(data))
  system.time({
    
    if (source_n == 2) ###########################################################################################
    {
      if (Means == TRUE) {
        sources1 <- data[-nrow(data), c(3:ncol(data))]
        id_sources <- (1:(nrow(data) - 1))
        id <- paste('S', id_sources, sep = '')
        sources <- cbind(id, sources1)
        mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
        
        n_tracers <- ((ncol(data) + 1) / 2) - 2
        pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
        results <- list()
        for (i0 in 0:(n_tracers - 1)) {
            if (triangles == "virtual") {
              results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
                pb = pb, pb_offset = i0 * nrow(mixtures) * points,
                tracer_current = i0 + 1L, tracer_total = n_tracers)
            } else if (triangles == "dataset") {
              results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
                pb = pb, pb_offset = i0 * nrow(mixtures) * points,
                tracer_current = i0 + 1L, tracer_total = n_tracers)
            } else {
              warning("Not applicable")
            }
        }
        close(pb)
      } else {
        ############################################
        sources  <- inputSource(data)
        mixtures <- inputSample(data)
        ###########################################
        n_tracers <- ncol(data) - 2
        pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
        results <- list()
        for (i0 in 0:(n_tracers - 1)) {
          if (triangles == "virtual") {
            results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
              pb = pb, pb_offset = i0 * nrow(mixtures) * points,
              tracer_current = i0 + 1L, tracer_total = n_tracers)
          } else if (triangles == "dataset") {
            results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
              pb = pb, pb_offset = i0 * nrow(mixtures) * points,
              tracer_current = i0 + 1L, tracer_total = n_tracers)
          } else {
            warning("Not applicable")
          }
        }
        close(pb)
      }
      for (i1 in 1:length(results)) {
        results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
        results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
      }
      
      
      # Prepare data for the CI calculation
      ################################################################
      # Removing those out of 0 - 1
      mn_CI <- list() 
      for (i2 in 1:length(results)) {
        new_CI <- results[[i2]][ 
          results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
            results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1  , ]

        mn_CI[[i2]] <-  colMeans(new_CI[, c(3:4)])
      }
      
      CI_new <- list()
      CI_new0 <- c()
      for (i3 in 1:length(results)) {
        CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
        CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
        CI_new0[3]   <- 1 #tmp
        CI_new0[4]   <- 1 #tmp
        
        CI_new[[i3]] <- CI_new0[, c(3, 4, 1, 2)] # forma cutre
        
        # Loop to prevent NA when no solutions are inside the range
        if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 4])) {
          CI_new[[i3]][, c(3, 4)] <- c(-1, -1)
        } else {
          CI_new[[i3]] <- CI_new[[i3]]
        }
      }
      
      # New CI index calculation (distance equation)
      ################################################################
      prct_in <- list()
      for (i4 in 1:length(results)) {
        tmp <- CI_new[[i4]][ 
            CI_new[[i4]][,3] > 0 & CI_new[[i4]][,3] < 1 &
            CI_new[[i4]][,4] > 0 & CI_new[[i4]][,4] < 1  , ]
        prct_in[[i4]] <-  nrow(tmp)/points*100
      }
      tmp1 <- do.call("rbind", prct_in)
      
      for (i3 in 1:length(results)) {
        CI_new0      <- results[[i3]][3]
        CI_new0[2]   <- results[[i3]][4]
        CI_new0[3]   <- 1 #tmp
        CI_new0[4]   <- 1 #tmp
        
        CI_new[[i3]] <- CI_new0[, c(3, 4, 1, 2)] # forma cutre
        
        # Loop to prevent NA when no solutions are inside the range
        if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 4])) {
          CI_new[[i3]][, c(3, 4)] <- c(-1, -1)
        } else {
          CI_new[[i3]] <- CI_new[[i3]]
        }
      }
      
    } else if (source_n == 3) ###########################################################################################
    {
      if (Means == TRUE) {
        sources1 <- data[-nrow(data), c(3:ncol(data))]
        id_sources <- (1:(nrow(data) - 1))
        id <- paste('S', id_sources, sep = '')
        sources <- cbind(id, sources1)
        mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
        
        n_tracers <- ((ncol(data) + 1) / 2) - 2
        pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
        results <- list()
        for (i0 in 0:(n_tracers - 1)) {
          if (triangles == "virtual") {
            results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
              pb = pb, pb_offset = i0 * nrow(mixtures) * points,
              tracer_current = i0 + 1L, tracer_total = n_tracers)
          } else if (triangles == "dataset") {
            results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
              pb = pb, pb_offset = i0 * nrow(mixtures) * points,
              tracer_current = i0 + 1L, tracer_total = n_tracers)
          } else {
            warning("Not applicable")
          }
        }
        close(pb)
      } else {
        ############################################
        sources  <- inputSource(data)
        mixtures <- inputSample(data)
        ###########################################
        n_tracers <- ncol(data) - 2
        pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
        results <- list()
        for (i0 in 0:(n_tracers - 1)) {
          if (triangles == "virtual") {
            results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
              pb = pb, pb_offset = i0 * nrow(mixtures) * points,
              tracer_current = i0 + 1L, tracer_total = n_tracers)
          } else if (triangles == "dataset") {
            results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
              pb = pb, pb_offset = i0 * nrow(mixtures) * points,
              tracer_current = i0 + 1L, tracer_total = n_tracers)
          } else {
            warning("Not applicable")
          }
        }
        close(pb)
      }
      for (i1 in 1:length(results)) {
        results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
        results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
        results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
      }
    
    
    # Prepare data for the CI calculation
    ################################################################
    # Removing those out of 0 - 1
    mn_CI <- list() 
    for (i2 in 1:length(results)) {
      new_CI <- results[[i2]][ 
          results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
          results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1 &
          results[[i2]][,5] > -0.1 & results[[i2]][,5] < 1.1  , ]

      mn_CI[[i2]] <-  colMeans(new_CI[, c(3:5)])
    }
    
    CI_new <- list()
    CI_new0 <- c()
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
      CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
      CI_new0[3]   <- results[[i3]][5] - mn_CI[[i3]][3] + 1 / source_n
      CI_new0[4]   <- 1 #tmp
      CI_new0[5]   <- 1 #tmp
      
      CI_new[[i3]] <- CI_new0[, c(4, 5, 1, 2, 3)] # forma cutre
      
      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 5])) {
        CI_new[[i3]][, c(3, 4, 5)] <- c(-1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
   
    # New CI index calculation (distance equation)
    ################################################################
    prct_in <- list()
    for (i4 in 1:length(results)) {
      tmp <- CI_new[[i4]][ 
          CI_new[[i4]][,3] > 0 & CI_new[[i4]][,3] < 1 &
          CI_new[[i4]][,4] > 0 & CI_new[[i4]][,4] < 1 &
          CI_new[[i4]][,5] > 0 & CI_new[[i4]][,5] < 1  , ]
      prct_in[[i4]] <-  nrow(tmp)/points*100
    }
    tmp1 <- do.call("rbind", prct_in)
    
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3]
      CI_new0[2]   <- results[[i3]][4]
      CI_new0[3]   <- results[[i3]][5]
      CI_new0[4]   <- 1 #tmp
      CI_new0[5]   <- 1 #tmp
      
      CI_new[[i3]] <- CI_new0[, c(4, 5, 1, 2, 3)] # forma cutre
      
      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 5])) {
        CI_new[[i3]][, c(3, 4, 5)] <- c(-1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
  } else if (source_n == 4)	###########################################################################################
    {
    if (Means == TRUE) {
      sources1 <- data[-nrow(data), c(3:ncol(data))]
      id_sources <- (1:(nrow(data) - 1))
      id <- paste('S', id_sources, sep = '')
      sources <- cbind(id, sources1)
      mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
      
      n_tracers <- ((ncol(data) + 1) / 2) - 2
      pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
      results <- list()
      for (i0 in 0:(n_tracers - 1)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else if (triangles == "dataset") {
          results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else {
          warning("Not applicable")
        }
      }
      close(pb)
    } else {
      ############################################
      sources  <- inputSource(data)
      mixtures <- inputSample(data)
      ###########################################
      n_tracers <- ncol(data) - 2
      pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
      results <- list()
      for (i0 in 0:(n_tracers - 1)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else if (triangles == "dataset") {
          results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else {
          warning("Not applicable")
        }
      }
      close(pb)
      }
      for (i1 in 1:length(results)) {
        results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
        results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
        results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
        results[[i1]][, 6] <- as.numeric(as.character(results[[i1]][, 6]))
      }

  # Prepare data for the CI calculation
  ################################################################
  # Removing those out of 0 - 1
      mn_CI <- list() 
      for (i2 in 1:length(results)) {
        new_CI <- results[[i2]][ 
            results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
            results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1 &
            results[[i2]][,5] > -0.1 & results[[i2]][,5] < 1.1 &
            results[[i2]][,6] > -0.1 & results[[i2]][,6] < 1.1  , ]

          mn_CI[[i2]] <-  colMeans(new_CI[, c(3:6)])
      }
    
  CI_new <- list()
  CI_new0 <- c()
    for (i3 in 1:length(results)) {
        CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
        CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
        CI_new0[3]   <- results[[i3]][5] - mn_CI[[i3]][3] + 1 / source_n
        CI_new0[4]   <- results[[i3]][6] - mn_CI[[i3]][4] + 1 / source_n
        CI_new0[5]   <- 1 #tmp
        CI_new0[6]   <- 1 #tmp
        
      CI_new[[i3]] <- CI_new0[, c(5, 6, 1, 2, 3, 4)] # forma cutre
      
      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6])) {
          CI_new[[i3]][, c(3, 4, 5, 6)] <- c(-1, -1, -1, -1)
      } else {
          CI_new[[i3]] <- CI_new[[i3]]
      }
    }

  # New CI index calculation (distance equation)
  ################################################################
    prct_in <- list()
    for (i4 in 1:length(results)) {
      tmp <- CI_new[[i4]][ 
             CI_new[[i4]][,3] > 0 & CI_new[[i4]][,3] < 1 &
             CI_new[[i4]][,4] > 0 & CI_new[[i4]][,4] < 1 &
             CI_new[[i4]][,5] > 0 & CI_new[[i4]][,5] < 1 &
             CI_new[[i4]][,6] > 0 & CI_new[[i4]][,6] < 1  , ]
      prct_in[[i4]] <-  nrow(tmp)/points*100
    }
    tmp1 <- do.call("rbind", prct_in)

    for (i3 in 1:length(results)) {
        CI_new0      <- results[[i3]][3]
        CI_new0[2]   <- results[[i3]][4]
        CI_new0[3]   <- results[[i3]][5]
        CI_new0[4]   <- results[[i3]][6]
        CI_new0[5]   <- 1 #tmp
        CI_new0[6]   <- 1 #tmp
        
      CI_new[[i3]] <- CI_new0[, c(5, 6, 1, 2, 3, 4)] # forma cutre
      
      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6])) {
          CI_new[[i3]][, c(3, 4, 5, 6)] <- c(-1, -1, -1, -1)
      } else {
          CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
    } else if (source_n == 5)	####################################################################################################
      {
    if (Means == TRUE) {
      sources1 <- data[-nrow(data), c(3:ncol(data))]
      id_sources <- (1:(nrow(data) - 1))
      id <- paste('S', id_sources, sep = '')
      sources <- cbind(id, sources1)
      mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
      
      n_tracers <- ((ncol(data) + 1) / 2) - 2
      pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
      results <- list()
      for (i0 in 0:(n_tracers - 1)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else if (triangles == "dataset") {
          results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else {
          warning("Not applicable")
        }
      }
      close(pb)
    } else {
      ############################################
      sources <- inputSource(data)
      mixtures <- inputSample(data)
      ###########################################
      n_tracers <- ncol(data) - 2
      pb <- txtProgressBar(min = 0, max = n_tracers * nrow(mixtures) * points, style = 3, width = 50, char = "=")
      results <- list()
      for (i0 in 0:(n_tracers - 1)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else if (triangles == "dataset") {
          results[[i0 + 1]] <- triangles_dataset(sources, mixtures, i0, points, seed,
            pb = pb, pb_offset = i0 * nrow(mixtures) * points,
            tracer_current = i0 + 1L, tracer_total = n_tracers)
        } else {
          warning("Not applicable")
        }
      }
      close(pb)
    }
    for (i1 in 1:length(results)) {
      results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
      results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
      results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
      results[[i1]][, 6] <- as.numeric(as.character(results[[i1]][, 6]))
      results[[i1]][, 7] <- as.numeric(as.character(results[[i1]][, 7]))
    }

  
  # Prepare data for the CI calculation
  ################################################################
  # Removing those out of 0 - 1
  mn_CI <- list() 
  for (i2 in 1:length(results)) {
    new_CI <- results[[i2]][ 
        results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
        results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1 &
        results[[i2]][,5] > -0.1 & results[[i2]][,5] < 1.1 &
        results[[i2]][,6] > -0.1 & results[[i2]][,6] < 1.1 &
        results[[i2]][,7] > -0.1 & results[[i2]][,7] < 1   , ]
    mn_CI[[i2]] <-  colMeans(new_CI[, c(3:7)])
  }
  
  CI_new <- list()
  CI_new0 <- c()
  for (i3 in 1:length(results)) {
    CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
    CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
    CI_new0[3]   <- results[[i3]][5] - mn_CI[[i3]][3] + 1 / source_n
    CI_new0[4]   <- results[[i3]][6] - mn_CI[[i3]][4] + 1 / source_n
    CI_new0[5]   <- results[[i3]][7] - mn_CI[[i3]][5] + 1 / source_n
    CI_new0[6]   <- 1 #tmp
    CI_new0[7]   <- 1 #tmp
    CI_new[[i3]] <- CI_new0[, c(6, 7, 1, 2, 3, 4, 5)] # forma cutre
    
    # Loop to prevent NA when no solutions are inside the range
    if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6]|is.na(CI_new[[i3]][1, 7]))) {
      CI_new[[i3]][, c(3, 4, 5, 6, 7)] <- c(-1, -1, -1, -1)
    } else {
      CI_new[[i3]] <- CI_new[[i3]]
    }
  }

  # New CI index calculation (distance equation)
  ################################################################
    prct_in <- list()
    for (i4 in 1:length(results)) {
      tmp <- CI_new[[i4]][ 
             CI_new[[i4]][, 3] > 0 & CI_new[[i4]][, 3] < 1 &
             CI_new[[i4]][, 4] > 0 & CI_new[[i4]][, 4] < 1 &
             CI_new[[i4]][, 5] > 0 & CI_new[[i4]][, 5] < 1 &
             CI_new[[i4]][, 6] > 0 & CI_new[[i4]][, 6] < 1 &
             CI_new[[i4]][, 7] > 0 & CI_new[[i4]][, 7] < 1  , ]
  
      prct_in[[i4]] <-  nrow(tmp) / points * 100
    }
    tmp1 <- do.call("rbind", prct_in)
    
  for (i3 in 1:length(results)) {
    CI_new0      <- results[[i3]][3]
    CI_new0[2]   <- results[[i3]][4]
    CI_new0[3]   <- results[[i3]][5]
    CI_new0[4]   <- results[[i3]][6]
    CI_new0[5]   <- results[[i3]][7]
    CI_new0[6]   <- 1 #tmp
    CI_new0[7]   <- 1 #tmp
    CI_new[[i3]] <- CI_new0[, c(6, 7, 1, 2, 3, 4, 5)] # forma cutre
    
    # Loop to prevent NA when no solutions are inside the range
    if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6]|is.na(CI_new[[i3]][1, 7]))) {
      CI_new[[i3]][, c(3, 4, 5, 6, 7)] <- c(-1, -1, -1, -1)
    } else {
      CI_new[[i3]] <- CI_new[[i3]]
    }
  }
    
   } else 
     {
       warning("Not implemented for this number of sources, please contact code developers || mvillalonga@ipe.csic.es")
      }
  
  # Prepare for output
  ################################################################
    # Prepare for output
    ################################################################
    if (Means == TRUE) {
      names_CI_R <- names(data[3:ceiling((ncol(data)/2))])
      names_CI_R <- gsub("Mean_","", names_CI_R)
    }
    if (Means == FALSE) {
      names_CI_R <- names(data[3:ncol(data)])
    }
    
    names_CI_R <- as.data.frame(names_CI_R)
    CI <- cbind(names_CI_R, tmp1)
    CI_new$DB <- data # add the database at the end
    CI_new$CI <- CI   # add the CI index at the end
    colnames(CI_new$CI) <- c("tracer", "CI_score")
    
    # NEW: Dynamic Theoretical Maximum Scaling
    # Benchmarks derived from perfect orthogonal standard simplexes
    ref_max <- switch(as.character(source_n),
                      "2" = 100,
                      "3" = 45.0,
                      "4" = 13.0,
                      "5" = 2.5,
                      100) # Fallback for unsupported dimensions
    # ref_max <- switch(as.character(source_n),
    #                   "2" = 100,
    #                   "3" = 50,
    #                   "4" = 16.67,
    #                   "5" = 4.17,
    #                   100) # Fallback for unsupported dimensions
    
    # # Check if the empirical data managed to beat the orthogonal geometry
    # # This prevents the squashing effect and preserves ranking granularity!
    # actual_max <- max(ref_max, max(CI_new$CI$CI_score, na.rm = TRUE))
    #     # Calculate the standardized score without a hard pmin cap
    # CI_new$CI$CI_score_norm <- round((CI_new$CI$CI_score / actual_max) * 100, 2)
    
    # Calculate the standardized score WITHOUT a ceiling cap!
    # Scores can exceed 100 if the natural dataset geometry is highly favorable.
    # How?
    # Because the benchmark is anchored to a theoretical orthogonal simplex, it 
    # is possible for tracers in natural datasets to score above 100. This occurs 
    # when the specific geometric covariance of the measured sources provides a 
    # more favorable, tightly bound mixing space than a purely orthogonal, generic system.
    CI_new$CI$CI_score_norm <- round((CI_new$CI$CI_score / ref_max) * 100, 2)
    
    # Create Qualitative Ratings using base R 'cut' function
    CI_new$CI$Category <- cut(CI_new$CI$CI_score_norm,
                              breaks = c(-Inf, 20, 40, 60, 80, Inf),
                              labels = c("Highly Non-Conservative", 
                                         "Poor", 
                                         "Moderate", 
                                         "Conservative", 
                                         "Highly Conservative"),
                              right = FALSE)
    
    # Add Visual Traffic Lights
    # CI_new$CI$Status <- cut(CI_new$CI$CI_score_norm,
    #                         breaks = c(-Inf, 20, 40, 60, 80, Inf),
    #                         labels = c("🔴", "🟠", "🟡", "🟢", "🌟"),
    #                         right = FALSE)
    # Add Visual Traffic Lights using Unicode escapes to appease CRAN
    CI_new$CI$Status <- cut(CI_new$CI$CI_score_norm,
                            breaks = c(-Inf, 20, 40, 60, 80, Inf),
                            labels = c("\U0001F534", "\U0001F7E0", "\U0001F7E1", "\U0001F7E2", "\U0001F31F"),
                            right = FALSE)
    # Sort by the normalized score (highest to lowest)
    # CI_new$CI <- CI_new$CI[order(-CI_new$CI$CI_score_norm), ]
    # rownames(CI_new$CI) <- NULL
    
    return(CI_new)      
  })
}


