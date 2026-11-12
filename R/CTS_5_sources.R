#' @title CTS engine for 5 sources
#'
#' @description Solves determined systems for 4 tracers and 5 sources to assess 
#' solution dispersion and discriminant capacity.
#'
#' @param source Data frame of sediment sources.
#' @param mixture Data frame of mixture samples.
#' @param iter Number of Monte Carlo iterations.
#' @param seed Integer for reproducibility.
#'
#' @return A data frame with one row per tracer quartet, sorted by
#'   \code{max_dw} (most discriminant first), with columns: \code{id}
#'   (tracer quartet), \code{w1} to \code{w5} (source contributions from
#'   mean values), \code{dw1} to \code{dw5} (dispersion of each contribution,
#'   half the 16-84\% quantile range), \code{cons} (proportion of physically
#'   feasible solutions) and \code{max_dw} (maximum dispersion).
#'
#' @export
#' 
CTS_5_sources <- function(source, mixture, iter = 1000, seed = 123456) {
  set.seed(seed)
  source_mat <- data.matrix(source[-1])
  mixture_mat <- data.matrix(mixture[-1])
  
  ntracers <- (ncol(source_mat) - 1) / 2
  tracer_names <- colnames(source_mat)[1:ntracers]
  n_col_idx <- ncol(source_mat)
  
  comb_matrix <- utils::combn(ntracers, 4)
  n_comb <- ncol(comb_matrix)
  
  df_results <- data.frame()
  pb <- txtProgressBar(min = 0, max = n_comb, style = 3, width = 50, char = "=")
  
  for (i in 1:n_comb) {
    t1 <- comb_matrix[1, i]; t2 <- comb_matrix[2, i]
    t3 <- comb_matrix[3, i]; t4 <- comb_matrix[4, i]
    
    w1_v <- numeric(iter); w2_v <- numeric(iter); w3_v <- numeric(iter); w4_v <- numeric(iter); w5_v <- numeric(iter)
    sc <- 0
    
    deg_f <- source_mat[, n_col_idx]
    
    for (l in 1:iter) {
      if (l == 1) {
        s11<-source_mat[1,t1]; s12<-source_mat[1,t2]; s13<-source_mat[1,t3]; s14<-source_mat[1,t4]
        s21<-source_mat[2,t1]; s22<-source_mat[2,t2]; s23<-source_mat[2,t3]; s24<-source_mat[2,t4]
        s31<-source_mat[3,t1]; s32<-source_mat[3,t2]; s33<-source_mat[3,t3]; s34<-source_mat[3,t4]
        s41<-source_mat[4,t1]; s42<-source_mat[4,t2]; s43<-source_mat[4,t3]; s44<-source_mat[4,t4]
        s51<-source_mat[5,t1]; s52<-source_mat[5,t2]; s53<-source_mat[5,t3]; s54<-source_mat[5,t4]
      } else {
        s11 <- source_mat[1,t1] + source_mat[1,t1+ntracers]*rt(1,deg_f[1])/sqrt(deg_f[1])
        s12 <- source_mat[1,t2] + source_mat[1,t2+ntracers]*rt(1,deg_f[1])/sqrt(deg_f[1])
        s13 <- source_mat[1,t3] + source_mat[1,t3+ntracers]*rt(1,deg_f[1])/sqrt(deg_f[1])
        s14 <- source_mat[1,t4] + source_mat[1,t4+ntracers]*rt(1,deg_f[1])/sqrt(deg_f[1])
        s21 <- source_mat[2,t1] + source_mat[2,t1+ntracers]*rt(1,deg_f[2])/sqrt(deg_f[2])
        s22 <- source_mat[2,t2] + source_mat[2,t2+ntracers]*rt(1,deg_f[2])/sqrt(deg_f[2])
        s23 <- source_mat[2,t3] + source_mat[2,t3+ntracers]*rt(1,deg_f[2])/sqrt(deg_f[2])
        s24 <- source_mat[2,t4] + source_mat[2,t4+ntracers]*rt(1,deg_f[2])/sqrt(deg_f[2])
        s31 <- source_mat[3,t1] + source_mat[3,t1+ntracers]*rt(1,deg_f[3])/sqrt(deg_f[3])
        s32 <- source_mat[3,t2] + source_mat[3,t2+ntracers]*rt(1,deg_f[3])/sqrt(deg_f[3])
        s33 <- source_mat[3,t3] + source_mat[3,t3+ntracers]*rt(1,deg_f[3])/sqrt(deg_f[3])
        s34 <- source_mat[3,t4] + source_mat[3,t4+ntracers]*rt(1,deg_f[3])/sqrt(deg_f[3])
        s41 <- source_mat[4,t1] + source_mat[4,t1+ntracers]*rt(1,deg_f[4])/sqrt(deg_f[4])
        s42 <- source_mat[4,t2] + source_mat[4,t2+ntracers]*rt(1,deg_f[4])/sqrt(deg_f[4])
        s43 <- source_mat[4,t3] + source_mat[4,t3+ntracers]*rt(1,deg_f[4])/sqrt(deg_f[4])
        s44 <- source_mat[4,t4] + source_mat[4,t4+ntracers]*rt(1,deg_f[4])/sqrt(deg_f[4])
        s51 <- source_mat[5,t1] + source_mat[5,t1+ntracers]*rt(1,deg_f[5])/sqrt(deg_f[5])
        s52 <- source_mat[5,t2] + source_mat[5,t2+ntracers]*rt(1,deg_f[5])/sqrt(deg_f[5])
        s53 <- source_mat[5,t3] + source_mat[5,t3+ntracers]*rt(1,deg_f[5])/sqrt(deg_f[5])
        s54 <- source_mat[5,t4] + source_mat[5,t4+ntracers]*rt(1,deg_f[5])/sqrt(deg_f[5])
      }
      m1 <- mixture_mat[t1]; m2 <- mixture_mat[t2]; m3 <- mixture_mat[t3]; m4 <- mixture_mat[t4]
      
      A <- matrix(c(
        1,   1,   1,   1,   1,
        s11, s21, s31, s41, s51,
        s12, s22, s32, s42, s52,
        s13, s23, s33, s43, s53,
        s14, s24, s34, s44, s54
      ), nrow = 5, ncol = 5, byrow = TRUE)
      
      B <- c(1, m1, m2, m3, m4)
      sol <- try(solve(A, B), silent = TRUE)
      
      if (!inherits(sol, "try-error")) {
        if (min(sol) >= 0 && max(sol) <= 1) sc <- sc + 1
        w1_v[l] <- sol[1]; w2_v[l] <- sol[2]; w3_v[l] <- sol[3]; w4_v[l] <- sol[4]; w5_v[l] <- sol[5]
      }
    }
    
    res <- data.frame(
      id = paste(tracer_names[t1], tracer_names[t2], tracer_names[t3], tracer_names[t4], sep=" "),
      w1=w1_v[1], w2=w2_v[1], w3=w3_v[1], w4=w4_v[1], w5=w5_v[1],
      dw1=(as.numeric(quantile(w1_v, 0.84)) - as.numeric(quantile(w1_v, 0.16)))/2,
      dw2=(as.numeric(quantile(w2_v, 0.84)) - as.numeric(quantile(w2_v, 0.16)))/2,
      dw3=(as.numeric(quantile(w3_v, 0.84)) - as.numeric(quantile(w3_v, 0.16)))/2,
      dw4=(as.numeric(quantile(w4_v, 0.84)) - as.numeric(quantile(w4_v, 0.16)))/2,
      dw5=(as.numeric(quantile(w5_v, 0.84)) - as.numeric(quantile(w5_v, 0.16)))/2,
      cons = sc / iter
    )
    df_results <- rbind(df_results, res)
    setTxtProgressBar(pb, i)
  }
  
  df_results$max_dw <- apply(df_results[, 7:11], 1, max)
  df_results <- df_results[order(df_results$max_dw), ]
  rownames(df_results) <- NULL
  
  close(pb); return(df_results)
}
# Alias
#' @rdname CTS_5_sources
#' @export
quartets <- CTS_5_sources

#' @rdname CTS_5_sources
#' @export
CTS_seeds_quartets <- CTS_5_sources