#' @title CTS engine for 4 sources
#'
#' @description Solves determined systems for 3 tracers and 4 sources to assess 
#' solution dispersion and discriminant capacity.
#'
#' @param source Data frame of sediment sources.
#' @param mixture Data frame of mixture samples.
#' @param iter Number of Monte Carlo iterations.
#' @param seed Integer for reproducibility.
#'
#' @return A data frame with one row per tracer triplet, sorted by
#'   \code{max_dw} (most discriminant first), with columns: \code{id}
#'   (tracer triplet), \code{w1} to \code{w4} (source contributions from
#'   mean values), \code{dw1} to \code{dw4} (dispersion of each contribution,
#'   half the 16-84\% quantile range), \code{cons} (proportion of physically
#'   feasible solutions) and \code{max_dw} (maximum dispersion).
#'
#' @export
#' 
CTS_4_sources <- function(source, mixture, iter = 1000, seed = 123456) {
  set.seed(seed)
  source_mat <- data.matrix(source[-1])
  mixture_mat <- data.matrix(mixture[-1])
  
  ntracers <- (ncol(source_mat) - 1) / 2
  tracer_names <- colnames(source_mat)[1:ntracers]
  n_col_idx <- ncol(source_mat)
  
  comb_matrix <- utils::combn(ntracers, 3)
  n_comb <- ncol(comb_matrix)
  
  df_results <- data.frame()
  pb <- txtProgressBar(min = 0, max = n_comb, style = 3, width = 50, char = "=")
  
  for (i in 1:n_comb) {
    t1 <- comb_matrix[1, i]; t2 <- comb_matrix[2, i]; t3 <- comb_matrix[3, i]
    w1 <- numeric(iter); w2 <- numeric(iter); w3 <- numeric(iter); w4 <- numeric(iter)
    sc <- 0
    
    for (l in 1:iter) {
      deg_f <- source_mat[, n_col_idx]
      if (l == 1) {
        s11<-source_mat[1,t1]; s12<-source_mat[1,t2]; s13<-source_mat[1,t3]
        s21<-source_mat[2,t1]; s22<-source_mat[2,t2]; s23<-source_mat[2,t3]
        s31<-source_mat[3,t1]; s32<-source_mat[3,t2]; s33<-source_mat[3,t3]
        s41<-source_mat[4,t1]; s42<-source_mat[4,t2]; s43<-source_mat[4,t3]
      } else {
        # Monte Carlo sampling using t-distribution
        s11 <- source_mat[1,t1] + source_mat[1,t1+ntracers]*rt(1,deg_f[1])/sqrt(deg_f[1])
        s12 <- source_mat[1,t2] + source_mat[1,t2+ntracers]*rt(1,deg_f[1])/sqrt(deg_f[1])
        s13 <- source_mat[1,t3] + source_mat[1,t3+ntracers]*rt(1,deg_f[1])/sqrt(deg_f[1])
        s21 <- source_mat[2,t1] + source_mat[2,t1+ntracers]*rt(1,deg_f[2])/sqrt(deg_f[2])
        s22 <- source_mat[2,t2] + source_mat[2,t2+ntracers]*rt(1,deg_f[2])/sqrt(deg_f[2])
        s23 <- source_mat[2,t3] + source_mat[2,t3+ntracers]*rt(1,deg_f[2])/sqrt(deg_f[2])
        s31 <- source_mat[3,t1] + source_mat[3,t1+ntracers]*rt(1,deg_f[3])/sqrt(deg_f[3])
        s32 <- source_mat[3,t2] + source_mat[3,t2+ntracers]*rt(1,deg_f[3])/sqrt(deg_f[3])
        s33 <- source_mat[3,t3] + source_mat[3,t3+ntracers]*rt(1,deg_f[3])/sqrt(deg_f[3])
        s41 <- source_mat[4,t1] + source_mat[4,t1+ntracers]*rt(1,deg_f[4])/sqrt(deg_f[4])
        s42 <- source_mat[4,t2] + source_mat[4,t2+ntracers]*rt(1,deg_f[4])/sqrt(deg_f[4])
        s43 <- source_mat[4,t3] + source_mat[4,t3+ntracers]*rt(1,deg_f[4])/sqrt(deg_f[4])
      }
      m1 <- mixture_mat[t1]; m2 <- mixture_mat[t2]; m3 <- mixture_mat[t3]
      
      det <- s21*(s32*s43-s33*s42)-s11*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+s11*(s22*s43-s23*s42)+s31*(s12*s43-s13*s42)-s21*(s12*s43-s13*s42)+(s22*s33-s23*s32)*s41-(s12*s33-s13*s32)*s41+(s12*s23-s13*s22)*s41-s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31
      
      if (det != 0) {
        x1 = (s21*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+(s22*s33-s23*s32)*s41 + m1*(-s32*s43+s22*s43+s33*s42-s23*s42-s22*s33+s23*s32) + m2*(s31*s43-s21*s43-s33*s41+s23*s41+s21*s33-s23*s31) + m3*(-s31*s42+s21*s42+s32*s41-s22*s41-s21*s32+s22*s31)) / det
        x2 = (-s11*(s32*s43-s33*s42)+s31*(s12*s43-s13*s42)-(s12*s33-s13*s32)*s41 + m1*(s32*s43-s12*s43-s33*s42+s13*s42+s12*s33-s13*s32) + m2*(-s31*s43+s11*s43+s33*s41-s13*s41-s11*s33+s13*s31) + m3*(s31*s42-s11*s42-s32*s41+s12*s41+s11*s32-s12*s31)) / det
        x3 = (s11*(s22*s43-s23*s42)-s21*(s12*s43-s13*s42)+(s12*s23-s13*s22)*s41 + m1*(-s22*s43+s12*s43+s23*s42-s13*s42-s12*s23+s13*s22) + m2*(s21*s43-s11*s43-s23*s41+s13*s41+s11*s23-s13*s21) + m3*(-s21*s42+s11*s42+s22*s41-s12*s41-s11*s22+s12*s21)) / det
        x4 = (-s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31 + m1*(s22*s33-s12*s33-s23*s32+s13*s32+s12*s23-s13*s22) + m2*(-s21*s33+s11*s33+s23*s31-s13*s31-s11*s23+s13*s21) + m3*(s21*s32-s11*s32-s22*s31+s12*s31+s11*s22-s12*s21)) / det
        
        if (min(x1, x2, x3, x4) >= 0 && max(x1, x2, x3, x4) <= 1) sc <- sc + 1
        w1[l] <- x1; w2[l] <- x2; w3[l] <- x3; w4[l] <- x4
      }
    }
    
    dw <- c(sd(w1), sd(w2), sd(w3), sd(w4)) # O usar el rango intercuartílico como en tus scripts
    res <- data.frame(id = paste(tracer_names[t1], tracer_names[t2], tracer_names[t3], sep=" "),
                      w1=w1[1], w2=w2[1], w3=w3[1], w4=w4[1],
                      dw1=(as.numeric(quantile(w1, 0.84)) - as.numeric(quantile(w1, 0.16)))/2,
                      dw2=(as.numeric(quantile(w2, 0.84)) - as.numeric(quantile(w2, 0.16)))/2,
                      dw3=(as.numeric(quantile(w3, 0.84)) - as.numeric(quantile(w3, 0.16)))/2,
                      dw4=(as.numeric(quantile(w4, 0.84)) - as.numeric(quantile(w4, 0.16)))/2,
                      cons = sc / iter)
    df_results <- rbind(df_results, res)
    setTxtProgressBar(pb, i)
  }
  df_results$max_dw <- apply(df_results[,6:9], 1, max)
  df_results <- df_results[order(df_results$max_dw), ]
  close(pb); return(df_results)
}

# --- Aliases ---
#' @rdname CTS_4_sources
#' @export
triplets <- CTS_4_sources

#' @rdname CTS_4_sources
#' @export
CTS_seeds_triplets <- CTS_4_sources