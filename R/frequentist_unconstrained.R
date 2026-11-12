#' Source apportionment via unconstrained least-squares with classical variability propagation
#'
#' Estimates the fractional contribution of each sediment source to a mixture
#' by solving an overdetermined mass-balance system with least-squares
#' minimisation. Source variability is propagated by Monte Carlo resampling:
#' at each iteration, synthetic source signatures are drawn from a normal
#' distribution defined by the measured mean and standard deviation, and the
#' system is solved independently for each draw.
#'
#' This function is the computational back-end called by \code{\link{unmixing}}.
#' For most users, calling \code{unmixing} directly is recommended.
#'
#' @param data Data frame in the standard FluvialTracer format (sources + mixture).
#' @param variability Character. \code{"SEM"} for Standard Error of the Mean,
#'   \code{"SD"} for Standard Deviation.
#' @param iter Integer. Number of Monte Carlo iterations.
#' @param means Logical. \code{TRUE} when the input is in averaged format
#'   (\code{Mean_}/\code{SD_} columns); \code{FALSE} for raw observations.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A data frame with one row per iteration containing the mixture ID,
#'   goodness-of-fit (GOF), and one column per source with the estimated
#'   fractional contribution.
#'
#' @export
#' 
frequentist_unconstrained <- function(data, variability = "SEM", iter = 1000, means = FALSE, seed = 123456L)
{
  set.seed(seed)
  
  sr_names <- unique(data[, 2])
  
  if(is.numeric(iter)){ iter <- iter} else {
    if(iter == "very short") iter <- 500
    if(iter == "short")      iter <- 1000
    if(iter == "medium")     iter <- 5000
    if(iter == "long")       iter <- 20000
    if(iter == "bananas")    iter <- 100000
  }
  
  if (means == TRUE) {
    source <- data[c(1:(nrow(data) - 1)), c(2:ncol(data))]
    x <- round((ncol(data) - 3) / 2 + 2)
    mixture <- inputMixture(data[nrow(data), 1:x])
    colnames(source)[1] <- "id"
    colnames(mixture) <- colnames(source[, c(1:(ncol(source) / 2))])
    
  } else {
    ###########################################
    source <- inputSource(data)
    mixture <- inputMixture(data)
    ###########################################
  }
  
  if (nrow(source)-1 >= ncol(mixture)) { #ntracer >= nsources - 1
    warning("As a minimum, n - 1 properties are required to discriminate rigorously between n sources. Additional properties are frequently required to increase the reliability of the results")
  }
  
  snames <- source[,1]
  
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])
  
  # normalize
  cols <- (ncol(source)-1)/2
  for (col in c(1:cols))
  {
    mx <- max(source[, col] + source[, cols + col])
    mn <- min(source[, col] - source[, cols + col])
    source[, col] <- (source[, col] - mn) / (mx - mn)
    source[, cols + col] <- source[, cols + col] / (mx - mn)
    mixture[, col] <- (mixture[, col] - mn) / (mx - mn)
  }
  
  nsource <- nrow(source)
  ntracer <- (ncol(source)-1)/2
  
  r <- matrix(, nrow = iter, ncol = nsource+2)
  
  csource <- matrix(, nrow = nsource, ncol = ntracer)
  
  
  message(crayon::cyan(crayon::bold("Summary of the model inputs:")))
  message(ntracer, " variables from ", nsource, " sources")
  
  
  # Introduce the progress bar
  pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=")
  
  for (i in c(1:iter))
  {
    y <- c()
    x <- matrix(, nrow = ntracer, ncol = nsource-1)
    
    if(i==2 || i==3) #central solution
    {
      for (j in c(1:ntracer))
      {
        y <- c(y, mixture[1,j][[1]]-source[nsource,j])
        for (k in c(1:(nsource-1)))
        {
          x[j, k] <- source[k,j]-source[nsource,j]
        }
      }
      for (j in c(1:ntracer))
      {
        for (k in c(1:nsource))
        {
          csource[k,j] <- source[k,j]
        }
      }
    }
    else
    {
      ls <- c()
      for (j in c(1:ntracer))
      {
        
        if (variability == "SEM") {
          x1 <- rt(1, source[nsource, ntracer*2+1]) / sqrt(source[nsource, ntracer*2+1])
          
        } else if (variability == "SD") {
          x1 <- rt(1, source[nsource, ntracer*2+1])
          
        } else  {
          warning("'variability' must be \"SEM\" or \"SD\". Please contact code developers.", call. = FALSE)
          
        }
        
        csource[nsource, j] <- source[nsource, j] + source[nsource, ntracer+j] * x1
        ls <- c(ls, csource[nsource, j])
      }
      
      for (j in c(1:ntracer))
      {
        y <- c(y, mixture[1,j][[1]] - ls[j])
        for (k in c(1:(nsource-1)))
        {
          if (variability == "SEM") {
            x1 <- rt(1, source[k, ntracer * 2 + 1]) / sqrt(source[k, ntracer * 2 + 1])
            
          } else if (variability == "SD") {
            x1 <- rt(1, source[k, ntracer * 2 + 1])
            
          } else  {
            warning("'variability' must be \"SEM\" or \"SD\". Please contact code developers.", call. = FALSE)
            
          }
          csource[k,j] <- source[k,j] + source[k,ntracer+j] * x1
          x[j, k] <- csource[k,j] - ls[j]
        }
      }
    }
    
    # least squares method
    model <- lm.fit(x, y)
    w <- as.vector(coef(model))
    w <- c(w, 1-sum(w))
    
    gof <- c()
    for (j in c(1:ntracer))
    {
      x1 <- 0.0
      for (k in c(1:nsource))
      {
        x1 <- x1 + w[k] * csource[k,j]
      }
      gof <- c(gof, (mixture[1,j][[1]]-x1)^2)
    }
    
    gof <- 1.0 - mean(gof)
    gof <- max(gof, 0)
    
    w <- c(1, gof, w)
    
    r[i, ] <- w
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  r <- as.data.frame(r)
  # colnames(r) <- c('id', 'gof', paste( "w.", snames[c(1:nsource)], sep=""))
  colnames(r) <- c('id', 'gof', paste( sr_names[c(1:nsource)], sep=""))
  
  r2 <- r
  # r2 <- r2[, c(setdiff(colnames(r2), "gof"))]
  
  message(crayon::green(crayon::bold("Summary of the model outputs:")))
  message("The results below show the outcome of the unmixing process, which incorporated a source variability analysis \nover ",
          iter, " iterations. Notice that the second and third rows of the result correspond to the solution for the \ncentral or mean value of the sources.")
  
  summary_tab <- aggregate(. ~ id, data = r2, function(x) c(mean = mean(x), SD = sd(x)))
  message(paste(utils::capture.output(print(summary_tab)), collapse = "\n"))
  
  invisible(r)
}