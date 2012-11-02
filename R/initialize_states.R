#' @title State initialization for iterative algorithms (randomly or variants of kmeans)
#'
#' @description 
#' Initializes the state/cluster assignment either uniformly at random from
#' \eqn{K} classes, or using initial \emph{kmeans++} (\code{\link{kmeanspp}}) 
#' clustering (in several variations on PLCs and/or FLCs).
#' 
#' @param nstates number of states 
#' @param nsamples number of samples.
#' @param method how to choose the labels: either uniformly at random from 
#' \eqn{\lbrace 1, \ldots, K \rbrace} or using K-means on PLCs and FLCs or 
#' a combination.  Default: \code{method = "random"}.  Other options are
#' \code{c("KmeansPLC","KmeansFLC","KmeansPLCFLC","KmeansFLCPLC")}
#' @param LCs (optional) a list of \code{PLC} (\eqn{N \times n_p} array) and 
#' \code{FLC} (\eqn{N \times n_f} array)
#' @keywords datagen distribution multivariate
#' @export
#' @examples
#' x1 = rnorm(1000)
#' x2 = rnorm(200, mean = 2)
#' yy = c(x1, x2)
#' ss = initialize_states(nstates = 2, nsamples = length(yy), method = "KmeansFLC", 
#'                        LCs = list(FLCs = yy))
#' plot(yy, col = ss, pch = 19)
#' points(x1, col = "blue")
#' 

initialize_states <- function(nstates = NULL, nsamples = NULL, 
                              method = "random", 
                              LCs = list(PLC = NULL, FLC = NULL)) {
  if (is.null(nstates)) {
    stop("You must provide the number of clusters.")
  }
  
  if (is.null(unlist(LCs)) & is.null(nsamples)){
    stop("You must either provide the total number of samples or data 'LCs'.")
  } 
  
  if (is.null(nsamples)) {
    stop("You must provide the total number of samples.")
  } 
  
  state_vector <- rep(NA, nsamples)
  if (method == "random") {
    state_vector <- sample.int(n = nstates, size = nsamples, replace = TRUE)
    # make sure that each state appears at least once
    state_vector[sample.int(nsamples, nstates, replace = FALSE)] <- 1:nstates
  } else if (method == "KmeansPLC") {
    state_vector <- kmeanspp(LCs$PLC, nstates, iter.max = 100, nstart = 10)$cluster
  } else if (method == "KmeansFLC") {
    state_vector <- kmeanspp(LCs$FLC, nstates, iter.max = 100, nstart = 10)$cluster
  } else if (any(method == c("KmeansPLCFLC", "KmeansFLCPLC"))) {
    # do two stage clustering 1) on the PLC (FLC) space, and then conditional
    # on the first stage cluster 2) do a clustering in each cluster but using
    # the FLCs (PLCs)
    first_stage_nstates <- ceiling(nstates^(1/3))
    second_stage_nstates <- floor(nstates/first_stage_nstates)
    second_stage_last_run_nstates <- nstates - first_stage_nstates * second_stage_nstates + second_stage_nstates
    
    
    if (method == "KmeansFLCPLC") {
      first_stage_data <- rbind(LCs$FLC)
      second_stage_data <- rbind(LCs$PLC)
    } else {
      first_stage_data <- rbind(LCs$PLC)
      second_stage_data <- rbind(LCs$FLC)
    }
        
    first_stage_labels <- kmeanspp(first_stage_data, first_stage_nstates, 
                                   iter.max = 100, nstart = 10)$cluster

    second_stage_labels <- rep(NA, nsamples)
    for (ll in 1:(first_stage_nstates-1)) {
      second_stage_labels[first_stage_labels == ll] <- (ll - 1) * first_stage_nstates + 
        kmeanspp(second_stage_data[first_stage_labels == ll,], 
                 second_stage_nstates, 
                 iter.max = 100, nstart = 10)$cluster
    }
    second_stage_labels[first_stage_labels == first_stage_nstates] = (first_stage_nstates - 1) * first_stage_nstates + 
      kmeanspp(second_stage_data[first_stage_labels == first_stage_nstates,], 
               second_stage_last_run_nstates, 
               iter.max = 100, nstart = 10)$cluster
    
    state_vector <- second_stage_labels
  }
  invisible(state_vector)
}

