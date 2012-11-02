#' @title Kmeans++
#'
#' @description 
#' \emph{kmeans++} clustering (see References) using R's built-in function 
#' \code{\link[stats]{kmeans}}.
#' 
#' @param data an \eqn{N \times d} matrix, where \eqn{N} are the samples and 
#' \eqn{d} is the dimension of space.
#' @param k number of clusters.
#' @param start first cluster center to start with
#' @param iter.max the maximum number of iterations allowed
#' @param nstart how many random sets should be chosen?
#' @param ... additional arguments passed to \code{\link[stats]{kmeans}}
#' @keywords multivariate cluster
#' @references
#' Arthur, D. and S. Vassilvitskii (2007). ``k-means++: The advantages of careful seeding.''
#'  In H. Gabow (Ed.), Proceedings of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms 
#'  [SODA07], Philadelphia, pp. 1027-1035. Society for Industrial and Applied Mathematics.
#'  
#' @export
#' @seealso \code{\link[stats]{kmeans}}
#' @examples
#' set.seed(1984)
#' nn <- 100
#' XX = matrix(rnorm(nn), ncol = 2)
#' YY = matrix(runif(length(XX)*2, -1, 1), ncol = ncol(XX))
#' ZZ = rbind(XX, YY)
#' 
#' cluster_ZZ = kmeanspp(ZZ, k=5, start = "random")
#' 
#' plot(ZZ, col = cluster_ZZ$cluster+1, pch = 19)
#' 

kmeanspp <- function(data, k = 1, start = "random", iter.max = 100, 
                     nstart = 10, ...) {
   
  kk <- k
  
  data <- cbind(data)
  
  nsamples <- nrow(data)
  ndim <- ncol(data)
  
  out <- list()
  out$tot.withinss <- Inf
  for (starts in 1:nstart) {  
    center_ids <- rep(0, length = kk)
    if (start == "random"){
      center_ids[1:2] = sample.int(nsamples, 1)
    } else if (start == "normal") { 
      center_ids[1:2] = which.min(dmvnorm(data, mean = colMeans(data), sigma = cov(data)))
    } else {
      center_ids[1:2] = start
    }
    
    for (ii in 2:kk){
      if (ndim == 1){
        dists <- apply(cbind(data[center_ids,]), 2, function(center) {rowSums((data - center)^2)})
      } else {
        dists <- apply(data[center_ids,], 1, function(center) {rowSums((data - center)^2)})
      }
      probs <- apply(dists, 1, min)
      probs[center_ids] <- 0
      center_ids[ii] <- sample.int(nsamples, 1, prob = probs)
    }
    
    out_temp <- kmeans(data, centers = data[center_ids,], iter.max = iter.max, ...)
    out_temp$inicial.centers <- data[center_ids,]
    if (out_temp$tot.withinss < out$tot.withinss){
      out <- out_temp
    }
  } 
  invisible(out)
}

