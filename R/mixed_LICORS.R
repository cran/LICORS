#' @title Mixed LICORS: An EM-like Algorithm for Predictive State Space Estimation
#'
#' @description 
#' \code{mixed_LICORS} is the core function of this package as it estimates the
#' ``parameters'' in the model for the spatio-temporal process.
#' \deqn{
#' P(X_1, \ldots, X_{\tilde{N}}) \propto \prod_{i=1}^{N} P(X_i \mid \ell^{-}_i) 
#' =  \prod_{i=1}^{N} P(X_i \mid \epsilon(\ell^{-}_i)) .
#' }
#' 
#' @param LCs list of PLCs and FLCs matrices (see output of 
#' \code{\link{data2LCs}} for details and formatting).
#' @param nstates_start number of states to start the EM algorithm
#' @param initialization a a) character string, b) vector, or c) matrix. a) 
#' results \code{nstates_start} many states initialized by passing the character
#'  string as \code{method} argument of \code{\link{initialize_states}}; 
#' if b) the vector will be taken as initial state labels;
#' if c) the matrix will be taken as initial weights. Note that for both b) and c)
#' \code{nstates_start} will be ignored.
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param max_iter maximum number of iterations in the EM
#' @param verbose logical; if \code{TRUE} it prints output in the console
#' as the EM is running
#' @param sparsity what type of sparsity (currently not implemented)
#' @param lambda penalization parameter; larger lambda gives sparser weights
#' @param loss an R function specifying the loss for cross-validation (CV). 
#' Default: mean squared error (MSE), i.e. 
#' \code{loss = function(x, xhat) mean((x-xhat)^2)}
#' @param alpha significance level to stop testing. Default: \code{alpha=0.01}
#' @param seed set seed for reproducibility. Default: \code{NULL}. If 
#' \code{NULL} it sets a random seed and then returns this seed in the output.
#' @param CV_train_ratio how much of the data should be training data.
#' Default: \code{0.75}, i.e. \eqn{75\%} of data is for training
#' @param CV_split_random indicator for splitting data randomly in training and
#'  test data (\code{TRUE}) or to use first part (in time) as training, rest as 
#'  test (\code{FALSE}; default).
#' @param method a list of length \eqn{2} with arguments \code{PLC} and 
#' \code{FLC} for the method of density estimation in each 
#' (either \code{"normal"} or \code{"nonparametric"}).
#' @return
#' An object of class \code{"LICORS"}.
#' @keywords nonparametric cluster multivariate distribution
#' @export
#' @seealso \code{\link{plot.mixed_LICORS}}, \code{\link{summary.mixed_LICORS}}
#' @examples
#' \dontrun{
#' data(contCA00)
#'  
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 2, FLC = 0), shape ="cone")
#' bb = data2LCs(t(contCA00$observed), LC_coords = LC_geom$coords)
#' 
#' mm = mixed_LICORS(bb, nstates_start = 10, init = "KmeansPLC", max_iter = 20)
#' plot(mm)
#' ff_new = estimate_LC_pdfs(bb$FLC, weight_matrix = mm$conditional_state_probs, 
#'                           method = "nonparametric")
#' matplot(bb$FLC, ff_new, pch = ".", cex = 2)
#'}
#'

mixed_LICORS <-
function(LCs = list(PLC = PLCs, FLC = FLCs, dim = list(original = NULL, truncated = NULL)),
         nstates_start=NULL,
         initialization=NULL, max_iter=500, 
         method=list(PLC = "normal", FLC = "nonparametric"), 
         alpha=0.01, 
         seed=NULL, verbose=TRUE, CV_train_ratio = 0.75, 
         CV_split_random = FALSE, 
         sparsity = "stochastic", 
         lambda = 0,
         loss = function(x, xhat) mean((x-xhat)^2)){
  
  FLCs = LCs$FLC
  PLCs = LCs$PLC
  NN = nrow(PLCs)
  
  if (is.null(LCs$dim)){
    LCs$dim = list(original = NULL, truncated = c(1, NN))  
  }
  
  
  if (is.null(seed)){
    seed = sample.int(10^6, 1)
  }
  set.seed(seed)
  
  if (is.null(initialization) & is.null(nstates_start)){
    stop("You must provide either 
         - an integer specifying the number of states
         - a vector with the initial state space labels
         - or a matrix with the initial weights_
         ")
  }
  state_vector = rep(NA, NN)
  
  if (CV_split_random){
    sel_train = sample.int(NN, size = floor(CV_train_ratio * NN), replace = FALSE)
  } else {
    sel_train = 1:(floor(CV_train_ratio * NN))
  }
  
  
  ntrain = length(sel_train)
  ntest = NN-ntrain
  
  PLC_train = PLCs[sel_train,]
  FLC_train = cbind(FLCs[sel_train,])
  
  PLC_test = PLCs[-sel_train,]
  FLC_test = cbind(FLCs[-sel_train,])
    
  #####################################################
  #### initialize variables to keep track of EM updates
  #####################################################
  
  # mixture weight penalty
  penalty = matrix(NA, nrow = max_iter, ncol = 2)
  colnames(penalty)  = c("1-L2", "entropy")
  
  loglik_train = rep(-Inf, max_iter)
  loglik_train_weighted = rep(-Inf, max_iter)
  loglik_train_weighted_only_update = rep(-Inf, max_iter)
  
  loglik_test = rep(-Inf, max_iter)
  loglik_test_weighted = rep(-Inf, max_iter)
  
  MSE_train = rep(NA, max_iter)
  MSE_train_weighted = rep(NA, max_iter)
  MSE_test = rep(NA, max_iter)
  MSE_test_weighted = rep(NA, max_iter)
  
  # negative entropy criterion (NEC) for choosing number of clusters
  NEC = rep(NA, max_iter)
  NEC_weighted = rep(NA, max_iter)  
  
  weights_train_temp = NULL
  ##########################
  #### Initialize states ###
  ##########################
  if (is.character(initialization)) {   
    if (is.null(nstates_start)){
      stop("You must provide the number of maximum states to start with.")
    } else {
      # assign states by kmeans or randomly
      state_vector[sel_train] <- initialize_states(nstates = nstates_start, 
                                                   nsamples = ntrain, 
                                                   method = initialization,
                                                   LCs = list(PLCs = PLC_train,
                                                              FLCs = FLC_train))
                                                   
      # state vector initialization for test data is irrelevant; but need
      # the full vector later
      state_vector[-sel_train] = initialize_states(nstates = nstates_start,
                                                   nsamples = ntest,
                                                   method = "random")
    }
  } else {
    if (is.vector(initialization)){
      state_vector <- initialization
    } else if (is.matrix(initialization)){
      weights_train_temp <- initialization[sel_train, ]
      state_vector <- weight_matrix2states(weights_train_temp)
    } else {
      stop("Initial states must be either a 
           i) character string describing the initialization method,
           ii) a vector with the labels (length = N)
           iii) or a weight matrix (dim = N x K)")
    }
    
  }
  state_vector_train = state_vector[sel_train] 
  state_vector_test = state_vector[-sel_train]
  
  if (is.null(weights_train_temp)){
    weights_train_temp = states2weight_matrix(state_vector = state_vector_train, 
                                              nstates_total = nstates_start)
  }
  nstates_start = max(state_vector)

  FLC_pdfs_train_one_cluster = estimate_LC_pdf_given_state(state_label = 1, 
                                                           LCs=FLC_train, 
                                                           state_vector = rep(1, nrow(FLC_train)),
                                                           method = method[["FLC"]])
  
  loglik_one_cluster = compute_LICORS_loglik(weight_matrix=1, 
                                             FLC_pdfs = cbind(FLC_pdfs_train_one_cluster))
  

  
  penalty[1,"1-L2"] = compute_mixture_penalty(weights_train_temp, "1-Lq", q = 2)
  penalty[1,"entropy"]= compute_mixture_penalty(weights_train_temp, "entropy", base = "nstates")
  
  nstates_start = ncol(weights_train_temp)

  FLC_pdfs_train = estimate_LC_pdfs(LCs = FLC_train, 
                                    weight_matrix = weights_train_temp,
                                    state_vector = NULL, #state_vector_train,
                                    method=method[["FLC"]])

  PLC_pdfs_train = estimate_LC_pdfs(LCs = PLC_train, 
                                    weight_matrix = weights_train_temp,
                                    state_vector = NULL, #state_vector_train,
                                    method=method[["PLC"]])
  
  predictions_train = predict_FLC_given_PLC(state_vector = NULL, 
                                            FLC_train = FLC_train, 
                                            FLCs_train_pdf = NULL,
                                            type = "mean", 
                                            weight_matrix_train = weights_train_temp)
  
  predictions_train_weighted = predict_FLC_given_PLC(state_vector = NULL, 
                                                     FLC_train = FLC_train,
                                                     FLCs_train_pdf = NULL,
                                                     type = "weightedmean", 
                                                     weight_matrix_train = weights_train_temp)
  
  MSE_train[1] = loss(FLC_train, predictions_train)
  MSE_train_weighted[1] = loss(FLC_train, predictions_train_weighted)
  
  PLC_pdfs_test = estimate_LC_pdfs(LCs = PLC_train, 
                                   weight_matrix = weights_train_temp,
                                   state_vector = NULL,
                                   method=method[["PLC"]],
                                   eval_LCs = PLC_test)

  weights_test = estimate_state_probs(weight_matrix = weights_train_temp,
                                      state_vector = NULL, 
                                      PLC_pdfs = PLC_pdfs_test, 
                                      FLC_pdfs = NULL,
                                      nstates_total = nstates_start)

  state_vector_test = weight_matrix2states(weights_test)
  
  predictions_test_weighted = predict_FLC_given_PLC(state_vector = NULL, 
                                                    FLC_train = FLC_train,
                                                    FLCs_train_pdf = NULL,
                                                    type = "weightedmean", 
                                                    weight_matrix_train = weights_train_temp,
                                                    weight_matrix_test = weights_test)
    
  
  predictions_test = predict_FLC_given_PLC(state_vector = NULL, 
                                           FLC_train = FLC_train,
                                           FLCs_train_pdf = NULL,
                                           type = "mean", 
                                           weight_matrix_train = weights_train_temp,
                                           weight_matrix_test = weights_test)
  
  MSE_test[1] = loss(FLC_test, predictions_test)
  MSE_test_weighted[1] = loss(FLC_test, predictions_test_weighted)
  
  loglik_train_weighted[1] = compute_LICORS_loglik(FLC_pdfs_train, weights_train_temp, lambda = lambda)
  loglik_train[1] = loglik_train_weighted[1]
  loglik_train_weighted_only_update[1] = loglik_train_weighted[1]
  
  FLC_pdfs_test = estimate_LC_pdfs(LCs = FLC_train, 
                                   weight_matrix = weights_train_temp,
                                   state_vector = NULL,
                                   method=method[["FLC"]],
                                   eval_LCs = FLC_test)
  
  loglik_test_weighted[1] = compute_LICORS_loglik(FLC_pdfs_test, weights_test, lambda = lambda)
  loglik_test[1] = compute_LICORS_loglik(FLC_pdfs_test, 
                                         states2weight_matrix(state_vector_test, 
                                                              nstates_total = nstates_start), 
                                         lambda = lambda)
  
  weights_final = weights_train_temp
  prediction_final = predictions_test
  nstates_final = nstates_start
  
  weights_best = weights_train_temp
  prediction_best = predictions_test
  prediction_best_weighted = predictions_test
  weights_best_test = weights_test
  iter_best = 1
  stop.iterations = FALSE
  temp_converged = FALSE
  overall_converged = FALSE
  
  merging_times = c()

  for (ii in 2:max_iter){
    merged = FALSE
    temp_merge = FALSE
    temp_more.sparse = TRUE
    
    iter_final = ii
    nstates_temp = ncol(weights_train_temp)

    effective_nsamples = colSums(weights_train_temp)
    min_nsamples_per_state = ncol(PLC_train) + 1
    
    # E step  
    weights_train_temp = estimate_state_probs(weight_matrix = weights_train_temp,
                                              state_vector = NULL, 
                                              PLC_pdfs = PLC_pdfs_train, 
                                              FLC_pdfs = FLC_pdfs_train,
                                              nstates_total = ncol(FLC_pdfs_train))

    if (any("deterministic" == sparsity) && lambda != 0){
      if (!merged && ii > 1000){#} && ii < max_iter){
        if (verbose){
          cat("sparsity enforced! \n")
          cat("Penalty before:", round(compute_mixture_penalty(weights_train_temp, "entropy", base = "nstates")*100,4), "% \n")
        }
        #weights_train_temp = LICORS_adjust_weights(weights_final, weights_train_temp, lambda = lambda)
        weights_train_temp = sparsify_weights(weights_train_temp,  NULL, lambda = lambda)
  	    if (verbose){
          cat("Penalty after:", round(compute_mixture_penalty(weights_train_temp, "entropy", base = "nstates")*100,4), "% \n")
        }
      }
    }

    if (merged) {
      PLC_pdfs_train = estimate_LC_pdfs(LCs = PLC_train, 
                                        weight_matrix = weights_train_temp,
                                        state_vector = NULL, #state_vector_train,
                                        method=method[["PLC"]])
      
      FLC_pdfs_train = estimate_LC_pdfs(LCs = FLC_train, 
                                        weight_matrix = weights_train_temp,
                                        state_vector = NULL, #state_vector_train,
                                        method=method[["FLC"]])
        
      # M step
      weights_train_temp = estimate_state_probs(weight_matrix = weights_train_temp,
                                                state_vector = NULL, 
                                                PLC_pdfs = PLC_pdfs_train, 
                                                FLC_pdfs = FLC_pdfs_train,
                                                nstates_total = ncol(FLC_pdfs_train))
    }
    
    if (!merged){
      # maximum weight difference
      max_weight_diff = max( sqrt( rowMeans( (weights_train_temp - weights_final)^2 ) ) )
    } else {
      max_weight_diff = 1
    }

    if (max_weight_diff < 10^(-2)){
      temp_converged = TRUE
    }
    
    if (temp_converged) {
      temp_merge = TRUE
      if (any("deterministic" == sparsity) && lambda != 0){
        temp_merge = FALSE
        temp_old_penalty = penalty[ii-1,"entropy"]
        if (verbose){
          cat("converged and sparsity enforced! \n")
          cat("Penalty before:", round(temp_old_penalty*100,4), "% \n")
        }
        # sparsity
        weights_train_temp = sparsify_weights(weights_train_temp, NULL, lambda=lambda)        
        temp_new_penalty = compute_mixture_penalty(weights_train_temp, "entropy")
        
        if (temp_old_penalty*0.999 > temp_new_penalty){
          temp_more.sparse = TRUE
        } else {
          temp_more.sparse = FALSE
        }    
        
        if (verbose){
          cat("Penalty after:", round(temp_new_penalty*100,4), "% \n")
        }
        
        max_weight_diff = max( sqrt( rowMeans( (weights_train_temp - weights_final)^2 ) ) )
        
        if (max_weight_diff < 10^(-2) 
            || ( any( abs(loglik_train_weighted[ii-1]-loglik_train_weighted[ii - 1:min(ii-1,10)]) < 10^(-6)) && loglik_train_weighted[ii-1] > -Inf)) { 
          temp_converged = TRUE
        }
      }
    }
    
    
    if (temp_more.sparse == FALSE){
      # if new sparser weights are actually not sparse, then merge clusters
      temp_merge = TRUE
    }
    if (verbose) {    
      cat("Solution more sparse", temp_more.sparse, "\n")
      cat("EM converged temporarily", temp_converged, "\n")
      cat("Should states be merged", temp_merge, "\n")
    }
    loglik_train_weighted[ii] = compute_LICORS_loglik(FLC_pdfs_train, weights_train_temp, lambda = lambda)
    loglik_train[ii] = compute_LICORS_loglik(FLC_pdfs_train, 
                                             states2weight_matrix(weight_matrix2states(weights_train_temp), 
                                                                nstates_total = nstates_temp),
                                             lambda = lambda )
    
    penalty[ii,"1-L2"]= compute_mixture_penalty(weights_train_temp, "1-Lq", q = 2)
    penalty[ii,"entropy"]= compute_mixture_penalty(weights_train_temp, "entropy", base = "nstates")    
    
    NEC[ii] = log(nstates_temp) * penalty[ii, "entropy"] / (loglik_train[ii] - loglik_one_cluster)
    NEC_weighted[ii] =  log(nstates_temp) * penalty[ii, "entropy"] / (loglik_train_weighted[ii] - loglik_one_cluster)
    
    if (temp_merge) {
        AA = estimate_state_adj_matrix(FLC_pdfs = FLC_pdfs_train, 
                                          alpha = NULL, 
                                          distance = function(f,g) {
                                            return(mean(abs(f-g))) })
        
        diag(AA) = 0 # for norm and testing (similarity)
        # if EM converged, and no merging is possible, then stop the iterations
        if (all(AA < alpha)){
          #if (TRUE){
          stop.iterations = TRUE
        } else {
          merged = TRUE
          merging_times = c(merging_times, ii)
          weights_train_temp = merge_states(which(AA == max(AA), arr.ind = TRUE)[1:2],
                                                weights_train_temp)
          
          temp_converged = FALSE
          if (verbose){
            cat("Merged two states into one after convergence with", nstates_temp ," states. \n")
          }
        }
      
    }
    
    PLC_pdfs_train = estimate_LC_pdfs(LCs = PLC_train, 
                                      weight_matrix = weights_train_temp,
                                      state_vector = NULL, #state_vector_train,
                                      method=method[["PLC"]])
    
    FLC_pdfs_train = estimate_LC_pdfs(LCs = FLC_train, 
                                      weight_matrix = weights_train_temp,
                                      state_vector = NULL, #state_vector_train,
                                      method=method[["FLC"]])    
    

    loglik_train_weighted_only_update[ii] = compute_LICORS_loglik(FLC_pdfs_train, 
                                                                  weights_train_temp, 
                                                                  lambda = lambda)

    nstates_temp = ncol(weights_train_temp)
    
    
    # evaluate in-sample MSE    
    predictions_train = predict_FLC_given_PLC(state_vector = NULL, 
                                              FLC_train = FLC_train, 
                                              FLCs_train_pdf = NULL,
                                              type = "mean", 
                                              weight_matrix_train = weights_train_temp)
    
    predictions_train_weighted = predict_FLC_given_PLC(state_vector = NULL, 
                                                       FLC_train = FLC_train,
                                                       FLCs_train_pdf = NULL,
                                                       type = "weightedmean", 
                                                       weight_matrix_train = weights_train_temp)
  
    
    MSE_train_weighted[ii] = loss(FLC_train, predictions_train_weighted)
    MSE_train[ii] =  loss(FLC_train, predictions_train)
    
    # evaluate out-of-sample MSE
    PLC_pdfs_test = estimate_LC_pdfs(LCs = PLC_train, 
                                     weight_matrix = weights_train_temp,
                                     state_vector = NULL,
                                     method=method[["PLC"]],
                                     eval_LCs = PLC_test)
    
    weights_test = estimate_state_probs(weight_matrix = weights_train_temp,
                                        state_vector = NULL, 
                                        PLC_pdfs = PLC_pdfs_test, 
                                        FLC_pdfs = NULL,
                                        nstates_total = ncol(PLC_pdfs_test))
    
    state_vector_test <- weight_matrix2states(weights_test)
    
    predictions_test_weighted = predict_FLC_given_PLC(state_vector = NULL, 
                                                      FLC_train = FLC_train,
                                                      FLCs_train_pdf = NULL,
                                                      type = "weightedmean", 
                                                      weight_matrix_train = weights_train_temp,
                                                      weight_matrix_test = weights_test)
    
    
    predictions_test = predict_FLC_given_PLC(state_vector = NULL, 
                                             FLC_train = FLC_train,
                                             FLCs_train_pdf = NULL,
                                             type = "mean", 
                                             weight_matrix_train = weights_train_temp,
                                             weight_matrix_test = weights_test)

    MSE_test_weighted[ii] =  loss(FLC_test, predictions_test_weighted)
    MSE_test[ii] = loss(FLC_test, predictions_test)

    FLC_pdfs_test = estimate_LC_pdfs(LCs = FLC_train, 
                                     weight_matrix = weights_train_temp,
                                     state_vector = NULL,
                                     method=method[["FLC"]],
                                     eval_LCs = FLC_test)
    
    loglik_test_weighted[1] = compute_LICORS_loglik(FLC_pdfs_test, weights_test, lambda = lambda)
    loglik_test[1] = compute_LICORS_loglik(FLC_pdfs_test, 
                                           states2weight_matrix(state_vector_test, 
                                                                nstates_total = ncol(FLC_pdfs_test)), 
                                           lambda = lambda)
    
    if (MSE_test_weighted[ii] < min(MSE_test_weighted[-ii], na.rm=TRUE)){
      weights_best = weights_train_temp
      prediction_best = predictions_test
      prediction_best_weighted = predictions_test_weighted
      weights_best_test = weights_test
      iter_best = ii
    }

    if (verbose){
      #cat(table(state_vector_train))
      cat("\n Finished step", ii, "with", nstates_temp, "states. \n")
      cat("Loglik:", loglik_train_weighted[ii], "\n")
      cat("Out-of-sample loglik:", loglik_test_weighted[ii], "\n" )
      cat("Penalty:", round(penalty[ii,]*100,1), "% \n")
      if (!merged){
        cat("Maximum difference in weights:", round(max_weight_diff, 4), "\n")
      } else {
  	    cat("Weights have been merged \n")
      }
      cat("In-sample MSE:", MSE_train[ii], "\n" )
      cat("Weighted In-sample MSE:", MSE_train_weighted[ii], "\n" )
      cat("Out-of-sample MSE:", MSE_test[ii], "\n" )
      cat("Weighted Out-of-sample MSE:", MSE_test_weighted[ii], "\n \n" )
      cat("In sample 'R^2':", round( 100* (1- MSE_train_weighted[ii] / var(FLC_train)),1), "%\n \n")
    }
    
    weights_final = weights_train_temp
   
    if (ncol(weights_train_temp) < 3) {
      stop.iterations = TRUE
    }

    if (stop.iterations){
      overall_converged = TRUE
      if (verbose) {
        cat("***************************************************** \n")
        cat("EM algorithm converged to (local) optimum at iteration", iter_final,".\n")
        cat("***************************************************** \n \n")
      }
      break
    }
    
  }
  
  if (!stop.iterations) {
    if (verbose) {
  	  cat("***************************************************** \n")
  	  cat("Finished all iterations without convergence. \n")
  	  cat("***************************************************** \n")
    }
  }

  weights_test_final = weights_test
  pred_state_test_final = weight_matrix2states(weights_test_final)
  
  weights_test_best = weights_best_test
  pred_state_test_best = weight_matrix2states(weights_test_best)
  
  FLC_pdfs_best = estimate_LC_pdfs(LCs = FLC_train, 
                                   weight_matrix = weights_best,
                                   state_vector = NULL, #state_vector_train,
                                   method=method[["FLC"]])
  
  iter_final = ii
  
  out = list()
  out$LCs = LCs
  out$sel_train = sel_train
  
  out$loglik_trace_train = loglik_train[1:iter_final]
  out$loglik_trace_test = loglik_test[1:iter_final]
  
  out$loglik_trace_train_weighted = loglik_train_weighted[1:iter_final]
  out$loglik_trace_test_weighted = loglik_test_weighted[1:iter_final]
  
  out$loglik_trace_train_weighted_only_update = loglik_train_weighted_only_update[1:iter_final]
  
  out$loglik_train = loglik_train[c(1, iter_best, iter_final)]
  
  out$loglik_test = loglik_test[c(1, iter_best, iter_final)]
  
  out$loglik_train_weighted = loglik_train_weighted[c(1, iter_best, iter_final)]
  
  out$loglik_test_weighted = loglik_test_weighted[c(1, iter_best, iter_final)]
  
  names(out$loglik_train) = names(out$loglik_test) = c("start", "best", "final")
  names(out$loglik_train_weighted) = names(out$loglik_test_weighted) = names(out$loglik_test)
  
  out$MSE_trace = cbind(MSE_train, MSE_train_weighted, MSE_test, MSE_test_weighted)
  colnames(out$MSE_trace) = c("train", "train weighted", "test", "test weighted")
  
  
  out$loglik_trace = cbind(out$loglik_trace_train, out$loglik_trace_train_weighted, 
                           out$loglik_trace_test, out$loglik_trace_test_weighted)
  colnames(out$loglik_trace) = c("train", "train weighted", "test", "test weighted")
  
  out$theta_train = weight_matrix2states(weights_best)
  
  out$theta = rep(NA, nrow(FLCs))
  out$theta[sel_train] = out$theta_train
  out$theta[-sel_train] = pred_state_test_best
  


  out$theta_mean_train = predict_FLC_given_PLC(state_vector = NULL, 
                                               FLC_train = FLC_train, 
                                               FLCs_train_pdf = NULL,
                                               type = "mean", 
                                               weight_matrix_train = weights_best)
  out$theta_weighted_mean_train = predict_FLC_given_PLC(state_vector = NULL, 
                                                        FLC_train = FLC_train, 
                                                        FLCs_train_pdf = NULL,
                                                        type = "weightedmean", 
                                                        weight_matrix_train = weights_best)  
  out$theta_median_train = predict_FLC_given_PLC(state_vector = NULL, 
                                                        FLC_train = FLC_train, 
                                                        FLCs_train_pdf = NULL,
                                                        type = "median", 
                                                        weight_matrix_train = weights_best)
  
  out$theta_final = rep(NA, nrow(FLCs))
  out$theta_final[-sel_train] = pred_state_test_final
  out$theta_final[sel_train] = weight_matrix2states(weights_final)
  
  out$MSE = c(loss(FLC_train, out$theta_mean_train), min(MSE_test, na.rm = TRUE))
  names(out$MSE) = c("in-sample", "out-of-sample")
  out$MSE_weighted = c(loss(FLC_train, out$theta_weighted_mean_train), min(MSE_test_weighted, na.rm = TRUE))
  names(out$MSE_weighted) = c("in-sample", "out-of-sample")
  
  out$theta_mean = rep(NA, nrow(FLCs))
  out$theta_mean[-sel_train] = prediction_best
  out$theta_mean[sel_train] = out$theta_mean_train  
  
  out$theta_weighted_mean = rep(NA, nrow(FLCs))
  out$theta_weighted_mean[-sel_train] = prediction_best_weighted
  out$theta_weighted_mean[sel_train] = out$theta_weighted_mean_train
  
  
  out$theta_median = rep(NA, nrow(FLCs))
  out$theta_median[-sel_train] = prediction_best
  out$theta_median[sel_train] = out$theta_median_train
  
  out$conditional_state_probs_final = matrix(NA, ncol = ncol(weights_final), 
                                          nrow = nrow(FLCs))
  out$conditional_state_probs_final[sel_train, ] = weights_final
  out$conditional_state_probs_final[-sel_train, ] = weights_test_final
  
  out$conditional_state_probs = matrix(NA, ncol = ncol(weights_best), nrow = nrow(FLCs))
  out$conditional_state_probs[-sel_train, ] = weights_best_test
  out$conditional_state_probs[sel_train, ] = weights_best
  
  out$marginal_state_probs = estimate_state_probs(weights_best)
  out$marginal_state_probs_final = estimate_state_probs(weights_final)
  out$marginal_state_probs_test = estimate_state_probs(weights_best_test)

  out$NEC = NEC
  out$NEC_weighted = NEC_weighted
  out$loglik_one_cluster = loglik_one_cluster

  out$nstates = c(nstates_start, length(out$marginal_state_probs), length(out$marginal_state_probs_final))
  names(out$nstates) = c("start", "best", "final")
  
  out$niter = c(max_iter, iter_best, iter_final)
  names(out$niter) = c("max", "best", "final")
  
  out$penalty = penalty
  
  out$alpha = alpha
  out$merging_times = merging_times
  out$nmerges = length(out$merging_times)
  out$CV_split_random = CV_split_random
  
  out$method = method
  out$sparsity = sparsity
  out$lambda = lambda
  out$seed = seed
  out$converged = overall_converged
  out$dim = LCs$dim
  class(out) = c("LICORS", "mixed_LICORS")
  
  invisible(out)
}
