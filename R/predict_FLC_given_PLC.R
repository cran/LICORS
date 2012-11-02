#' @title Predict FLCs given new PLCs
#'
#' @description 
#' This function predicts FLCs given new PLCs based on the estimated 
#' \eqn{\epsilon} mappings and estimated conditional distributions.
#' 
#' @param FLC_train the matrix of training FLCs. This matrix has \eqn{N} rows, 
#' and \eqn{n_f} columns
#' @param FLCs_train_pdf \eqn{N \times K} matrix with density estimates of 
#' the training FLCs.
#' @param weight_matrix_train \eqn{N \times K} weight matrix from the training 
#' data
#' @param weight_matrix_test \eqn{\tilde{N} \times K} weight matrix from the test data
#' @param state_vector vector of length \eqn{\tilde{N}} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param type type of prediction: \code{'mean'}, \code{'median'}, 
#' \code{'weightedmean'}.
#' @return
#' \eqn{N \times K} matrix
#' @keywords methods
#' @export
#'

predict_FLC_given_PLC <- function(state_vector = NULL, 
                                  FLC_train, FLCs_train_pdf = NULL, 
                                  type = "mean", 
                                  weight_matrix_train = NULL, 
                                  weight_matrix_test = weight_matrix_train) {
  
  if (!is.null(weight_matrix_train)) {
    state_vector_test <- weight_matrix2states(weight_matrix_test)
    state_vector <- weight_matrix2states(weight_matrix_train)
  } else {
    state_vector_test <- state_vector
  }
  
  pred_test <- rep(NA, length(state_vector_test))
  unique_train_states <- sort(unique(state_vector))
  nstates <- length(unique_train_states)
  # if state_vector has lower than weight_matrix_test states (because none of the
  # rows has its maximum in say column j), then adjust weight_matrix_test to this
  # size
  weight_matrix_test <- weight_matrix_test[, unique_train_states]
  
  FLC_train <- cbind(FLC_train)
  pred_states_train <- rep(NA, nstates)
  
  # estimate pred_tests from the training data
  for (ii in 1:nstates) {
    if (type == "mean") {
      pred_states_train[ii] <- mean(FLC_train[state_vector == unique_train_states[ii], 
                                              ])
    } else if (type == "weightedmean") {
      temp_state_weights_norm <- normalize(weight_matrix_train[, ii])
      pred_states_train[ii] <- sum(FLC_train * temp_state_weights_norm)
    } else if (type == "median") {
      pred_states_train[ii] <- median(FLC_train[state_vector == unique_train_states[ii], ])
    }
  }
  if (type == "mode") {
    pred_states_train <- FLC_train[apply(FLCs_train_pdf, 2, which.max)]
  }
  # make pred_tests for test data by weighting with weight matrix
  if (type == "weightedmean") {
    pred_test <- rowSums(sweep(weight_matrix_test, 2, pred_states_train, "*"))
  } else {
    for (ii in 1:nstates) {
      pred_test[state_vector_test == unique_train_states[ii]] <- pred_states_train[ii]
    }
  }
  invisible(pred_test)
} 