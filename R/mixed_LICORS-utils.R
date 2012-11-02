#' @title Utilities for ``LICORS'' class
#' @name mixed_LICORS-utils
#' @aliases plot.mixed_LICORS summary.mixed_LICORS
#' @description 
#' 
#' The \code{"mixed_LICORS"} class is the output from the 
#' \code{\link{mixed_LICORS}} estimator.
#'
NULL

#' @rdname mixed_LICORS-utils
#'
#' @description 
#' \code{plot.mixed_LICORS} gives a visual summary of the estimates
#' such as marginal state probabilities, conditional state probabilities 
#' (= weight matrix), predictive state densities, trace plots for 
#' log-likelihood/loss/penalty.
#'
#' @param x object of class \code{"mixed_LICORS"} 
#' @param type should only \code{"training"}, \code{"test"}, or \code{"both"} be
#' plotted. Default: \code{"both"}.
#' @param cex.axis The magnification to be used for axis annotation relative to 
#' the current setting of \code{cex}.
#' @param cex.lab The magnification to be used for x and y labels 
#' relative to the current setting of \code{cex}.
#' @param cex.main The magnification to be used for main titles relative to the 
#' current setting of \code{cex}.
#' @param line on which margin line should the labels be ploted, 
#' starting at 0 counting outwards (see also \code{\link[graphics]{mtext}}).
#' @param ... optional arguments passed to \code{plot} or to \code{summary} 
#' @keywords hplot
#' @import fields
#' @method plot mixed_LICORS
#' @export
#' @examples
#' # see examples of LICORS-package
#'

plot.mixed_LICORS  =  function(x, type = "both", cex.axis = 1.5, cex.lab = 1.5, 
                         cex.main = 2, line = 1.5, ...){
  out = x
  op = par(no.readonly=TRUE)
  dev.off()
  if (type == "train"){
    state_probs_matrix = out$conditional_state_probs[out$sel_train,]
    #state_probs = out$pi
    loglik = out$loglik_trace[, c("train", "train weighted")][1:out$niter['final'],]
    mse = out$MSE_trace[, c("train", "train weighted")][1:out$niter['final'],]
  } else if (type == "test"){
    state_probs_matrix = out$conditional_state_probs[-out$sel_train,]
    loglik = out$loglik_trace[, c("test", "test weighted")][1:out$niter['final'],]
    mse = out$MSE_trace[, c("test", "test weighted")][1:out$niter['final'],]
  } else {
    state_probs_matrix = out$conditional_state_probs
    #state_probs = LICORS_pi_eps(out$theta)
    loglik = out$loglik_trace[1:out$niter['final'],]
    mse = out$MSE_trace[1:out$niter['final'],]
  }
  
  ordered_states = order(colSums(state_probs_matrix), decreasing = TRUE)
  
  state_probs_matrix = state_probs_matrix[, ordered_states]
  state_probs = estimate_state_probs(state_probs_matrix)
  names(state_probs) = paste(1:length(state_probs))
  names(state_probs) = NULL
  
  nPLCs = nrow(state_probs_matrix)
  
  states_over_time = rep(NA, out$niter['final'])
  states_over_time[1] = out$nstates["start"]
  
  for (ii in 1:out$nmerges){
    states_over_time[out$merging_times[ii]] = out$nstates["start"] - ii
  }  
  
  FLCs_train_order = order(rbind(out$LCs[["FLC"]])[out$sel_train, ])
  FLCs_train_ordered = out$LCs[["FLC"]][out$sel_train,][FLCs_train_order]

  FLC_pdfs_train = estimate_LC_pdfs(LCs = cbind(out$LCs[["FLC"]][out$sel_train,]), 
                                    weight_matrix = out$conditional_state_probs[out$sel_train, ][, ordered_states],
                                    state_vector = NULL, #state_vector_train,
                                    method=out$method[["FLC"]])
  
  FLC_pdfs_train = FLC_pdfs_train[FLCs_train_order,]
  
  image_colors = two.colors(n = 100, "darkblue", "red", "green")
  #image_colors = tim_colors(100)
  #  image_colors = colorRampPalette(brewer.pal(9, name="YlOrRd"))(100)
  #  image_colors = colorRampPalette(brewer.pal(9, name="RdGy"))(100)
  
  #  image_colors = rev(colorRampPalette(brewer.pal(9, name="RdBuYl"))(100))
  #state_colors = two_colors(out$nstates["best"], "darkgreen", "darkred", "gray")
  state_colors = colorRampPalette(brewer.pal(9, name="Set1"))(out$nstates['best'])
  
  state_tickmarks = pretty(1:out$nstates["best"])
  state_probs_tickmarks = pretty(state_probs)
  LC_tickmarks = pretty(1:nPLCs, 5)
  iter_tickmarks = pretty(1:out$niter['final'], out$niter['final'])
  FLC_tickmarks = pretty(out$FLCs_train[FLCs_train_ordered])
  pdf_FLC_tickmarks = pretty(FLC_pdfs_train)
  ###########################
  ### Start the plot
  ###########################
  
  par(mfrow = c(1,1), mar = c(5,5,2,5), cex.axis = cex.axis, cex.lab = cex.lab)
  layout(matrix(c(1,2,3, 3, 4,5,6, 7), byrow=FALSE, ncol = 2))
  # conditional FLC distributions
  par(mar = c(1,4,1,1))
  matplot(FLCs_train_ordered, FLC_pdfs_train, type="l", lwd=2, lty = 1,
          col = state_colors, axes = FALSE, ylab = "", xlab = "", main ="", cex.main = cex.main)
  box() 
  axis(1, at = FLC_tickmarks, labels = paste(FLC_tickmarks))
  axis(4, at = pdf_FLC_tickmarks, labels=paste(pdf_FLC_tickmarks))
  
  abline(v = FLCs_train_ordered[apply(FLC_pdfs_train, 2, which.max)], col=state_colors, lty =2)
  #mtext("FLC", 1, line = line, cex = cex.lab)
  mtext("p(x|S)", 2, line = line, cex = cex.lab)
  legend("topright", "training only")
  
  
  # distribution over states
  par(mar = c(1,4,1,1))
  barplot(state_probs, space = 0, main="" , xlab="", ylab="", col = state_colors, 
          cex.axis = cex.axis, cex.main = cex.main, ylim = c(0, max(state_probs)*1.05),
          axes = FALSE, xlim = c(0+0.25, out$nstates["best"]-0.25))
  if (type == "both"){
    legend("topright", "training & test")
  } else {
    legend("topright", type)
  }
  
  #mtext("state id", 1, line = line, cex = cex.lab)
  mtext("P(S)", 2, line = line, cex = cex.lab)
  box()
  #axis(1, at = state_tickmarks, labels = paste(state_tickmarks), col = state_colors)
  axis(4, at = state_probs_tickmarks, labels=paste(state_probs_tickmarks))
  
  par(mar = c(5,4,0.1,1))
  #image(1:out$nstates["best"], 1:nPLCs, t(state_probs_matrix), axes=FALSE, 
  #      col=image_colors, xlab = "", ylab = "", main = "", cex.main = cex.main)
  image2(state_probs_matrix, axes=FALSE, legend = FALSE,
        col=image_colors, xlab = "", ylab = "", main = "", cex.main = cex.main)
  abline(v = 0:out$nstates["best"] + 0.5, col = "gray", lwd = 1)
  if (type == "both"){
    abline(h = length(out$theta) - length(-out$sel_train), col = "white", lwd = 2)
  }
  box(col = "black")
  mtext("state id", 1, line = line+2, cex = cex.lab)
  mtext("P(S|x, PLC)", 2, line  = line, cex=cex.lab)
  mtext("PLC id", 4, line = line+2, cex = cex.lab, at = nrow(state_probs_matrix)/3)
  
  if (type == "both"){
    mtext("test", 2, line = 0.25, at = (length(out$theta) - length(-out$sel_train))*0.5)
    mtext("training", 2, line = 0.25, at = nrow(state_probs_matrix)*4/5 )
  }
  
  
  #mtext(expression("P(S| x, L"^-"-glyphosate line1"), 4, line = line, cex = cex.lab)
  #image.plot(t(state_probs_matrix), axes=F, col=image_colors, horizontal=F, legend.only=T, legend.width=.25, legend.shrink=T, add=T, axis.args=list(lwd=0, mgp=c(2,.5,0), tck=-.5), smallplot= c(.91,.94, 0.1,.99))
  
  axis(1, at = state_tickmarks, labels = paste(state_tickmarks), col = state_colors)
  axis(4, at = LC_tickmarks, labels=rev(paste(LC_tickmarks)))
  
  
  # Right column of the plot
  par(mar = c(0,5,1,3))
  if (type != "both"){
    #matplot(loglik, ylab = "", main = "MSE", type = "l", lty = c(1,2,1,2), col = c(1,1,2,2), lwd = 2, axes = FALSE)
    #legend("bottomleft", colnames(out$MSE_trace), lty = c(1,2,1,2), col=c(1,1,2,2), lwd=2)
    #axis(1, at = iter_tickmarks, labels = paste(iter_tickmarks) )
    #axis(2)
    #box()
    
    #abline(v = out$niter["best"], lty = 1)
    #mtext("Iteration", 1, line = line, cex = cex.lab)
    #mtext("MSE", 2, line = line, cex = cex.lab)
    
    ts.plot(ts(loglik), pch = 1, main = "", ylab = "", xlab = "", axes = FALSE)
    abline(v = c(0, out$merging_times-0.5))
    axis(2)
    axis(1, at = iter_tickmarks, labels = paste(iter_tickmarks) )
    box()
    mtext("Iteration", 1, line = line, cex = cex.lab)
    mtext(paste("Log-likelihood (", type, ")", sep=""), 2, line = line, cex = cex.lab)
    
    par(new=T)
    plot(states_over_time, pch = 19, xlab = "", ylab = "", axes = FALSE)
    states_over_time = na.locf(states_over_time)
    lines(states_over_time, lwd=2)
    axis(2, at = pretty(out$nstates["start"]:out$nstates["final"], n = out$nmerges))
    mtext("Number of states", 4, line=2, cex=cex.lab)
    
  } else {
    matplot(loglik, ylab = "", main = "", type = "l", lty = c(1,2,1,2), col = c(1,1,2,2), lwd = 2, axes = FALSE, cex.main = cex.main)
    #axis(1, at = iter_tickmarks, labels = paste(iter_tickmarks) )
    box()
    axis(4)
    abline(v = out$niter["best"], lty = 3, lwd = 2)
    #mtext("Iteration", 1, line = line, cex = cex.lab)
    mtext("log-likelihood", 2, line = line, cex = cex.lab)
  }
  
  par(mar = c(0,5,0,3))
  # MSE comparison
  matplot(mse, ylab = "", main = "", type = "l", lty = c(1,2,1,2), 
          col = c(1,1,2,2), lwd = 2, axes = FALSE, cex.main = cex.main)
  #legend("bottomleft", colnames(out$MSE_trace), lty = c(1,2,1,2), col=c(1,1,2,2), lwd=2)
  box()  
  axis(4)
  
  abline(v = out$niter["best"], lty = 3, lwd = 2)
  #mtext("Iteration", 1, line = line, cex = cex.lab)
  mtext("MSE", 2, line = line, cex = cex.lab)
  
  par(mar = c(4,5,0,3))
  matplot(out$penalty[1:out$niter['final'], "entropy"], type = "l", lty = 1:2, col=1:2, 
          lwd=2,main = "", ylab = "", cex.main = cex.main,
          axes = FALSE)
  box()
  axis(1, at = iter_tickmarks, labels = paste(iter_tickmarks) )
  axis(4)
  
  #legend("topright", colnames(out$penalty), lty = 1:2, col = 1:2, lwd=2, cex = cex.lab)
  abline(v = out$niter["best"], lty = 3, lwd = 2)
  mtext("Iteration", 1, line = line+1, cex = cex.lab)  
  mtext("penalty", 2, line = line, cex = cex.lab)
  
  par(mar = c(1,7,1,3))
  plot.new()
  legend("topright", colnames(mse), lty = c(1,2,1,2), col=c(1,1,2,2), lwd=2, cex = cex.lab)
  plot.window(xlim=c(0,1), ylim=c(0,1))
  rect(0, seq(0, 1, length=100)[-100],
       .1, seq(0, 1, length=100)[-1],
       col=image_colors, border=NA)
  #text(0.1, 0.5, "probability", las = 1)
  
  #box()
  axis(2, at=pretty(seq(0, 1, length=100)), las=0)
  mtext("probability", side = 2, cex = cex.lab*0.75, line = 2, at = 0.25)
  
  par(op)
}

#' @rdname mixed_LICORS-utils
#' @description 
#' \code{summary.mixed_LICORS} prints out a summary of the estimated LICORS model.
#'
#' @param object object of class \code{"mixed_LICORS"} 
#' @keywords model nonparametric
#' @export
#' @method summary mixed_LICORS
#' @examples
#' # see examples in LICORS-package
#'

summary.mixed_LICORS  =  function(object, ...){
  
  cat(rep("*", 20))
  cat("\n")
  
  if (object$converged){
    cat("Mixed LICORS converged after", object$niter["final"], "iterations
        (", object$nstates["final"], "states).")
  } else {
    cat("Mixed LICORS did NOT converge. \n")
    cat("Check trace plots to see if the solution is an interior optimum (vertical dashed bar in trace plots).\n \n")
  }
  cat("\n")
  cat("The best model W* was achieved at iteration", object$niter["best"], 
      " with", object$nstates["start"], "remaining states (starting from", object$nstates["start"],"). \n")
  
  
  entropies = rep(NA,3)
  names(entropies) = c("both", "training", "test")
  entropies["both"] = compute_mixture_penalty(object$conditional_state_probs, type = "entropy", base = "nstates")
  entropies["training"] = compute_mixture_penalty(object$conditional_state_probs[object$sel_train,], type = "entropy", base = "nstates")
  entropies["test"] = compute_mixture_penalty(object$conditional_state_probs[-object$sel_train,], type = "entropy", base = "nstates")
  cat("\n")  
  cat("The average entropy penalty for W* equals", round(object$penalty[object$niter["best"],"entropy"]*100, 1), "%, where 0% is unique cluster assignment, and 100% is uniformly at random.\n\n") 
  cat("Entropy for both vs training vs test (in %): \n")
  cat(round(cbind(entropies*100),1))
  cat("\n")
  cat(rep("*", 20))
}




