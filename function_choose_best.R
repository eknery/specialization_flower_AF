choose_best = function (fit_results){
  # calculate delta_aicc 
  delta_aicc = as.numeric(fit_results$fit_metrics$aicc) - min(as.numeric(fit_results$fit_metrics$aicc))
  fit_results$fit_metrics = data.frame(fit_results$fit_metrics, delta_aicc)
  # find first and second lowest delta aicc
  first_delta = fit_results$fit_metrics[fit_results$fit_metrics$aicc == min(as.numeric(fit_results$fit_metrics$aicc)),]
  minus_first = fit_results$fit_metrics[-which( fit_results$fit_metrics$aicc == min(as.numeric(fit_results$fit_metrics$aicc)) ),]
  second_delta =  minus_first[minus_first$aicc == min(as.numeric(minus_first$aicc)),]
  # index to check complexity
  model_names = fit_results$fit_metrics[,1]
  first_index = which(model_names == first_delta$model)
  second_index = which(model_names == second_delta$model)
  # compare delta aicc and pick best-fit model
  if(first_index > second_index){
    if (second_delta$delta_aicc - first_delta$delta_aicc > 2){
      best_fit_model = first_delta
    } else {
      best_fit_model = second_delta
    }
  } else {
    best_fit_model = first_delta
  }
  # finding best-model estimates
  best_fit_estimates = fit_results$model_estimates[names(fit_results$model_estimates) == best_fit_model$model]
  best_fit_estimates = best_fit_estimates[[1]]
  # pick best-fit estimates
  theta_se = c()
  for(j in 1:nrow(best_fit_estimates$theta) ){
    thse = c(best_fit_estimates$theta[j,])
    names(thse)[1] = paste("theta", as.character(j), sep="_")
    names(thse)[2] = paste("se", as.character(j), sep="_")
    theta_se = c(theta_se, thse)
  }
  alp_sig = c()
  for(j in 1:ncol(best_fit_estimates$solution) ){
    alsi=best_fit_estimates$solution[,j]
    names(alsi)[1] = paste("alpha", as.character(j), sep="_")
    names(alsi)[2] = paste("sigma_sq", as.character(j), sep="_")
    alp_sig = c(alp_sig, alsi)
  }
  best_estimate = c(theta_se,alp_sig)
  # packing things up
  best_fit_metric = best_fit_model[1,]
  best_list = list(best_fit_metric[-5], best_estimate)
  names(best_list) = c("best_fit", "best_estimates")
  return(best_list)
}

