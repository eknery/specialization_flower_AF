# dependencies:
# library(OUwie)

fit_evo_models = function(tree, regimes, models_to_fit){
  #setting fitting tables
  model_fit_table = data.frame(matrix(NA, nrow= length(models_to_fit), ncol=3))
  colnames(model_fit_table) = c("model","llik","aicc")
  #setting estimate tables
  model_estimate_list = vector("list", length(models_to_fit))
  for (i in 1:length(models_to_fit)){
    # fitting models
    fit=OUwie(phy=tree, data=regimes, model=models_to_fit[i], lb=0, ub=Inf) 
    # picking fitting metrics
    model_fit_table[i,] = c(models_to_fit[i],fit$loglik,fit$AICc)
    # picking model estimates
    model_estimate = vector("list", 2)
    model_estimate[[1]] = fit$theta
    model_estimate[[2]] = fit$solution
    names(model_estimate) = c("theta","solution")
    model_estimate_list[[i]] = model_estimate
    names(model_estimate_list)[i] = models_to_fit[i]
  }
  fit_results = vector("list", 2)
  fit_results[[1]] = model_fit_table
  fit_results[[2]] = model_estimate_list
  names(fit_results) = c("fit_metrics","model_estimates")
  return(fit_results)
}
