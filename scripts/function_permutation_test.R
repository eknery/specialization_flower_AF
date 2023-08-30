
permutation_test = function(factor, response, stats = "mean", iter=99){
  if (stats == "mean"){
    obs_means = aggregate(response, by=list(factor), FUN=mean)
    obs_diff = obs_means$x[1] - obs_means$x[2]
    # loop over random factor distribution
    rand_diff_distribution = c()
    for (i in 1:iter){
      rand_factor = sample(factor)
      rand_means = aggregate(response, by=list(rand_factor), FUN=mean)
      rand_diff = rand_means$x[1] - rand_means$x[2]
      rand_diff_distribution = c(rand_diff_distribution, rand_diff)
    }
  }
  if (stats == "median"){
    obs_medians = aggregate(response, by=list(factor), FUN=median)
    obs_diff = obs_medians$x[1] - obs_medians$x[2]
    # loop over random factor distribution
    rand_diff_distribution = c()
    for (i in 1:iter){
      rand_factor = sample(factor)
      rand_medians = aggregate(response, by=list(rand_factor), FUN=median)
      rand_diff = rand_medians$x[1] - rand_medians$x[2]
      rand_diff_distribution = c(rand_diff_distribution, rand_diff)
    }
  }
  # test tail
  if(obs_diff < median(rand_diff_distribution) ){
    out = sum(obs_diff < rand_diff_distribution)
  } 
  if(obs_diff > median(rand_diff_distribution) ){
    out = sum(obs_diff > rand_diff_distribution)
  }
  # p-value
  pvalue = 1 - (out/ (iter + 1))
  names(pvalue) ="p_value"
  return(pvalue)
  
}



