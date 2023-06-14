
permutation_test = function(factor, response, stats = "mean", iter=99, out_dir=getwd(), name=".tiff"){
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
    # test tail
    if(obs_diff < median(rand_diff_distribution)){
      out = sum(obs_diff < rand_diff_distribution)
    } else {
      out = sum(obs_diff > rand_diff_distribution)
    }
    # p-value
    pvalue = 1 - (out/ (iter + 1))
    names(pvalue) ="p_value"
    # graph
    tiff(paste(out_dir, name, sep="/"), units="in", width=3.5, height=3, res=600)
    hist(rand_diff_distribution, main="")
    abline(v=obs_diff, col="red", lty=2, lwd=2)
    dev.off()
    # otuput
    return(pvalue)
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
    # test tail
    if(obs_diff < median(rand_diff_distribution)){
      out = sum(obs_diff < rand_diff_distribution)
    } else {
      out = sum(obs_diff > rand_diff_distribution)
    }
    # p-value
    pvalue = 1 - (out/ (iter + 1))
    names(pvalue) ="p_value"
    # graph
    tiff(paste(out_dir, name, sep="/"), units="in", width=3.5, height=3, res=600)
    hist(rand_diff_distribution, main="")
    abline(v=obs_diff, col="red", lty=2, lwd=2)
    dev.off()
    # otuput
    return(pvalue)
  }
}



