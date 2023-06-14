# requires hypervolume

hypervolume_tss= function(data, absences, k, thresholds=0.5){
  # folds
  folds=kfold(data, k=k)
  # preparing performance table
  performance = data.frame(matrix(0, nrow= 1, ncol=5))
  colnames(performance) = c('threshold', 'accuracy','sensitivity', 'specificity', 'tss')
  for (i in 1:k){
    test = data[folds == i,]
    train = data[folds != i,]
    ### estimating species bandwidth and training hypervolume
    for (ths in thresholds){
        tryCatch({
        data_hv = hypervolume_gaussian(train,  samples.per.point = ceiling((10^(3 + sqrt(ncol(train))))/nrow(train)), 
                                     kde.bandwidth = estimate_bandwidth(train), sd.count = 3, chunk.size = 100,
                                     quantile.requested = ths, quantile.requested.type = "probability")
      
        ### predicting true presences
        data_test=hypervolume_inclusion_test(data_hv, test, reduction.factor = 1, fast.or.accurate ="fast")
        a=sum(data_test, na.rm = TRUE)
        c=length(data_test)-sum(data_test, na.rm = TRUE)
        ### predicting virtual absences
        abs_test = hypervolume_inclusion_test(data_hv, absences, reduction.factor = 1, fast.or.accurate ="fast")
        b=sum(abs_test, na.rm = TRUE)
        d=length(abs_test)-sum(abs_test, na.rm = TRUE)
        ### performance statistics
        n=nrow(test)+nrow(absences)
        accuracy= (a + d)/n
        sensitivity= a/(a+c)
        specificity= d/(b+d)
        tss= sensitivity + (specificity -1)
        performance = rbind(performance, c(ths, accuracy, sensitivity, specificity, tss) )
        },error=function(e){} )
        }
    }
  performance = performance[-1,]
  return(performance)
}
