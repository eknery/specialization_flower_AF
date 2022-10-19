
pca_evaluation = function(df, iter=99, dir = getwd() ){
  test_results = list()
  ### observed values 
  pca = prcomp(df, center = F)
  stdev = pca$sdev / sum(pca$sdev)
  load = pca$rotation
  pc_numbers = paste("PC", as.character(1:length(pca$sdev)), sep='')
  names(stdev) = pc_numbers
  ### random values
  # set result objects
  rand_stdev_df = data.frame(matrix(0, nrow=1, ncol=length(stdev) ) )
  rand_load_list = vector('list', iter)
  #loop!
  for (n in 1:iter){
    rand_num = sample(x=1:nrow(df), size = nrow(df), replace = F)
    rand_data = df[rand_num, 1]
    for (j in 2:ncol(df)){
      rand_num = sample(x=1:nrow(df), size = nrow(df), replace = F)
      rand_data = cbind(rand_data, df[rand_num, j])
    }
    rand_pca = prcomp(rand_data, center = F)
    rand_stdev = rand_pca$sdev/ sum(rand_pca$sdev)
    rand_load_list[[n]] = rand_pca$rotation
    rand_stdev_df = rbind(rand_stdev_df, rand_stdev)
  }
  rand_stdev_df = rand_stdev_df[-1,]
  ### analyzing stdev values
  pvalues = c()
  for(i in 1:length(stdev) ){
    obs_stdev = stdev[i]
    rand_distibution = rand_stdev_df[,i]
    greater = sum(obs_stdev > rand_distibution)
    total =  length(rand_distibution) +1 
    p = 1 - (greater / total)
    pvalues = c(pvalues, p)
  }
  names(pvalues) = pc_numbers
  test_results[[1]] = pvalues
  # exporting test values
  file_name = "stdev_pca_test.csv"
  write.table(pvalues, paste(dir,file_name, sep="/"), sep=",", quote=F, col.names=T)
  # summary stdev
  rand_min = apply(rand_stdev_df, MARGIN = 2, min)
  rand_median = apply(rand_stdev_df, MARGIN = 2, median)
  rand_max = apply(rand_stdev_df, MARGIN = 2, max)
  stdev_summary = data.frame(pc_numbers, rand_min, rand_median, rand_max, stdev) 
  # ploting
  one_plot = ggplot(data= stdev_summary, aes(x=pc_numbers, y=stdev ) ) +
    geom_errorbar(aes(ymin=rand_min, ymax=rand_max), col="gray", size=1, width=0.01)+
    geom_point(aes(y=rand_median), color="gray", size = 1.5) +
    geom_point(aes(y=stdev), color="red", size = 0.75) +
    xlab("principal component")+ ylab("standard deviation")+
    theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=12))
  # exporting plot
  file_name = "stdev_pca.tiff"
  tiff(paste(dir,file_name, sep="/"), units="in", width=3.5, height=3, res=600)
    print(one_plot)
  dev.off()
  ### selecting loadings from significant PCs
  # select significant pcs
  sign_pcs = which(pvalues < 0.05)
  # setting result objects
  sign_load_list = vector('list', length(sign_pcs))
  names(sign_load_list) = names(sign_pcs)
  # loop!
  for (j in sign_pcs){
    position = which(j == sign_pcs)
    rand_pc_loads = rand_load_list[[1]][,j]
    for (n in 2:iter){
      rand_pc_loads = cbind(rand_pc_loads, rand_load_list[[n]][,j])
    }
    sign_load_list[[position]] = rand_pc_loads
  }
  ### analysing loading values
  plot_list = list()
  for (j in sign_pcs){
    # indexing
    position = as.numeric(which(j == sign_pcs))
    # take null distribution by pc axis
    load_distribution = data.frame(sign_load_list[[position]])
    # observed loadings
    obs_load = load[,j] # !!!
    load_df = cbind(obs_load, load_distribution)
    colnames(load_df)[1] = "observed_load"
    # summary load
    variables = rownames(load)
    rand_min = apply(load_df[,-1], MARGIN = 1, min)
    rand_median = apply(load_df[,-1], MARGIN = 1, median)
    rand_max = apply(load_df[,-1], MARGIN = 1, max)
    obs_load = load_df[,1]
    load_summary = data.frame(variables, rand_min, rand_median, rand_max, obs_load) 
    # plotting
    one_plot = ggplot(data= load_summary, aes(x=variables, y=obs_load ) ) +
      geom_errorbar(aes(ymin=rand_min, ymax=rand_max), col="gray", size=1, width=0.05)+
      geom_point(aes(y=rand_median), color="gray", size = 1) +
      geom_point(aes(y=obs_load), color="red", size = 0.75) +
      xlab("variable")+ ylab("loadings")+
      theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=12), axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1))
    plot_list[[position]] = one_plot
    # testing
    pvalues = c()
    for (i in 1:nrow(load)){
      load_median = median( as.numeric(load_df[i,-1]) )
      if (load_df[i,1] >  load_median){
        greater = sum(load_df[i,1] > load_df[i,-1])
        total =  length(load_df[i,-1]) + 1 
        p = 1 - (greater / total)
        names(p) = rownames(load)[i]
        pvalues = c(pvalues, p)
      } else {
        smaller = sum(load_df[i,1] <= load_df[i,-1])
        total =  length(load_df[i,-1]) + 1 
        p = 1 - (smaller / total)
        names(p) = rownames(load)[i]
        pvalues = c(pvalues, p)
      }
    }
    file_name = paste("pc_loading_test_", as.character(j),".csv", sep="")
    write.table(pvalues, paste(dir,file_name, sep="/"), sep=",", quote=F, col.names=T)
  }
  
  for (j in sign_pcs){
    # indexing
    position = as.numeric(which(j == sign_pcs))
    file_name = paste("pc_loading_", as.character(j),".tiff", sep="")
    tiff(paste(dir,file_name, sep="/"), units="in", width=3.5, height=3, res=600)
    print(plot_list[[position]])
    dev.off()
  }
  return(test_results)
}