### plotting libraries
library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)
library(ggplot2)

### create directory for regressions
# check if dir exists
dir_check = dir.exists(paths="3_hypervolume_inference/regression_analyses")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "3_hypervolume_inference/regression_analyses", showWarnings = , recursive = FALSE, mode = "0777")
}

### graphical params
mycols = c( "#1E88E5", "#D81B60")
names(mycols) = c("AF", "other")

### no slope
sister_no_metrics = read.table("3_hypervolume_inference/sister_no_metrics.csv", sep=',', h=T)
# no_slope vector
angular_no = sister_no_metrics$angular_no
af_no_slope = angular_no[sister_no_metrics$state == "AF"]
non_no_slope = angular_no[sister_no_metrics$state == "other"]
no_slope = c(af_no_slope, non_no_slope)

### listing traits 
all_trait_names = list.files("2_comparative_analyses/OUWIE")

### regression parameters
reg_params = c()

### looping over traits
for (trait_name in all_trait_names){
  ### setting dir and axis name  
  dir_name = paste(trait_name,"/best_estimates.csv", sep="")
  axis_name =  paste("theta", trait_name, sep=" ")
  ### estimates of floral trait evolution
  best_estimates = read.table(paste("2_comparative_analyses/OUWIE/",dir_name, sep=""), sep=",", h=T)
  # take trait estimates
  theta = c(best_estimates$theta_1, best_estimates$theta_2)
  # states
  state = c( rep("AF", length(best_estimates$theta_1)) , rep("other", length(best_estimates$theta_2) ) )
  #### reliable estimate
  # bounds
  lw_bound = median(theta) - IQR(theta)*1.5
  up_bound = median(theta) + IQR(theta)*1.5
  index = which(theta > lw_bound & theta < up_bound)
  # select reliable theta
  theta_clean = theta[index]
  # select corresponding no values
  no_slope_clean = no_slope[index]
  # select corresponding state values
  state_clean = state[index]
  ### linear regression
  # fit line
  linear = lm(no_slope_clean ~ theta_clean)
  reg_intercept = linear$coefficients["(Intercept)"]
  reg_slope = linear$coefficients["theta_clean"]
  # test
  summary_linear = summary(linear)
  r2 = summary_linear$adj.r.squared
  pvalue = round(summary_linear$coefficients[2,4], 3)
  # take all reg param
  one_reg = c(trait_name, reg_intercept, reg_slope, r2, pvalue)
  ### update regression params
  reg_params = rbind(reg_params, one_reg)
  ### plot linear relationship
  # theta dataframe
  theta_df = data.frame(state_clean, theta_clean, no_slope_clean)
  # plotting
  plot_linear = ggplot(data= theta_df, aes(x=theta_clean, y=no_slope_clean)) +
    geom_point(aes(color=state_clean),position = position_jitter(width = 0.07), size = 1, alpha = 0.65) +
    geom_abline(intercept = reg_intercept, slope = reg_slope) +
    scale_fill_manual(values=mycols)+
    scale_colour_manual(values=mycols)+
    xlab(axis_name) + ylab("NO slope")+
    theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL), panel.border=element_rect(fill=NA,colour="black"), axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=6), legend.position = "none") 
  # export
  tiff(paste("3_hypervolume_inference/regression_analyses/",trait_name, ".tiff", sep=""), units="in", width=3.5, height=3, res=600)
    print(plot_linear)
  dev.off()
}

write.table(reg_params, "3_hypervolume_inference/regression_analyses/reg_params.csv", sep=",", quote = F, row.names = F)
