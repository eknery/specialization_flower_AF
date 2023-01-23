
### no slope
sister_no_metrics = read.table("3_hypervolume_inference/sister_no_metrics.csv", sep=',', h=T)
# no vector
angular_no = sister_no_metrics$angular_no
af_no_slope = angular_no[sister_no_metrics$state == "AF"]
non_no_slope = angular_no[sister_no_metrics$state == "other"]
no = c(af_no_slope, non_no_slope)


### listing traits 
all_trait_names = list.files("2_comparative_analyses/OUWIE")


### best estimates fro floral traits
best_estimates = read.table("2_comparative_analyses/OUWIE/pore_size/best_estimates.csv", sep=",", h=T)

# take trait estimates
theta = c(best_estimates$theta_1, best_estimates$theta_2)
# reliable estimates
lw_bound = median(theta) - IQR(theta)/2
up_bound = median(theta) + IQR(theta)/2
index = which(theta > lw_bound & theta < up_bound)
# select reliable theta
theta_clean = theta[index]
# select corresponding no values
no_slope = no[index]

# retrive common numbers 
linear = lm(no_slope ~ theta_clean)
summary(linear)
plot(theta_clean, no_slope)
abline(linear)

plot_param = ggplot(data= param_df, aes(x=state, y=parameter, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1, alpha = 0.65) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.25)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab(param_name)+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "non-endemic")) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL), panel.border=element_rect(fill=NA,colour="black"), axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=6), legend.position = "none") 
