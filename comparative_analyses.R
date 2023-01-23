### require packages
library(ape)
library(phytools)
library(geiger)
library(OUwie)
library(nlme)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

### load flower pc scores
center_flower_df = read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)
# choose flower proxy!
flower_proxy = center_flower_df
# sampled species
sampled_species = center_flower_df$species
# geographic state
state = center_flower_df$state
names(state) = center_flower_df$species

### load mcc phylogeentic tree
mcc_phylo = read.tree("2_comparative_analyses/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("2_comparative_analyses/pruned_phylos"))

### plotting libraries
library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)
library(ggplot2)
# my colors
mycols = c( "#1E88E5", "#D81B60")
names(mycols) = c("AF", "other")

############################## biogeographic reconstruction ###################

### create directory for DEC results
# check if dir exists
dir_check = dir.exists(paths="2_comparative_analyses/DEC_ancestral_reconstructions")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "2_comparative_analyses/DEC_ancestral_reconstructions", showWarnings = , recursive = FALSE, mode = "0777")
}

# reading range data
geog_fn = ("0_data/spp_distribution_af.data")
moref(geog_fn)

# converting phylip format to tipranges
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geog_fn)
tipranges

# setting maximum number of areas occupied for reconstructions
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))

# Initialize DEC model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# location of the geography text file
BioGeoBEARS_run_object$geogfn = geog_fn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$min_branchlength = 0.001    

# set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$include_null_range = FALSE    

# computing options
BioGeoBEARS_run_object$num_cores_to_use = 1

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

#### DEC over trees 
for (i in 1:n_phylo){ 
  # phylogeny tree location
  trfn = paste("2_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(i), sep="")
  tr = read.tree(trfn)
  # n tips and nodes
  n_tips = Ntip(tr)
  n_nodes = tr$Nnode
  # inputting tree into DEC
  BioGeoBEARS_run_object$trfn = trfn
  # fitting DEC
  res_DEC = bears_optim_run(BioGeoBEARS_run_object)
  # node marginal ML
  relprobs_matrix = res_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
  # node states
  state_labels=c("AF", "other", "AFother")
  node_states = get_ML_states_from_relprobs(relprobs=relprobs_matrix, statenames=state_labels, returnwhat = "states", if_ties = "takefirst")
  # getting only ancestral nodes
  anc_node = (n_tips+1):(n_tips+n_nodes)
  state = node_states[anc_node]
  anc_node_states = data.frame(anc_node, state)
  write.table(anc_node_states , paste("2_comparative_analyses/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_"), sep=",", row.names=F, quote=F)
}

################################# current altitude analysis ###########################

### load packages
library(raster)
library(sp)
library(sf)

### loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

### load altitude raster
ras_alt = raster("0_data/rasters/altitude.gri")

### altitude values per species
# extarct values
alt_values = raster::extract(ras_alt ,spp_points[,2:3])
# center to species
center_alt_values = aggregate(alt_values, by=list(spp_points$species), function(x){median(x, na.rm=T)} )
# keep only species with sampled flowers
spp_altitude = center_alt_values[center_alt_values$Group.1 %in% sampled_species,]
# name columns
colnames(spp_altitude) = c("species", "altitude")
# add states
spp_altitude = data.frame(state, spp_altitude)
# export
write.table(spp_altitude, "2_comparative_analyses/spp_altitude.csv", sep=",", row.names= F, quote = F)

### testing altitude difference
# difference per distribution
center_sp_alt = aggregate(spp_altitude$altitude, by= list(spp_altitude$state), median)
dispersion_sp_alt = aggregate(spp_altitude$altitude, by= list(spp_altitude$state), IQR)

# source permutation test
source("function_permutation_test.R")
permutation_test(factor= spp_altitude$state, response= spp_altitude$altitude,iter=999, out_dir = "2_comparative_analyses", name= "altitude.tiff")

### plotting
# text size
axis_title_size = 10
x_text_size = 8
# plot
tiff("2_comparative_analyses/altitude_per_distribution.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_altitude, aes(x=state, y=altitude, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  ylim(c(0,1500))+
  xlab("geographic distribution")+ ylab("species' elevation (m)")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "non-endemic"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),axis.text.y = element_text(angle = 90),legend.position = "none")
dev.off()

####################### ancestral reconstruction of altitude ###################

### load species' altitude
spp_altitude = read.table("2_comparative_analyses/spp_altitude.csv", sep=",", h=T)
# altitude vector
altitude = spp_altitude$altitude
names(altitude) = spp_altitude$species

### acestral reconstruction over trees
anc_data = c()
for (i in 1:n_phylo){ 
  # phylogeny tree location
  trfn = paste("2_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(i), sep="")
  tr = read.tree(trfn)
  # n tips and nodes
  n_tips = Ntip(tr)
  n_nodes = tr$Nnode
  # node ages
  node_ages = round(node.depth.edgelength(tr),5)
  present = round(max(node_ages), 5)
  node_ages = node_ages - present
  # ancestral node ages
  anc_node_ages = node_ages[(n_tips+1):(n_tips+n_nodes)]
  # ancestral biogeographic state
  dec_fn = paste("2_comparative_analyses/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_")
  anc_node_states = read.table(dec_fn, sep=",", h=T)
  anc_node_states$state[anc_node_states$state == "AFother"] = "other"
  # ancestral reconstruction of altitude
  anc_altitude = as.numeric( fastAnc(tree = tr, x = altitude) )
  # organizing into one dataframe
  one_rec = data.frame(anc_node_states, anc_node_ages, anc_altitude)
  # update!
  anc_data = rbind(anc_data, one_rec)
}

### dividing into time intervals
# time boundaries
old_age = min(anc_data$anc_node_ages)
new_age = round(max(anc_data$anc_node_ages), 3)
# intervals
intervals = anc_data$anc_node_ages
breaks = seq(new_age, old_age, by= (old_age - new_age)/15)
for (i in 1:length(breaks)){
  intervals[which(anc_data$anc_node_ages > breaks[i+1] & anc_data$anc_node_ages < breaks[i])] = (breaks[i] + breaks[i+1])/2
}
anc_data = data.frame(anc_data, intervals)

### summarizing reconstruction by state across intervals
list_anc_data = split(anc_data, f= anc_data$state)
all_summary_anc = data.frame()
for (i in 1:length(list_anc_data)){
  central = aggregate(list_anc_data[[i]]$anc_altitude, by= list(list_anc_data[[i]]$intervals),  function(x){median(x, na.rm=T)})
  dispersion  = aggregate(list_anc_data[[i]]$anc_altitude, by= list(list_anc_data[[i]]$intervals), function(x){IQR(x, na.rm=T)} )
  state = rep(names(list_anc_data)[i], nrow(central) )
  summary_anc = cbind(state, central, dispersion[,-1])
  all_summary_anc = rbind(all_summary_anc, summary_anc)
}
colnames(all_summary_anc) = c("state", "age", "central", "dispersion")

# drop negative ancestral values
all_summary_anc$central - (all_summary_anc$dispersion/2)

### plotting
# text size
axis_title_size = 10
x_text_size = 8

tiff("2_comparative_analyses/altitude_reconstruction.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= all_summary_anc, aes(x=age, y=central, group= state, color=state) ) +
  geom_point(size = 1, alpha = 1) +
  geom_line(size=1)+
  geom_errorbar(size=0.75, width=0, aes(ymin=central-(dispersion/2), ymax=central+(dispersion/2)))+
  scale_colour_manual(values=mycols)+
  xlim(c(-12,-0.5))+ ylim(c(0,1500))+
  xlab("time before present (m.y.a.)")+ ylab("species' elevation (m)")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),axis.text.y = element_text(angle = 90),legend.position = "none")
dev.off()

################### fitting evolutionary models to traits over trees ################

# check if dir exists
dir_check = dir.exists(paths="2_comparative_analyses/OUWIE")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "2_comparative_analyses/OUWIE", showWarnings = , recursive = FALSE, mode = "0777")
}

### model fitting and selection functions
source("function_fit_evo_models.R")
source("function_choose_best.R")

### setting species, regimes, and trait dataframe
species = flower_proxy$species 
regime = flower_proxy$state
trait_df = flower_proxy[,-c(1:2)]

# setting max parameter numbers
max_param_num = length(unique(regime))*4

### loop over all traits
for (j in 1:ncol(trait_df)){ # 
  ### choose trait
  trait = trait_df[,j] ## j !
  trait_name = colnames(trait_df)[j]
  spp_trait_regimes = data.frame(species, regime, trait)
  ### setting output dir
  # check if output dir exists
  dir_check = dir.exists(paths=paste("2_comparative_analyses/OUWIE/",trait_name, sep="") )
  # create output dir if not created yet
  if (dir_check == FALSE){
    dir.create(path= paste("2_comparative_analyses/OUWIE/",trait_name, sep=""), showWarnings = , recursive = FALSE, mode = "0777")
  }
  # output dir
  out_dir = paste("2_comparative_analyses/OUWIE/",trait_name, sep="")
  ### setting result objects
  all_best_models = data.frame(matrix(NA, nrow= n_phylo, ncol=4))
  colnames(all_best_models) = c(c("model","llik","aicc","delta_aicc"))
  all_best_estimates = vector("list" , n_phylo)
  all_models = c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
  ### fitting all models to all phylogenetic trees
  for (i in 1:n_phylo){ 
    # phylogenetic tree
    tr_fn = paste("2_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(i), sep="")
    tr = read.tree(tr_fn)
    # DEC ancestral states
    dec_fn = paste("2_comparative_analyses/DEC_ancestral_reconstructions/anc_node_states_", as.character(i), sep="")
    anc_node_states = read.table(dec_fn, sep=",", h=T)
    # replacing 'AFother' per 'other'
    anc_node_states$state[anc_node_states$state == "AFother"] = "other"
    # setting ancestral states into phylogenetic tree
    tr$node.label = anc_node_states$state
    # fitting evolutionary models
    all_fits = fit_evo_models(tree=tr, regimes=spp_trait_regimes, models_to_fit = all_models)
    best_choice = choose_best(all_fits)
    all_best_models[i,] = best_choice$best_fit
    all_best_estimates[[i]] = best_choice$best_estimates
    print(i)
  }
  # organizing estimates into dataframe
  best_estimates =c()
  for (i in 1:length(all_best_estimates)){
    param_num = length(all_best_estimates[[i]])
    if ( param_num < max_param_num ){
      diff = max_param_num - param_num
      na_param = rep(NA, diff)
      all_best_estimates[[i]] = c(all_best_estimates[[i]], na_param)
    } 
    best_estimates = rbind(best_estimates,all_best_estimates[[i]])
  }
  ### exporting results
  write.table(all_best_models, paste(out_dir,"/best_models.csv", sep=""), sep=",",quote=F,row.names=F)
  write.table(best_estimates, paste(out_dir,"/best_estimates.csv", sep=""), sep=",",quote=F,row.names=F)
}

############################## describing best-fit model #########################

### loading libraries
library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)

### listing traits 
all_trait_names = list.files("2_comparative_analyses/OUWIE")

# my colors
mycols = c( "#1E88E5", "#D81B60")

### loop over all traits
for (trait_name in all_trait_names){
  ### setting output dir
  dir = paste("2_comparative_analyses/OUWIE/",trait_name, sep="")
  ### load model fit and estimates
  all_best_models = read.table(paste(dir,"/best_models.csv", sep=""), sep=",", h=T)
  best_estimates = read.table(paste(dir,"/best_estimates.csv", sep=""), sep=",", h=T)
  ### find model with best-fit across trees
  best_fit_count = table(all_best_models$model)
  max_count = max(best_fit_count)
  best_model = best_fit_count[best_fit_count == max_count]
  best_model_name = names(best_model)
  # export
  write.table(best_fit_count, paste(dir,"/best_model_count.csv", sep=""), sep=",", quote=F, row.names=F) 
  ### retrieve best estimates
  model_estimates = best_estimates[all_best_models$model == best_model_name,]
  n_rows = nrow(model_estimates)
  # take sigma
  sigma_sq_df = model_estimates[,which(colnames(model_estimates) %in% c("sigma_sq_1","sigma_sq_2", "sigma_sq_3"))]
  sigma_sq = c()
  if( is.vector(sigma_sq_df) ){
    sigma_sq = sigma_sq_df
  } else {
    for (j in 1:ncol(sigma_sq_df)){ sigma_sq = c(sigma_sq, sigma_sq_df[,j]) }
  }
  sigma_sq = as.numeric(sigma_sq)
  sigma = sqrt(sigma_sq)
  # take theta
  theta_df = model_estimates[,which(colnames(model_estimates) %in% c("theta_1","theta_2", "theta_3"))]
  theta = c()
  if( is.vector(theta_df) ){
    theta = theta_df
  } else {
    for (j in 1:ncol(theta_df)){ theta = c(theta, theta_df[,j]) }
  }  
  # take alpha
  alpha_df = model_estimates[,which(colnames(model_estimates) %in% c("alpha_1","alpha_2", "alpha_3"))]
  alpha = c()
  if( is.vector(alpha_df) ){
    alpha = alpha_df
  } else {
    for (j in 1:ncol(alpha_df)){ alpha = c(alpha, alpha_df[,j]) }
  }
  # state
  state = c(rep("AF", n_rows), rep("other", n_rows))
  # organize into dataframe
  estimates_df = data.frame(state, sigma, theta, alpha)
  ### describing parameter estimates
  mean_estimates = aggregate(estimates_df[,-1], by=list(estimates_df[,1]), function(x){mean(x, na.rm=T)})
  sd_estimates = aggregate(estimates_df[,-1], by=list(estimates_df[,1]), function(x){sd(x, na.rm=T)})
  # export
  write.table(mean_estimates, paste(dir,"/mean_estimates.csv", sep=""), sep=",", row.names=F, quote=F)
  write.table(sd_estimates, paste(dir,"/sd_estimates.csv", sep=""), sep=",", row.names=F, quote=F)
  ### plotting parameter estimates
  # looping over parameters
  for (param_name in c("sigma", "theta", "alpha")){
    # picking one parameter
    param_index = which(param_name == colnames(estimates_df)) - 1 
    param_df = data.frame(estimates_df[,"state"], estimates_df[,param_name])
    colnames(param_df) = c("state", "parameter")
    # checking outliers
    med = median(param_df$parameter, na.rm = T) 
    bound= IQR(param_df$parameter, na.rm = T)*1.5
    up_bound = med + bound
    out_sum = sum(param_df$parameter > up_bound, na.rm = T)
    if (out_sum != 0){  param_df = param_df[param_df$parameter < up_bound,] }
    # plot parameter
    plot_param = ggplot(data= param_df, aes(x=state, y=parameter, fill=state)) +
      geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1, alpha = 0.65) +
      geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.25)+
      scale_fill_manual(values=mycols)+
      scale_colour_manual(values=mycols)+
      xlab("geographic distribution")+ ylab(param_name)+
      scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "non-endemic")) +
      theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL), panel.border=element_rect(fill=NA,colour="black"), axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=6), legend.position = "none") 
   # export plot
    tiff(paste(dir, "/",param_name, ".tiff", sep=""), units="in", width=2.5, height=2, res=600)
      print(plot_param)
    dev.off()
  }
}

