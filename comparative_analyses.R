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
geo_state = center_flower_df$state
names(geo_state) = center_flower_df$species

### load species' occupied elevation
spp_altitude = read.table("2_hypervolume_inference/spp_altitude.csv", sep=",", h=T)
# sampled indexes
index_sampled = which(spp_altitude$species %in% sampled_species)
altitude = spp_altitude$altitude[index_sampled ]
names(altitude) = spp_altitude$species[index_sampled]

### load mcc phylogeentic tree
mcc_phylo = read.tree("3_comparative_analyses/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("3_comparative_analyses/pruned_phylos"))

############################### fitting pgls models ###########################

### create directory for pgls models
# check if dir exists
dir_check = dir.exists(paths="3_comparative_analyses/PGLS")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "3_comparative_analyses/PGLS", showWarnings = , recursive = FALSE, mode = "0777")
}

### setting common variables
# choose phylogenetic tree
phylo = mcc_phylo

# set predictor variable
predictor = geo_state

# set dataframe with response variables
response_df = data.frame(flower_proxy[,-c(1,2)])
naming_vector = flower_proxy$species

### on the relationship
#flower_size = linear
#herkogamy = linear
#pore_size = log
#sitgma_size = log

log_these = c("altitude", "pore_size", "sitgma_size")

for (j in 1:ncol(response_df) ){
  ### picking one response variable
  response_name = colnames(response_df)[j]
  ### setting output dir
  # check if output dir exists
  dir_check = dir.exists(paths=paste("3_comparative_analyses/PGLS/", response_name, sep="") )
  # create output dir if not created yet
  if (dir_check == FALSE){
    dir.create(path= paste("3_comparative_analyses/PGLS/", response_name, sep=""), showWarnings = , recursive = FALSE, mode = "0777")
  }
  # ouput dir
  out_dir = paste("3_comparative_analyses/PGLS/", response_name, sep="")
  ### setting result objects
  best_evo_model = c("model", "aicc")
  intercept = c()
  angular = c()
  t_value = c()
  p_value = c()
  shapiro_w = c()
  shapiro_p = c()
  r2 = c()
  ### applying leave-one-out
  for (i in 1:nrow(response_df)){
    # droping n observatio in the predictor
    iter_predictor = predictor[-i]
    # selecting response variable
    iter_response = response_df[-i,j]
    iter_names = naming_vector[-i]
    names(iter_response) = iter_names
    # need transformation?
    if (response_name %in% log_these){iter_response = log(iter_response)}
    # initial values 
    init_bm = sd(iter_response)
    init_ou = mean(iter_response) - min(iter_response)
    # drop species from phylogenetic tree
    to_drop = naming_vector[i]
    iter_phylo = drop.tip(phylo, to_drop)
    # fitting models to flower traits
    fit_bm = fitContinuous(phy= iter_phylo, iter_response,  model="BM")
    fit_ou = fitContinuous(phy= iter_phylo, iter_response,  model="OU")
    # chosing best-fit model for response variable by aicc
    if (fit_bm$opt$aicc < fit_ou$opt$aicc){
      cor_mtx = corBrownian(value=init_bm, phy= iter_phylo, form=~iter_names)
      best_evo_model = rbind(best_evo_model, c("BM", fit_bm$opt$aicc))
    } else {
      if (fit_bm$opt$aicc - fit_ou$opt$aicc < 2){
        cor_mtx = corBrownian(value=init_bm, phy= iter_phylo, form=~iter_names)
        best_evo_model = rbind(best_evo_model, c("BM", fit_bm$opt$aicc))
      } else {
        cor_mtx = corMartins(value= init_ou, phy = iter_phylo, fixed = T, form=~iter_names)
        best_evo_model = rbind(best_evo_model, c("OU", fit_ou$opt$aicc))
      }
    }
    # fitting pgls
    fit_gls = gls(iter_response ~ iter_predictor, correlation=cor_mtx,  method = "REML")
    summary_gls = summary(fit_gls)
    # taking coefficients and test
    gls_test_table = summary_gls$tTable
    intercept = c(intercept, summary_gls$coefficients[1])
    angular = rbind(angular, summary_gls$coefficients[-1])
    t_value = c(t_value, gls_test_table[,3][2])
    p_value = rbind(p_value, gls_test_table[,4])
    # checking residuals
    res = resid(fit_gls)[1:length(iter_response)]
    shapiro = shapiro.test(res)
    shapiro_w = c(shapiro_w, shapiro$statistic)
    shapiro_p = c(shapiro_p, shapiro$p.value)
    # calculating R model fit
    ssr = sum(res^2)
    sst =  sum((iter_response- mean(iter_response))^2)
    r2 = c(r2, c(1 - (ssr/sst)) )
    print(paste("round:", as.character(i)))
    }
  print(paste("response done:", response_name))
  # organizing model evalution and pgls fit
  pgls_stats = data.frame(intercept, angular, t_value, p_value, r2, shapiro_w, shapiro_p)
  # exporting
  write.table(best_evo_model, paste(out_dir,"/best_evo_model.csv", sep=""), row.names=F, quote=F, sep=",")
  write.table(pgls_stats, paste(out_dir,"/pgls_stats.csv", sep=""), row.names=F, quote=F, sep=",")
}

############################## biogeographic reconstruction ###################

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
  trfn = paste("3_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(i), sep="")
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
  write.table(anc_node_states , paste("3_comparative_analyses/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_"), sep=",", row.names=F, quote=F)
}

################### fitting evolutionary models to traits over trees ################

### create directory for ouwie models
dir.create(path= "3_comparative_analyses/OUWIE", showWarnings = , recursive = FALSE, mode = "0777")

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
for (j in 1:ncol(trait_df)){
  ### choose trait
  trait = trait_df[,j]
  trait_name = colnames(trait_df)[j]
  spp_trait_regimes = data.frame(species, regime, trait)
  ### setting output dir
  # check if output dir exists
  dir_check = dir.exists(paths=paste("3_comparative_analyses/OUWIE/",trait_name, sep="") )
  # create output dir if not created yet
  if (dir_check == FALSE){
    dir.create(path= paste("3_comparative_analyses/OUWIE/",trait_name, sep=""), showWarnings = , recursive = FALSE, mode = "0777")
  }
  # output dir
  out_dir = paste("3_comparative_analyses/OUWIE/",trait_name, sep="")
  ### setting result objects
  all_best_models = data.frame(matrix(NA, nrow= n_phylo, ncol=4))
  colnames(all_best_models) = c(c("model","llik","aicc","delta_aicc"))
  all_best_estimates = vector("list" , n_phylo)
  all_models = c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
  ### fitting all models to all phylogenetic trees
  for (i in 1:n_phylo){ 
    # phylogenetic tree
    tr_fn = paste("3_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(i), sep="")
    tr = read.tree(tr_fn)
    # DEC ancestral states
    dec_fn = paste("3_comparative_analyses/DEC_ancestral_reconstructions/anc_node_states_", as.character(i), sep="")
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

### laoding libraries
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
trait_names = list.files("3_comparative_analyses/OUWIE")

# my colors
mycols = c( "#1E88E5", "#D81B60")

### loop over all traits
for (trait_name in trait_names){
  ### setting output dir
  dir = paste("3_comparative_analyses/OUWIE/",trait_name, sep="")
  ### load model fit and estimates
  all_best_models = read.table(paste(dir,"/best_models.csv", sep=""), sep=",", h=T)
  best_estimates = read.table(paste(dir,"/best_estimates.csv", sep=""), sep=",", h=T)
  ### find model with best-fit across trees
  best_fit_count = table(all_best_models$model)
  max_count = max(best_fit_count)
  best_model = best_fit_count[best_fit_count == max_count]
  best_model_name = names(best_model)
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

