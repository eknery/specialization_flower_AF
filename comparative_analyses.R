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
mean_pc_df = read.table("1_flower_analyses/mean_pc_df.csv", sep=",", h=T)
# sampled species
sampled_species = mean_pc_df$species
# flower data
flower_traits = mean_pc_df$pc1_score
names(flower_traits) = mean_pc_df$species

# geographic data
geo_states = mean_pc_df$state
names(geo_states) = mean_pc_df$species
  
### load hypervolume data
spp_hvolumes = read.table("2_hypervolume_inference/spp_hvolumes.csv",  sep=",", h=T)
# selecting only sampled species
samp_hvolumes = spp_hvolumes[spp_hvolumes$species %in% sampled_species,]
hvolumes = samp_hvolumes$hvolume
names(hvolumes) = samp_hvolumes$species

### load altitude data
spp_altitude = read.table("2_hypervolume_inference/spp_altitude.csv",  sep=",", h=T)
# selecting only sampled species
samp_altitude = spp_altitude[which(spp_altitude$species %in% sampled_species),]
altitude = samp_altitude$altitude
names(altitude) = samp_altitude$species

### counting pruned phylognetic trees
n_phylos = length(list.files("3_comparative_analyses/pruned_phylos"))

######################## niche divergence ~ flower traits ####################

predictor = flower_traits

plot(predictor, response)

n = 1

# load one phylogenetic tree
phylo_fn = paste("3_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(n), sep="")
phylo = read.tree(phylo_fn)

# load response variable
sister_hv_comparison = read.table(paste("2_hypervolume_inference/sister_hv_comparisons/sister_hv_comparison_", as.character(n), ".csv", sep=""), sep=',', h=T)
no = sister_hv_comparison$intersection/sister_hv_comparison$minimal_hv
names(no) = sister_hv_comparison$species
response = no


# fitting models to flower traits
fit_bm = fitContinuous(phy= phylo, response,  model="BM")
fit_ou = fitContinuous(phy= phylo, response,  model="OU")
# chosing best-fit model for response variable by aicc
if (fit_bm$opt$aicc - fit_ou$opt$aicc < 0 ){
  cor_mtx = corBrownian(value=1, phy= phylo, form=~1)
  print("best-fit BM")
} else {
  if (fit_bm$opt$aicc - fit_ou$opt$aicc < 2){
    cor_mtx = corBrownian(value=1, phy= phylo, form=~1)
    print("best-fit BM")
  } else {
    cor_mtx = corMartins(value= 1, phy = phylo, fixed = T)
    print("best-fit OU")
  }
}

# fitting pgls
fit_gls = gls(response ~ predictor, correlation=cor_mtx,  method = "REML")
# taking coefficients
intercept = fit_gls$coefficients[1]
angular = fit_gls$coefficients[-1]
# checking residuals
res = resid(fit_gls)[1:length(predictor)]
shapiro.test(res)
# calculating R model fit
ssr = sum(res^2)
sst =  sum((response- mean(response))^2)
r2 = 1 - (ssr/sst)





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
######### fitting evolutionary models to traits over trees #########################

### model fitting and selection functions
source("function_fit_evo_models.R")
source("function_choose_best.R")

### setting trait regimes
df = data.frame(mean_pc_df[,2], mean_pc_df[,1], mean_pc_df[,3])
colnames(df) = c("species", "state", "trait")
spp_trait_regimes = df

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
  tr$node.label = anc_node_states$state
  # fitting evolutionary models
  all_fits = fit_evo_models(tree=tr, regimes=spp_trait_regimes, models_to_fit = all_models)
  best_choice = choose_best(all_fits)
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
  print(i)
}

best_estimates =c()
for (i in 1:length(all_best_estimates)){
  best_estimates = rbind(best_estimates,all_best_estimates[[i]])
}

### exporting results
write.table(all_best_models,"3_comparative_analyses/best_models.csv",sep=",",quote=F,row.names=F)
write.table(best_estimates,"3_comparative_analyses/best_estimates.csv",sep=",",quote=F,row.names=F)

table(all_best_models$model)

best_estimates[all_best_models$model=="BMS",]

