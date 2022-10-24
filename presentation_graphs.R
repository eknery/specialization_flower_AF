### require packages
library(ape)
library(phytools)
library(geiger)
library(OUwie)
library(nlme)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

### create directory for pgls models
# check if dir exists
dir_check = dir.exists(paths="4_presentation_graphs")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "4_presentation_graphs", showWarnings = , recursive = FALSE, mode = "0777")
}

### phylogenetic tree location
trfn ="3_comparative_analyses/pruned_mcc_phylo.nwk"
tr = read.tree(trfn)

# n tips and nodes
n_tips = Ntip(tr)
n_inner_nodes = tr$Nnode

### load flower pc scores
center_flower_df = read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)

# sampled species
sampled_species = center_flower_df$species

# geographic state
state = center_flower_df$state

### load species' occupied elevation
spp_altitude = read.table("2_hypervolume_inference/spp_altitude.csv", sep=",", h=T)
# sampled indexes
index_sampled = which(spp_altitude$species %in% sampled_species)
altitude = spp_altitude$altitude[index_sampled ]
names(altitude) = spp_altitude$species[index_sampled]


############################### fitting DEC #################################
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
# inputting tree into DEC
BioGeoBEARS_run_object$trfn = trfn
# fitting DEC
res_DEC = bears_optim_run(BioGeoBEARS_run_object)
# node marginal ML
relprobs_matrix = res_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
colnames(relprobs_matrix) =c("AF", "other", "AFother")

write.table(relprobs_matrix,"4_presentation_graphs/relprobs_matrix.csv", sep=",", quote=F, row.names=F)

################### plotting phylogenetic reconstruction and tratis #############

### loading relative probabilities at notes
relprobs_matrix = read.table("4_presentation_graphs/relprobs_matrix.csv", sep=",", h=T)

### setting states
# tip states probs
tip_states_probs = relprobs_matrix [1:n_tips, ]
# ancestral state probs
inner_node_probs = relprobs_matrix [(1+n_tips):(n_tips+n_inner_nodes),]
# state colors
state_cols=c( "#1E88E5","#D81B60", "#eb4683")
names(state_cols)=c("AF",  "other", "AFother")

tiff("4_presentation_graphs/dec_mcc_ranges.tiff", units="in", width=4, height=6, res=600)
  plotTree(tree=tr,fsize=0.75, ftype="i")
  tiplabels(pie=tip_states_probs, piecol=state_cols, cex=0.4)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), pie= inner_node_probs, piecol=state_cols, cex=0.8)
  axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()

########################## plotting spatial features ##########################

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

### my colors
mycols = c( "#1E88E5", "#D81B60")

### creating a altitude dataframe
alt_df = data.frame(state, altitude)

# plot
tiff("4_presentation_graphs/altitude_per_geography.tiff", units="in", width=2, height=2, res=600)
ggplot(data= alt_df, aes(x=state, y=altitude, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1, alpha = 0.65) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.25, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("altitude (m)")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "non-endemic")) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL), panel.border=element_rect(fill=NA,colour="black"), axis.title=element_text(size=8,face="bold"), axis.text=element_text(size=6), legend.position = "none") 
dev.off()
