### packages for PCM
library(ape)
library(phytools)
library(geiger)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

### other libraries
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

### create directory for pgls models
# check if dir exists
dir_check = dir.exists(paths="3_graphs")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "4_presentation_graphs", showWarnings = , recursive = FALSE, mode = "0777")
}


### phylogenetic tree location
trfn ="2_comparative_analyses/pruned_mcc_phylo.nwk"
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

#################################  altitude analysis ###########################

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

################### plotting phylogenetic reconstruction and traits #############

### loading relative probabilities at notes
relprobs_matrix = read.table("3_graphs/relprobs_matrix.csv", sep=",", h=T)

### setting states
# tip states probs
tip_states_probs = relprobs_matrix [1:n_tips, ]
# ancestral state probs
inner_node_probs = relprobs_matrix [(1+n_tips):(n_tips+n_inner_nodes),]
# state colors
state_cols=c( "#1E88E5","#D81B60", "#eb4683")
names(state_cols)=c("AF",  "other", "AFother")
# bar colrors
bar_colors = c()
for (i in 1:n_tips){
  bool = relprobs_matrix [i, ] == 1
  one_col =c("#1E88E5","#D81B60", "#eb4683")[bool]
  bar_colors = c(bar_colors, one_col)
}
names(bar_colors) = tr$tip.label

### trait vectors
trait_name = "anther_rel_diff"
trait= center_flower_df[[trait_name]]
names(trait) = center_flower_df$species

range(trait)

tiff("3_graphs/dec_mcc_ranges_3.tiff", units="in", width=3, height=6, res=600)
  plotTree(tree=tr,fsize=0.75, ftype="i")
  tiplabels(pie=tip_states_probs, piecol=state_cols, cex=0.75)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), pie= inner_node_probs, piecol=state_cols, cex=1.5)
  axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()

tiff(paste0("3_graphs/",trait_name,".tiff"), units="in", width=3.5, height=6, res=600)
  plotTree.wBars(tree=tr, x=trait,fsize=0.75, col=bar_colors, scale=2, ftype="i", method="plotTree")
  tiplabels(pie=tip_states_probs, piecol=state_cols, cex=0.5)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), pie= inner_node_probs, piecol=state_cols, cex=1)
  axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()



