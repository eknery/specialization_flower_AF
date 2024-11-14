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
  dir.create(path= "3_graphs")
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
#species and geo states
spp_state = center_flower_df %>% 
  select(species, geo_state)

tr_order = c("cyathanthera","cipoensis","cubatanensis","capixaba", "discolor",
             "budlejoides", "racemifera", "fasciculata", "pepericarpa", 
             "willdenowii", "flammea", "formosa", "cinerascens", "sclerophylla",
             "corallina", "shepherdii", "petroniana", "castaneiflora", "lymanii",
             "hyemalis", "multispicata", "valtheri", "ruficalyx", "lepidota",
             "polyandra", "macrothyrsa", "heliotropoides", "dispar", "ferruginata",
             "burchellii", "tiliifolia", "rufescens", "punctata", "dichrophylla",
             "chrysophylla", "cowanii", "amnicola", "longispicata", "fallax",
             "stenostachya", "pterocaulon", "macuxi", "mayarae", "argyrophylla",
             "alborufescens", "serialis", "secundiflora", "navioensis", "hypoleuca",
             "lourteigiana", "albicans")

#################################  altitude analysis ###########################

### read altitude values
spp_alt = read.table("0_data/spp_alt_values.csv", sep=",", h=T)

### summarize by species
med_alt = spp_alt %>% 
  group_by(spp_points.species) %>% 
  reframe(alt = median(spp_alt_values),
          q1 = quantile(spp_alt_values, prob = 0.25),
          q3 = quantile(spp_alt_values, prob = 0.75)
          ) %>% 
  dplyr::rename(species = spp_points.species) %>% 
  mutate(species = factor(species, levels = rev(tr_order) ))

### 
med_alt = med_alt %>% 
  plyr::join(y = spp_state, type = "left", by = "species") %>% 
  filter(!is.na(geo_state))

med_alt %>% 
  group_by(geo_state) %>% 
  reframe(median(alt), IQR(alt))

### ploting
alt_plot = ggplot(data= med_alt, 
                  aes(x=species, 
                      y= alt)
                  ) +
  
  geom_point(aes(color=geo_state),
             size = 1, 
             alpha = 0.65,
             position = position_jitter(width = 0.07)
             ) +
  
  geom_linerange(aes(ymin = q1, 
                     ymax = q3,
                     color=geo_state)
                 )+
  ylim(0, 1600)+
  
  scale_color_manual(values = c("#1E88E5","#D81B60") )+
  
  coord_flip() +
  
  xlab("") +
  ylab("elevation (m)") +
  
  theme(panel.background=element_rect(fill="white"),
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=6, face="bold"),
        axis.text.x = element_text(size= 4, angle = 45),
        axis.text.y = element_text(size= 5, angle = 0),
        legend.position = "none"
        )

# export plot
tiff("3_graphs/alt_plot.tiff",
     units="cm", width=3, height=9.50, res=600)
  print(alt_plot)
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

write.table(relprobs_matrix,"3_graphs/relprobs_matrix.csv", sep=",", quote=F, row.names=F)

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
  plotTree(tree=tr,
           type="fan",
           fsize=0.75, 
           ftype="i")
  tiplabels(pie=tip_states_probs, 
            piecol=state_cols, 
            cex=0.75)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), 
             pie= inner_node_probs, 
             piecol=state_cols, 
             cex=1.5)
  axisPhylo(pos=c(0.5), 
            font=3, 
            cex.axis=0.5)
dev.off()

tiff(paste0("3_graphs/",trait_name,".tiff"), units="in", width=3.5, height=6, res=600)
  plotTree.wBars(tree=tr, 
                 x=trait,
                 fsize=0.75, 
                 col=bar_colors, 
                 scale=2, 
                 ftype="i", 
                 method="plotTree")
  tiplabels(pie=tip_states_probs, 
            piecol=state_cols, 
            cex=0.5)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes),
             pie= inner_node_probs, 
             piecol=state_cols,
             cex=1)
  axisPhylo(pos=c(0.5), 
            font=3, 
            cex.axis=0.5)
dev.off()



