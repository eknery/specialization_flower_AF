
### packages
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(phytools)
library(geiger)
library(OUwie)
library(nlme)

### load pc scores
mean_pc_df = read.table("1_flower_analyses/mean_pc_df.csv", sep=",", h=T)
# sampled species
sampled_species = mean_pc_df$species

### load phylogenetic trees
# mcc tree
mcc = read.tree("0_data/mcc_phylo.nwk")
# random sample of trees

