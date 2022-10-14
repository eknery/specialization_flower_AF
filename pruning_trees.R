### require packages
library(ape)
library(phytools)
library(geiger)

### load flower pc scores
center_flower_df = read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)
# sampled species
sampled_species = center_flower_df$species


### load mcc phylogenetic tree
mcc_phylo = read.tree("0_data/mcc_phylo.nwk")

# prunning mcc phylogenetic tree
pruned_mcc_phylo = drop.tip(mcc_phylo, mcc_phylo$tip.label[-match(sampled_species, mcc_phylo$tip.label)])
# exporting prunned mcc
write.tree(pruned_mcc_phylo, file = "3_comparative_analyses/pruned_mcc_phylo.nwk", append = FALSE, digits = 10, tree.names = FALSE)


### loading random sample of trees
rand_phylos = read.tree("0_data/100_rand_phylos.nwk")

# pruning trees to sampled species
pruned_phylos = rand_phylos
for (i in 1:length(rand_phylos)){
  pruned_phylos[[i]] = drop.tip(rand_phylos[[i]], rand_phylos[[i]]$tip.label[-match(sampled_species, rand_phylos[[i]]$tip.label)])
}
#exporting phylogenetic trees
for (i in 1:length(pruned_phylos)){
  write.tree(pruned_phylos[[i]], file = paste("3_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(i), sep=""), append = FALSE, digits = 10, tree.names = FALSE)
}