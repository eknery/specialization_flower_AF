
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

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

## define geopraphic states 
high_ths = 0.90
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "AF"
geo_states[af_percentage <= low_ths] = "other"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "AFother"
names(geo_states) = spp_count_domain$species

### loading all flower traits
ftraits = read.table("0_data/flower_trait_matrix.csv", sep=",", h=T)
summary(ftraits)

## selfing-diagnostic traits
# flower size
flower_size = ftraits$hypathium_height
# herkogamy
minor_stamen_height = ftraits$minor_filet_height + ftraits$minor_anther_height
herkogamy = ftraits$style_height - minor_stamen_height
# pollen receipt & deposition
style_height = ftraits$style_height
stigma_size = ftraits$stigma_width
pore_size = ftraits$pore_long_section
# flower dataframe
flower_df = data.frame(flower_size, herkogamy, style_height, stigma_size, pore_size)

### sourcing other functions
source("function_pca_evaluation.R")

############################## flower differences by geographic state ###################

### Trait correlation analyses 
pca_evaluation(df= flower_df, iter=99, dir= paste(getwd(), "1_flower_analyses", sep="/") )
# observed patterns
pca = prcomp(flower_df, center = F)
stdev = pca$sdev / sum(pca$sdev)
load = pca$rotation
# export observed pca
write.table(stdev, paste(getwd(), "1_flower_analyses/observed_pca_stdev.csv", sep="/"), sep=",", quote=F, col.names=T)
write.table(load, paste(getwd(), "1_flower_analyses/observed_loadings.csv", sep="/"), sep=",", quote=F, col.names=T)

# sampling per group
table(samp_geo_states)

### dataframe with sp names
species = ftraits$species
flower_df = data.frame(species, flower_df)

# mean sp traits
mean_flower_df = aggregate(flower_df[,-1], by=list(flower_df$species), mean)

# sampled geographic states
sampled_bool = names(geo_states) %in% mean_flower_df$Group.1
samp_geo_states = geo_states[sampled_bool]

# geographic groups into flower data
mean_flower_df = data.frame(samp_geo_states, mean_flower_df)

aggregate(mean_flower_df[,-c(1,2)], by=list(mean_flower_df$samp_geo_states), mean)
