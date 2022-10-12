
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
flower_size = ftraits$hypathium_height + ftraits$style_height
# herkogamy
minor_stamen_height = ftraits$minor_filet_height + ftraits$minor_anther_height
herkogamy = ftraits$style_height - minor_stamen_height
# pollen presentation & receipt
pore_size = ftraits$pore_long_section
stigma_size = ftraits$stigma_width
# flower dataframe
flower_df = data.frame(flower_size, herkogamy, pore_size, stigma_size)


################################ flower trait structure #################### 

### sourcing other functions
source("function_pca_evaluation.R")
pca_evaluation(df= flower_df, iter=99, dir= paste(getwd(), "2_flower_analyses", sep="/") )

### observed patterns
pca = prcomp(flower_df, center = F)
stdev = pca$sdev / sum(pca$sdev)
load = pca$rotation
# export observed pca
write.table(stdev, paste(getwd(), "1_flower_analyses/observed_pca_stdev.csv", sep="/"), sep=",", quote=F, col.names=T)
write.table(load, paste(getwd(), "1_flower_analyses/observed_loadings.csv", sep="/"), sep=",", quote=F, col.names=T)

### mean sp score
# pc scores
pc_df = data.frame(flower_df$species, pca$x[,1])
colnames(pc_df) = c("species", "pc1_score")
# mean sp traits
flower_pc_df = aggregate(pc_df[,-1], by=list(pc_df$species), mean)
# sampled geographic states
sampled_bool = names(geo_states) %in% flower_pc_df$Group.1
samp_geo_states = geo_states[sampled_bool]
# geographic groups into pc data
flower_pc_df = data.frame(samp_geo_states, flower_pc_df)
colnames(flower_pc_df) = c("state","species", "pc1_score")
write.table(flower_pc_df, paste(getwd(), "1_flower_analyses/flower_pc_df.csv", sep="/"), sep=",", quote=F, row.names=F, col.names=T)

tiff("1_flower_analyses/pca_by_geography.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= flower_pc_df, aes(x=state, y=pc1_score, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("PC1 score")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()

############################## flower traits by geographic state ###################

### mean sp trait
species = ftraits$species
flower_df = data.frame(species, flower_df)
# mean sp traits
mean_flower_df = aggregate(flower_df[,-1], by=list(flower_df$species), mean)

### traits by geographic group
# sampled geographic states
sampled_bool = names(geo_states) %in% mean_flower_df$Group.1
samp_geo_states = geo_states[sampled_bool]
# geographic groups into flower data
mean_flower_df = data.frame(samp_geo_states, mean_flower_df)
# descriptive statistics
geo_center = aggregate(mean_flower_df[,-c(1,2)], by=list(mean_flower_df$samp_geo_states), median)
geo_disper = aggregate(mean_flower_df[,-c(1,2)], by=list(mean_flower_df$samp_geo_states), IQR)
# export observed pca
write.table(geo_center, paste(getwd(), "1_flower_analyses/trait_center_per_geography.csv", sep="/"), sep=",", quote=F, row.names=F, col.names=T)
write.table(geo_disper, paste(getwd(), "1_flower_analyses/trait_dispersion_per_geography.csv", sep="/"), sep=",", quote=F, row.names=F, col.names=T)

### plotting traits by geography
# my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")
# plot list
plot_list = list()
# loop over variables
for(i in 3:ncol(mean_flower_df) ){
  trait_df = mean_flower_df[,c(1,2,i)]
  trait_name =  colnames(trait_df)[3]
  colnames(trait_df)[3] = "trait"
  one_plot = ggplot(data= trait_df, aes(x=samp_geo_states, y=trait, fill=samp_geo_states)) +
    geom_point(aes(color=samp_geo_states),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
    geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
    scale_fill_manual(values=mycols)+
    scale_colour_manual(values=mycols)+
    xlab("geographic distribution")+ ylab(trait_name)+
    scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
    theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
  plot_list[[i-2]] = one_plot
  names(plot_list)[i-2] = trait_name
}

for (i in 1:length(plot_list) ){
  trait_name = names(plot_list)[i] 
  file_name = paste(trait_name,".tiff", sep="")
  tiff(paste("1_flower_analyses",file_name, sep="/"), units="in", width=3.5, height=3, res=600)
    print(plot_list[[i]])
  dev.off()
}
