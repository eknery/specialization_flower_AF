
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
geo_state = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_state[af_percentage >= high_ths] = "AF"
geo_state[af_percentage <= low_ths] = "other"
geo_state[af_percentage > low_ths & af_percentage < high_ths] = "other"
names(geo_state) = spp_count_domain$species

### loading all flower traits
ftraits = read.table("0_data/flower_trait_matrix.csv", sep=",", h=T)
## selfing-diagnostic traits
# flower size
flower_size = ftraits$hypathium_height + ftraits$style_height
# herkogamy
stamen_height = ftraits$minor_filet_height + ftraits$minor_anther_height
herkogamy = ftraits$style_height - stamen_height
# pollen presentation & receipt
pore_size = ftraits$pore_long_section
stigma_size = ftraits$stigma_width
# flower dataframe
flower_df = data.frame(flower_size, herkogamy, pore_size, stigma_size)

################################ flower trait structure #################### 

### sourcing other functions
source("function_pca_evaluation.R")
pca_pvalues = pca_evaluation(df= flower_df, iter=999, dir= paste(getwd(), "1_flower_analyses", sep="/") )

# picking significant
pca_significant = which(pca_pvalues[[1]] < 0.05)

### observed patterns
pca = prcomp(flower_df, center = F)
stdev = pca$sdev / sum(pca$sdev)
load = pca$rotation
# export observed pca
write.table(stdev, paste(getwd(), "1_flower_analyses/observed_pca_stdev.csv", sep="/"), sep=",", quote=F, col.names=T)
write.table(load, paste(getwd(), "1_flower_analyses/observed_loadings.csv", sep="/"), sep=",", quote=F, col.names=T)

### mean sp score
# pc scores
pc_df = data.frame(ftraits$species, pca$x[,pca_significant])
colnames(pc_df)[1] = "species"
# mean sp traits
pc_flower_df = aggregate(pc_df[,-1], by=list(pc_df$species), median)
# sampled geographic states
sampled_bool = names(geo_state) %in% pc_flower_df$Group.1
samp_geo_state = geo_state[sampled_bool]
# geographic groups into pc data
pc_flower_df = data.frame(samp_geo_state, pc_flower_df)
colnames(pc_flower_df)[1:2] = c("state","species")
write.table(pc_flower_df, paste(getwd(), "1_flower_analyses/pc_flower_df.csv", sep="/"), sep=",", quote=F, row.names=F, col.names=T)

### plotting
# my colors
mycols = c( "#1E88E5", "#D81B60")
names(mycols) = c("AF","other")
# plot
tiff("1_flower_analyses/pc_by_geography.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= pc_flower_df, aes(x=state, y=PC1, fill=state)) +
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

### center sp trait
species = ftraits$species
flower_df = data.frame(species, flower_df)
# center sp traits
center_flower_df = aggregate(flower_df[,-1], by=list(flower_df$species), median)

### traits by geographic group
# sampled geographic states
sampled_bool = names(geo_state) %in% center_flower_df$Group.1
samp_geo_state = geo_state[sampled_bool]
# geographic groups into flower data
center_flower_df = data.frame(samp_geo_state, center_flower_df)
# descriptive statistics
geo_center = aggregate(center_flower_df[,-c(1,2)], by=list(center_flower_df$samp_geo_state), mean)
geo_disper = aggregate(center_flower_df[,-c(1,2)], by=list(center_flower_df$samp_geo_state), sd)
# name
colnames(center_flower_df)[1:2] = c("state", "species")
# exporting
write.table(center_flower_df, "1_flower_analyses/center_flower_df.csv", sep=",", quote=F, row.names=F, col.names=T)
write.table(geo_center, "1_flower_analyses/trait_center_per_geography.csv", sep=",", quote=F, row.names=F, col.names=T)
write.table(geo_disper,"1_flower_analyses/trait_dispersion_per_geography.csv", sep=",", quote=F, row.names=F, col.names=T)

######################## plotting traits by geography #######################

### loading species' flower traits
center_flower_df= read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)

### plotting
# my colors
mycols = c( "#1E88E5","#D81B60")
# plot list
plot_list = list()
# loop over variables
for(i in 3:ncol(center_flower_df) ){
  trait_df = center_flower_df[,c(1,2,i)]
  trait_name =  colnames(trait_df)[3]
  colnames(trait_df)[3] = "trait"
  # checking outliers
  med = median(trait_df$trait)
  bound= IQR(trait_df$trait)*1.5
  up_bound = med + bound
  out_sum = sum(trait_df$trait > up_bound)
  if (out_sum != 0){
    trait_df = trait_df[trait_df$trait < up_bound,]
  }
  # plot
  one_plot = ggplot(data= trait_df, aes(x=state, y=trait, fill=state)) +
    geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
    geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
    scale_fill_manual(values=mycols)+
    scale_colour_manual(values=mycols)+
    xlab("geographic distribution")+ ylab(trait_name)+
    scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "non-endemic"))+
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

