
### laoding libraries
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(PupillometryR)
library(RColorBrewer)
library(reshape2)

###
set.seed(42)

### loading occurrence count per domain
spp_count_domain = read.table("./0_data/spp_count_domain.csv", h=T, sep=",")

## define geopraphic states 
high_ths = 0.90
geo_state = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_state[af_percentage >= high_ths] = "AF"
geo_state[af_percentage < high_ths] = "other"
species = spp_count_domain$species
geo_state = data.frame(cbind(species, geo_state))

### loading all flower traits
ftraits = read.table("0_data/flower_trait_matrix.csv", sep=",", h=T)

# flower size
flower_size = ftraits$hypathium_height + ftraits$style_height

## pore size
pore_size = ftraits$pore_long_section

## relative pore size
rel_pore_size = ftraits$pore_long_section/ftraits$minor_anther_height

## stamen relative difference
major_stamen = ftraits$major_anther_height + ftraits$major_filet_height
minor_stamen = ftraits$minor_anther_height + ftraits$minor_filet_height
stamen_dim = (major_stamen - minor_stamen) / minor_stamen

## anther relative difference
anther_rel_diff = (ftraits$major_anther_height - ftraits$minor_anther_height )/ ftraits$minor_anther_height

## herkogamy
herkogamy = (ftraits$style_height) - (ftraits$hypathium_height + ftraits$minor_filet_height + ftraits$minor_anther_height)

# flower dataframe
species = ftraits$species
flower_df = data.frame(species,  
                       flower_size,
                       rel_pore_size, 
                       stamen_dim,  
                       herkogamy)

############################## flower traits by geographic state ###################

### merging geographic states
flower_df = flower_df %>% 
  merge(y= geo_state, by="species",all.x=TRUE)

### centering by geographic state and species
center_flower_df = flower_df |> 
  group_by(geo_state, species) |> 
  reframe(flower_size = median(flower_size),
          rel_pore_size = median(rel_pore_size),
          stamen_dim = median(stamen_dim), 
          herkogamy = median(herkogamy)
          )

summary_traits = center_flower_df  %>% 
  group_by(geo_state) |>
  reframe(median(flower_size),
          IQR(flower_size),
          median(rel_pore_size),
          IQR(rel_pore_size),
          median(stamen_dim), 
          IQR(stamen_dim),
          median(herkogamy),
          IQR(herkogamy)
  )


# exporting
write.table(center_flower_df, "1_flower_analyses/center_flower_df.csv", sep=",", quote=F, row.names=F, col.names=T)
write.table(summary_traits,"1_flower_analyses/summary_traits.csv", sep=",", quote=F, row.names=F, col.names=T)

######################## Testing floral trait differences  #######################

### loading species' flower traits
center_flower_df= read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)

### source function
source("scripts/function_permutation_test.R")

# set factor
factor = center_flower_df$geo_state

# loop over traits
all_traits = colnames(center_flower_df)[3:ncol(center_flower_df)]
all_pvalues = c()
for(j in 3:ncol(center_flower_df)){
  trait_name = colnames(center_flower_df)[j]
  response = as.data.frame(center_flower_df)[,j]
  pvalue = permutation_test(factor = factor, response = response, stats="median", iter=9999)
  all_pvalues = c(all_pvalues, pvalue)
}
all_pvalues = data.frame(all_traits, all_pvalues)
all_pvalues


write.table(all_pvalues,"1_flower_analyses/all_pvalues.csv", sep=",", quote=F, row.names=F, col.names=T)

############################## Assessing allometry ############################

### loading species' flower traits
center_flower_df= read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)

###### regarding whole plante
### height by species
center_height = aggregate(ftraits$plant_height, by=list(ftraits$species), function(x){median(x, na.rm = T) })

### drop non-height species
drop_it = which(is.na(center_height$x))
height = center_height$x[-drop_it]
flower = center_flower_df[-drop_it,]

### test relationship
# flower trait
trait_name = "rel_pore_size"
trait = flower[,trait_name]
# boundaries
bound = sd(trait)*1.96
max_bound = mean(trait) + bound
min_bound = mean(trait) - bound
# remove outliers
if(sum(trait > max_bound) != 0){
  drop= which(trait > max_bound)
  trait = trait[-drop]
  height = height[-drop]
}
if(sum(trait < min_bound) != 0){
  drop= which(trait < min_bound)
  trait = trait[-drop]
  height = height[-drop]
}
# fit model
model = lm(trait ~ height)
# test significance
test_model = summary(model)
r = round(test_model$r.squared, 2) # extract R square
pvalue = round(test_model$coefficients[,4][2],3)
# plot
tiff(paste("1_flower_analyses/allometry_whole_plant/",trait_name, ".tiff", sep=""), units="in", width=3.5, height=4, res=600)
  plot(x= height, y= trait, xlab= "height (m)", ylab= trait_name, col= "black")
  abline(model)
  text (x=8, y =0.15, labels = paste("R2: ",r))
  text (x=8, y =0.12, labels = paste("p-value: ",pvalue))
dev.off()

###### within flower
# flower trait
flower_size = center_flower_df$flower_size
trait_name = "rel_pore_size"
trait = center_flower_df[,trait_name]
# boundaries
bound = sd(trait)*1.96
max_bound = mean(trait) + bound
min_bound = mean(trait) - bound
# remove outliers
if(sum(trait > max_bound) != 0){
  drop= which(trait > max_bound)
  trait = trait[-drop]
  flower_size = flower_size[-drop]
}
if(sum(trait < min_bound) != 0){
  drop= which(trait < min_bound)
  trait = trait[-drop]
  flower_size = flower_size[-drop]
}
# fit model
model = lm(trait ~ flower_size)
# test significance
test_model = summary(model)
r = round(test_model$r.squared, 2) # extract R square
pvalue = round(test_model$coefficients[,4][2],3)
# plot
tiff(paste("1_flower_analyses/allometry_within_flower/",trait_name, ".tiff", sep=""), units="in", width=3.5, height=4, res=600)
  plot(x= flower_size, y= trait, xlab= "flower size (mm)", ylab= trait_name, col= "black")
  abline(model)
  text (x=8, y =0.15, labels = paste("R2: ",r))
  text (x=8, y =0.13, labels = paste("p-value: ",pvalue))
dev.off()


############################### flower trait structure ########################## 

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