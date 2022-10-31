### packages
library(raster)
library(sp)
library(sf)
library(spatialEco)
library(rgeos)
library(hypervolume)
library(dismo)

### load flower pc scores
center_flower_df = read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)
# sampled species
sampled_species = center_flower_df$species
# geographic data
state = center_flower_df$state
names(state) = center_flower_df$species

### loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

### loading bg coordinates
bg_points=read.table("0_data/bg_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

### loading raster layers
ras1 = raster("0_data/rasters/temperature_diurnal_range")
ras2 = raster("0_data/rasters/precipitation_seasonality")
ras3 = raster("0_data/rasters/solar_radiation")
ras4 = raster("0_data/rasters/soil_pH.gri")
env_ras= stack(ras1,ras2,ras3,ras4)

### counting pruned phylognetic trees
n_phylos = length(list.files("2_comparative_analyses/pruned_phylos"))

################################## raster preparation ###########################

### croping rasters
#set extent
ext=extent(-100, -33, -34, 24)
allras_crop=crop(env_ras[[1]], ext)
for (i in 2:length(env_ras@layers) ){
  ras_crop =  crop(env_ras[[i]], ext)
  allras_crop= raster::stack(allras_crop, ras_crop)
}
crop_env_ras = allras_crop

### keeping mean and sd values
ras_means= cellStats(crop_env_ras, stat='mean', na.rm=TRUE)
ras_sds = cellStats(crop_env_ras, stat='sd', na.rm=TRUE)

### converting rasters to z-scores
scale_env_ras = crop_env_ras[[c(1:nlayers(env_ras))]]
for (i in 1:nlayers(scale_env_ras)) { 
  scale_env_ras[[i]] = (scale_env_ras[[i]] - cellStats(scale_env_ras[[i]], 'mean')) / cellStats(scale_env_ras[[i]], 'sd') 
}

### treating spp occurrences
# extracting spp z-values
scale_spp_env = raster::extract(scale_env_ras,spp_points[,2:3])
scale_spp_env = data.frame(spp_points$species, scale_spp_env)
colnames(scale_spp_env)[1] = "species"

# removing spp NAs
spp_nas = c()
for (i in 2:length(scale_spp_env[1,])){
  spp_nas= c(spp_nas, which(is.na(scale_spp_env[,i])) )}
if (length(spp_nas) > 0){
  scale_spp_env = scale_spp_env[-spp_nas,]}

### treating background
# extracting bg z-values
scale_bg_env = raster::extract(scale_env_ras,bg_points)

# removing bg NAs
bg_nas = c()
for (i in 1:length(scale_bg_env[1,])){
  bg_nas= c(bg_nas, which(is.na(scale_bg_env[,i])) )}
if (length(bg_nas) > 0){
  scale_bg_env = scale_bg_env[-bg_nas,]}

########################### evaluating hypervolume prediction ########################

source("funtion_hypervolume_tss.R")

# formatting background to spatial points
bg_sp = SpatialPoints(bg_points)
crs(bg_sp) = '+proj=longlat +datum=WGS84 +no_defs'  

# species names
all_spp_names = sort(unique(spp_points$species))

# list of performance tables
all_perf_tables = vector('list', length(all_spp_names))
names(all_perf_tables) = all_spp_names

# thresholds to test
my_ths = c(0.05, 0.25, 0.5, 0.75, 0.95)

# ssp undergoing leave-one-out validation
narrow_spp = c("angelana","capixaba" ,"dura","kollmannii","kriegeriana", "penduliflora","suberosa")

### looping over species
for (sp_name in all_spp_names){
  index = which(all_spp_names == sp_name)
  sp_data = scale_spp_env[scale_spp_env$species == sp_name, -1]
  # removing sp_points from background
  sp_point = spp_points[spp_points$species == sp_name,2:3]
  circ_around <- circles(sp_point, d=7000, lonlat=TRUE)
  poly_around <- polygons(circ_around)
  abs_points = erase.point(bg_sp, poly_around, inside = TRUE)
  # extracting z-scores for background
  one_scale_abs = raster::extract(scale_env_ras, abs_points)
  # removing absence NAs
  one_abs_nas = c()
  for (i in 1:length(one_scale_abs[1,])){
    one_abs_nas = c(one_abs_nas, which(is.na(one_scale_abs[,i])) )}
  if (length(one_abs_nas) > 0){
    one_scale_abs = one_scale_abs[-one_abs_nas,]}
  # evaluating models
  if (sp_name %in% narrow_spp){
    all_perf_tables[[index]] = hypervolume_tss(data=sp_data, absences=one_scale_abs, k=nrow(sp_data), thresholds = my_ths)
  } else {
    all_perf_tables[[index]] = hypervolume_tss(data=sp_data, absences=one_scale_abs, k= 3, thresholds = my_ths)
  }
}

# organizing into single table
hv_performance = data.frame(matrix(0,nrow=1, ncol=6,))
colnames(hv_performance) = c("species","threshold","accuracy","sensitivity","specificity","tss")

for (i in 1:length(all_perf_tables)){
  sp_name = names(all_perf_tables)[i]
  species = rep(sp_name, nrow(all_perf_tables[[i]]) )
  one_performance = data.frame(species, all_perf_tables[[i]])
  hv_performance = rbind(hv_performance, one_performance)
}

write.table(hv_performance, "2_hypervolume_inference/hv_performance.csv", sep=',', quote=F, row.names=F)

# mean performance by threshold
mean_hv_performance = data.frame(matrix(0,nrow=1, ncol=6,))
colnames(mean_hv_performance) = c("species","Group.1","accuracy","sensitivity","specificity","tss")

for(sp_name in all_spp_names){
  sp_performance = hv_performance[hv_performance$species==sp_name,]
  mean_sp_performance = aggregate(sp_performance[3:6], by=list(sp_performance$threshold), mean)
  species = rep(sp_name, nrow(mean_sp_performance))
  mean_sp_performance = data.frame(species, mean_sp_performance)
  mean_hv_performance = rbind(mean_hv_performance, mean_sp_performance)
}

mean_hv_performance = mean_hv_performance[-1,]
colnames(mean_hv_performance)[2] = "threshold"
write.table(mean_hv_performance, "2_hypervolume_inference/mean_hv_performance.csv", sep=',', quote=F, row.names=F)

############################ fitting hypervolumes ######################

# loading threshold values
mean_hv_performance = read.table("2_hypervolume_inference/mean_hv_performance.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# getting best thresholds for each species
best_ths = rep(NA, length(all_spp_names))
for(sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  mean_sp_performance = mean_hv_performance[mean_hv_performance$species==sp_name,]
  sp_th = mean_sp_performance$threshold[mean_sp_performance$tss == max(mean_sp_performance$tss)] 
  best_ths[index] = sp_th[1]
}
names(best_ths) = all_spp_names

# fitting models per species
all_spp_hv = list()
for (sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  sp_data = scale_spp_env[scale_spp_env$species== sp_name, 2:5]
  sp_ths = best_ths[names(best_ths) == sp_name]
  sp_band = estimate_bandwidth(sp_data)
  sp_hv = hypervolume_gaussian(sp_data,  samples.per.point = ceiling((10^(3 + sqrt(ncol(sp_data))))/nrow(sp_data)), 
                               kde.bandwidth = sp_band, sd.count = 3, chunk.size = 100,
                               quantile.requested = sp_ths, quantile.requested.type = "probability")
  all_spp_hv[[index]] = sp_hv
}
names(all_spp_hv) = all_spp_names

############################# sister hv comparison ############################

# sourcing other functions
source("function_sister_pairs.R")

### setting sister taxa for each tree
for (n in 1:n_phylos){
  # phylogeny tree location
  phylo_fn = paste("2_comparative_analyses/pruned_phylos/pruned_phylo_", as.character(n), sep="")
  phylo = read.tree(phylo_fn)
  # phylogenetic distance
  phylo_distance = cophenetic(phylo)
  # sister taxa
  sister_taxa_list = sister_pairs(phylo_distance)
  ### sister divergence time
  sister_divergence = c()
  for(focal_sp in sampled_species){
    focal_sp_dists= round(phylo_distance[focal_sp,],5)
    min_phylo_dist = round( min(phylo_distance[focal_sp,][phylo_distance[focal_sp,] != 0]), 5)
    sister_divergence = c(sister_divergence, min_phylo_dist)
  }
  ### sister hv comparison
  sister_hv_comparison = matrix(0, nrow=length(sampled_species), ncol=4)
  for (sp_name in sampled_species){
    index = which(sp_name == sampled_species)
    # sp hv
    sp_hv = all_spp_hv[[sp_name]]
    # sister hv
    sister_name = sister_taxa_list[[sp_name]]
    if (length(sister_name) > 1){
      one_sister_name = sister_name[1]
      sister_hv = all_spp_hv[[one_sister_name]]
      for (i in 2:length(sister_name) ){
        one_sister_name = sister_name[i]
        one_sister_hv = all_spp_hv[[one_sister_name]]
        sister_set = hypervolume_set(sister_hv, one_sister_hv, check.memory = F)
        sister_hv = sister_set[[4]]
      }
    } else{
      sister_hv = all_spp_hv[[sister_name]]
    }
    # getting sp and sister volumes
    sister_hv_comparison[index,1] = sp_hv@Volume
    sister_hv_comparison[index,2] = sister_hv@Volume
    # getting intersection and union volumes 
    hv_set = hypervolume_set(sp_hv, sister_hv, check.memory = F)
    sister_hv_comparison[index,3] = hv_set[[3]]@Volume # taking intersection
    sister_hv_comparison[index,4] = hv_set[[4]]@Volume # taking union
  }
  # minimal hypervolume per sister pair
  min_hv = apply(sister_hv_comparison[,1:2], MARGIN = 1, FUN= min)
  # dataframe
  sister_hv_comparison = data.frame(sampled_species, sister_hv_comparison, min_hv, sister_divergence)
  colnames(sister_hv_comparison) =c("species", "sp_hv", "sister_hv", "intersection", "union", "minimal_hv","divergence_time")
  #exporting
  write.table(sister_hv_comparison, paste("2_hypervolume_inference/sister_hv_comparisons/sister_hv_comparison_", as.character(n), ".csv", sep="" ), sep=',', quote=F, row.names=F)
  # update me!
  print(paste("Time:", Sys.time(), "Loop iterarion:", as.character(n) ) )
}

######################## calculating sister NO metrics #########################

sister_no_metrics = c()
for (n in 1:n_phylos ){ 
  sister_hv_comparison = read.table(paste("2_hypervolume_inference/sister_hv_comparisons/sister_hv_comparison_", as.character(n), ".csv", sep=""), sep=',', h=T)
  geo_groups = split(sister_hv_comparison, state)
  for (i in 1:length(geo_groups) ){
    group_name = names(geo_groups)[i]
    one_group = geo_groups[[i]]
    no = one_group$intersection/one_group$union
    mean_no = mean(no)
    linear_model = lm(no ~ one_group$divergence_time)
    intercept_no = linear_model$coefficients["(Intercept)"]
    angular_no = linear_model$coefficients[2]
    one_group_metric = c(group_name, mean_no, intercept_no, angular_no)
    sister_no_metrics = rbind(sister_no_metrics, one_group_metric)
  }
}

sister_no_metrics = data.frame(sister_no_metrics)
colnames(sister_no_metrics) = c("state", "mean_no", "intercept_no", "angular_no")
rownames(sister_no_metrics) = NULL
sister_no_metrics$mean_no = as.numeric(sister_no_metrics$mean_no)
sister_no_metrics$intercept_no = as.numeric(sister_no_metrics$intercept_no)
sister_no_metrics$angular_no = as.numeric(sister_no_metrics$angular_no)

str(sister_no_metrics)

#exporting
write.table(sister_no_metrics, "2_hypervolume_inference/sister_no_metrics.csv", sep=',', quote=F, row.names=F)

############################## analyzing RO metrics ############################

sister_no_metrics = read.table("2_hypervolume_inference/sister_no_metrics.csv", sep=',', h=T)

### summarizing metrics 
means = aggregate(sister_no_metrics$angular_no, by = list(sister_no_metrics$state), mean )
sds = aggregate(sister_no_metrics$angular_no, by = list(sister_no_metrics$state), sd )
summary_no_group = cbind(means, sds[,2])
colnames(summary_no_group) = c("state", "mean", "sd")
# export
write.table(summary_no_group, "2_hypervolume_inference/summary_no_group.csv", sep=',', row.names=F, quote=F)

### permutation test
# set comparison
i = 1
j = 2
diff = means$x[i] -  means$x[j]
# random comparisons
iterations = 1000
rand_diff =  c()
for(n in 1:iterations){
  rand_state = sample(sister_no_metrics$state)
  rand_means = aggregate(sister_no_metrics$angular_no, by = list(rand_state), mean )
  rand_diff[n] = rand_means$x[i] - rand_means$x[j]
}

if (diff < median(rand_diff)){
  pvalue = 1 - ( sum(diff < rand_diff) / (iterations +1) )
} else{
  pvalue = 1 - ( sum(diff >= rand_diff) / (iterations +1) )
}
print(paste("p-value:", pvalue))

################################## plotting ##################################
### packages
library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)
library(ggplot2)

### graphical parameters
# my colors
mycols = c( "#1E88E5", "#D81B60")
names(mycols) = c("AF", "other")
# text size
axis_title_size = 10
x_text_size = 8

tiff("2_hypervolume_inference/angular_no_per_distribution.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= sister_no_metrics, aes(x=state, y=angular_no, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("angular NO")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "non-endemic"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),axis.text.y = element_text(angle = 90),legend.position = "none")
dev.off()
