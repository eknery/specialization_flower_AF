setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

### packages
library(raster)
library(sp)
library(sf)
library(spatialEco)
library(rgeos)
library(hypervolume)
library(dismo)

#loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# loading bg coordinates
bg_points=read.table("0_data/bg_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# loading raster layers
ras1 = raster("0_data/rasters/temperature_diurnal_range")
ras2 = raster("0_data/rasters/precipitation_seasonality")
ras3 = raster("0_data/rasters/solar_radiation")
ras4 = raster("0_data/rasters/soil_pH.gri")
env_ras= stack(ras1,ras2,ras3,ras4)


################################## data preparation ###########################

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
scale_spp_env = extract(scale_env_ras,spp_points[,2:3])
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

##################### evaluating hypervolume prediction ####################

source("3_hypervolume_inference/funtion_hypervolume_tss.R")

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

write.table(hv_performance, "3_hypervolume_inference/hv_performance.csv", sep=',', quote=F, row.names=F)

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
write.table(mean_hv_performance, "3_hypervolume_inference/mean_hv_performance.csv", sep=',', quote=F, row.names=F)

########################### inferring environmental niches #############################

# loading threshold values
mean_hv_performance=read.table("3_hypervolume_inference/mean_hv_performance.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# getting best thresholds for each species
best_ths = rep(NA, length(all_spp_names))
for(sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  mean_sp_performance = mean_hv_performance[mean_hv_performance$species==sp_name,]
  sp_th = mean_sp_performance$threshold[mean_sp_performance$tss == max(mean_sp_performance$tss)] 
  best_ths[index] = sp_th[1]
}
names(best_ths) = all_spp_names

# getting spp niches
hv_spp_list = vector('list', length(all_spp_names))
hv_scale_positions = data.frame(matrix(NA, nrow=length(all_spp_names), ncol=ncol(scale_spp_env) ) )
spp_hvolumes = rep(NA, length(all_spp_names))

for (sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  sp_data = scale_spp_env[scale_spp_env$species== sp_name, 2:5]
  sp_ths = best_ths[index]
  sp_band = estimate_bandwidth(sp_data)
  sp_hv = hypervolume_gaussian(sp_data,  samples.per.point = ceiling((10^(3 + sqrt(ncol(sp_data))))/nrow(sp_data)), 
                               kde.bandwidth = sp_band, sd.count = 3, chunk.size = 100,
                               quantile.requested = sp_ths, quantile.requested.type = "probability")
  hv_spp_list[[index]] = sp_hv
  spp_hvolumes[index] = sp_hv@Volume
  max_prob = which(sp_hv@ValueAtRandomPoints == max(sp_hv@ValueAtRandomPoints) )
  niche_position = sp_hv@RandomPoints[max_prob,]
  hv_scale_positions[index,] = c(sp_name, niche_position)
}

# organizing hv volumes
spp_hvolumes = data.frame(all_spp_names, spp_hvolumes)
colnames(spp_hvolumes) = c("species", "hvolume")

# organizing scale values
names(hv_spp_list) = all_spp_names
hv_scale_values = rep(0, 5)
for (sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  species = rep(sp_name, 100)
  sp_values = hv_spp_list[[index]]@RandomPoints
  rand_n = sample(nrow(sp_values), 100)
  rand_sp_values = sp_values[rand_n,]
  hv_scale_values= rbind(hv_scale_values,data.frame(species,rand_sp_values) )
}
hv_scale_values = hv_scale_values[-1,]

# exporting
write.table(hv_scale_values, "3_hypervolume_inference/hv_scale_values.csv", sep=',', quote=F, row.names=F)
write.table(hv_scale_positions, "3_hypervolume_inference/hv_scale_positions.csv", sep=',', quote=F, row.names=F)
write.table(spp_hvolumes, "3_hypervolume_inference/spp_hvolumes.csv", sep=',', quote=F, row.names=F)

### setting regimes
spp_hvolumes = read.table("3_hypervolume_inference/spp_hvolumes.csv", sep=',', h=T)
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)
regimes = data.frame(spp_geographic_distribution, spp_hvolumes$hvolume)
colnames(regimes)[3] = "hvolume"

# exporting
write.table(regimes, "3_hypervolume_inference/regimes_hvolumes.csv", sep=',', quote=F, row.names=F)

############################# calculating differences and effects #################
regimes_hvolumes = read.table("3_hypervolume_inference/regimes_hvolumes.csv", header =T, sep=",",  na.strings = "NA", fill=T)

state_overall_hv = aggregate(regimes_hvolumes$hvolume, list(regimes_hvolumes$state), median)
state_iqr_hv = aggregate(regimes_hvolumes$hvolume, list(regimes_hvolumes$state), IQR )

(state_overall_hv$x[1] - state_overall_hv$x[2])/ state_overall_hv$x[2]
(state_overall_hv$x[1] - state_overall_hv$x[3])/ state_overall_hv$x[3]

state_iqr_rao

############################## back transforming values ###########################

#loading z-scale values
hv_scale_values = read.table("3_hypervolume_inference/hv_scale_values.csv", header =T, sep=",",  na.strings = "NA", fill=T)
scale_values= hv_scale_values[,2:5]

#back transforming each column
env_values = c()
for (j in 1:ncol(scale_values)){
  env_values = cbind(env_values, c((ras_sds[j]* scale_values[,j]) + ras_means[j]) )
}

# organizing
hv_env_values = data.frame(hv_scale_values$species, env_values)
colnames(hv_env_values) = colnames(hv_scale_values)

# getting iqr scale values
iqr = function(x){
  q= quantile(x, probs = c(0.25, 0.75) )
  iqr = as.numeric(q[2] - q[1])
  return(iqr)
}
iqr_env_values = aggregate(hv_env_values[,2:5], by=list(hv_env_values[,1]), iqr)
colnames(iqr_env_values)[1] = "species"

#exporting
write.table(hv_env_values, "3_hypervolume_inference/hv_env_values.csv", sep=',', quote=F, row.names=F)
write.table(iqr_env_values, "3_hypervolume_inference/iqr_env_values.csv", sep=',', quote=F, row.names=F)

