
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

### loading all flower traits
ftraits = read.table("0_data/flower_trait_matrix.csv", sep=",", h=T)

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### selfing-diagnostic traits
# flower size
flower_size = ftraits$hypathium_openning
# herkogamy
major_stamen_height = ftraits$major_filet_height + ftraits$major_anther_height
minor_stamen_height = ftraits$minor_filet_height + ftraits$minor_anther_height
herkogamy_major = ftraits$style_height - major_stamen_height
herkogamy_minor = ftraits$style_height - minor_stamen_height
# stigma e pore size
stigma_size = ftraits$stigma_width
pore_size = ftraits$pore_long_section

### define geopraphic states 
high_ths = 0.90
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "AF"
geo_states[af_percentage <= low_ths] = "other"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "AFother"
names(geo_states) = spp_count_domain$species


### flower dataframe
flower_df = data.frame(flower_size, herkogamy_minor, stigma_size, pore_size)


################################# PCA analyses #################################

data_df = flower_df

### observed patterns 
pca = prcomp(data_df, center = F)
stdev = pca$sdev / sum(pca$sdev)
pc_numbers = paste("pc", as.character(1:length(pca$sdev)), sep='')
names(stdev) = pc_numbers

### boots values
boot_stdev_df = data.frame(matrix(0, nrow=1, ncol=4))
boot_rotation_list = vector('list', 99)

for (i in 1:99){
  boot_num = sample(x=1:nrow(data_df), size = nrow(data_df), replace = T)
  boot_data = data_df[boot_num,]
  boot_pca = prcomp(boot_data, center = F)
  boot_rotation_list[[i]] = boot_pca$rotation
  boot_stdev = boot_pca$sdev/ sum(boot_pca$sdev)
  boot_stdev_df = rbind(boot_stdev_df, boot_stdev)
}
boot_stdev_df = boot_stdev_df[-1,]

### testing if observed values are greater than random values
pvalues = c()
for(i in 1:length(stdev) ){
  obs_stdev = stdev[i]
  boot_distibution = boot_stdev_df[,i]
  greater = sum(obs_stdev > boot_distibution)
  total =  length(boot_distibution) +1 
  p = 1 - (greater / total)
  pvalues = c(pvalues, p)
}
names(pvalues) = pc_numbers

# summary
boot_min = apply(boot_stdev_df, MARGIN = 2, min)
boot_median = apply(boot_stdev_df, MARGIN = 2, median)
boot_max = apply(boot_stdev_df, MARGIN = 2, max)
stdev_summary = data.frame(pc_numbers, boot_min, boot_median, boot_max, stdev) 
summary(stdev_summary)

# ploting
tiff("1_flower_analyses/stdev_pca.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= stdev_summary, aes(x=pc_numbers, y=stdev ) ) +
  geom_errorbar(aes(ymin=boot_min, ymax=boot_max), col="gray", size=1, width=0.01)+
  geom_point(aes(y=boot_median), color="gray", size = 1.5) +
  geom_point(aes(y=stdev), color="red", size = 1) +
  xlab("principal component")+ ylab("standard deviation")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=12))
dev.off()

# exporting test values
write.table(pvalues, "1_flower_analyses/stdev_boot_test.csv", sep=",", quote=F, col.names=T)
