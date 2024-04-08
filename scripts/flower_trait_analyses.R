### laoding libraries
if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("dplyr")) install.packages("dplyr"); library("dplyr")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr"); library("ggpubr")
if (!require("PupillometryR")) install.packages("PupillometryR"); library("PupillometryR") 
if (!require("reshape2")) install.packages("reshape2"); library("reshape2")

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

######################## plotting floral traits ################################

### loading species' flower traits
center_flower_df= read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)

### my colors
mycols = c( "#1E88E5", "#D81B60")
names(mycols) = c("AF", "other")


trait_plot = ggplot(data= center_flower_df, aes(x=geo_state, 
                                   y= herkogamy, 
                                   fill=geo_state)) +
  geom_point(aes(color=geo_state),
             size = 0.75, 
             alpha = 0.65,
             position = position_jitter(width = 0.07)) +
  
  geom_boxplot(width = 0.4, 
               alpha = 0.25,
               linewidth = 0.2,
               outlier.shape = NA)+
  
  scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "non-endemic")) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  
  xlab("geographic distribution")+ 
  ylab("Herkogamy")+
  
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"), 
        axis.title=element_text(size=6,face="bold"), 
        axis.text=element_text(size=4), 
        legend.position = "none") 

# export plot
tiff("3_graphs/herkogamy_plot.tiff", 
     units="cm", width=3.75, height=3, res=600)
  print(trait_plot)
dev.off()

