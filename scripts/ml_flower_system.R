############################## loading libraries ##############################

if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("randomForest")) install.packages("randomForest"); library("randomForest")
if (!require("caret")) install.packages("caret"); library("caret")
if (!require("Metrics")) install.packages("Metrics"); library("Metrics")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("phytools")) install.packages("phytools"); library("phytools")

################################ loading data ###############################

### load mcc phylogeentic tree
mcc_phylo = read.tree("2_comparative_analyses/pruned_mcc_phylo.nwk")

### pollination and mating system data
pmd = read.table("0_data/pollination_mating_data.csv", sep= ",", header= T, fill= T)
full_names = unique(paste0("Miconia ", pmd$species))

### trait data
trait_path = "0_data/MelastomaTraits (Jun 30, 2023).xlsx"
trait  = readxl::read_xlsx(path= trait_path, sheet= 3 )

############################## evaluating some patterns #######################

ms_summary = pmd %>% 
  reframe(
    n_ms = sum( !is.na(hand_self_fruiting_perc) & !is.na(bagged_fruiting_perc) ),
    n_out = sum( hand_self_fruiting_perc ==0  & 
                   bagged_fruiting_perc==0 , na.rm = T),
    n_self_apo = sum( hand_self_fruiting_perc >  0 &
                        bagged_fruiting_perc >  0, na.rm = T),
    n_self = sum( hand_self_fruiting_perc >  0 &
                    bagged_fruiting_perc ==  0, na.rm = T), 
    n_apo = sum( hand_self_fruiting_perc ==  0 & 
                   bagged_fruiting_perc >  0, na.rm = T)
  )


################# processing pollinating and mating systems ####################

large_pores = c("cremanium", "chaenopleura","glossocentrum", "hypoxanthus")

### classificando sistemas de polinização e acasalamento
system = pmd  %>% 
  mutate(
    species = paste0("Miconia ", species),
    ps = case_when(
      bees == 1 & other_insects == 0 & vertebrates == 0   ~ "specialist",
      bees == 1 & (other_insects == 1 | vertebrates == 1) ~ "generalist",
      TRUE                                                ~ NA
    ),
    ms1 = case_when(
      bagged_fruiting_perc >  0 &
        hand_self_fruiting_perc > 0  ~ bagged_fruiting_perc,
      bagged_fruiting_perc ==  0 &
        hand_self_fruiting_perc > 0  ~ NA
    ),
    ms2 = case_when(
      hand_self_fruiting_perc == 0 & bagged_fruiting_perc == 0 ~ "outcrosser",
      hand_self_fruiting_perc >  0 ~ "selfer",
      TRUE                         ~ NA
    ),
    pore_size = case_when(
      Cogniaux %in% large_pores    ~ 1,
      !Cogniaux %in% large_pores   ~ 0,
      TRUE                         ~ NA
    ),
  ) %>% 
  select(species, pore_size, nectar, ps, ms1, ms2)


########################## processing floral traits ############################

### Melastom traits
flw1 = trait %>% 
  filter(Genus == "Miconia" & Species %in% full_names) %>% 
  mutate(
    species = Species,
    flower_size = Style_length_mean,
    stamen_dim = (Stamen_length_antesepal_mean - Stamen_length_antepetal_mean)/Stamen_length_antesepal_mean,
    herkogamy = Style_length_mean - Stamen_length_antepetal_mean
  ) %>% 
  select(species, flower_size, stamen_dim, herkogamy) %>% 
  filter(!is.na(flower_size))


####################### join floral traits and system data #####################

### join floral traits and system data
data = system %>% 
  plyr::join(y = flw1, type = "left", by = "species") 

### ML vars
vars = c("nectar","flower_size", "pore_size", "stamen_dim", "herkogamy", "target")

### data for pollination system
data_ps = data %>% 
  filter(!is.na(ps) & !is.na(flower_size) ) %>% 
  mutate(target = factor(ps, levels = c("generalist", "specialist") )) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(vars) )
 
### data for mating system 1
data_ms1 = data %>% 
  filter(!is.na(ms1) & !is.na(flower_size)) %>% 
  mutate(target = ms1) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(vars ) ) 

### data for mating system 2
data_ms2 = data %>% 
  filter(!is.na(ms2) & !is.na(flower_size)) %>% 
  mutate(target = factor(ms2, levels = c("selfer", "outcrosser"))) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(vars) ) 

### shuffling function
shuffle = function(data){
  set.seed(41)
  sdata = data[sample(1:nrow(data)), ]
  return(sdata)
}

### shuffling
data_ms1 = shuffle(data = data_ms1)
data_ms2 = shuffle(data = data_ms2)

################################## Random Forest ##############################

### function to set train and test
set_train_test = function(data){
  list = list()
  set.seed(41) 
  ind = sample(2, size = nrow(data), replace = T, prob = c(0.6, 0.4))
  list$train = data[ind==2,]
  list$test = data[ind==1,]
  return(list)
}

target = "ms1"
nectar_exclude = F

### choose which dataset
if(target == "ps"){
  dt = data_ps
}
if(target == "ms1"){
  dt = data_ms1
}
if(target == "ms2"){
  dt = data_ms2
}

### exclude nectar?
if(nectar_exclude){
  dt = dt %>%
    select(!nectar)
}

### separating train and test
list = set_train_test(data = dt)
train = list$train
test = list$test

### training model 
rf = randomForest(formula = target ~ . , 
                  data = train, 
                  proximity = T) 

### predictions
pred = predict(rf, test)

### confusion matrix
if(target == "ps"){
  confusionMatrix(pred, test$target)
}
if(target == "ms1"){
  mae(pred, test$target)
}

### importance list
importance(rf)

### improtance plot
varImpPlot(rf,
           sort = F,
           main = paste0("Variable Importance ", target))

############################## dummy model ###################################

dummy_model = function(train, test, target){
  values = train[[target]]
  set.seed(41) 
  pred = values[sample(nrow(test), replace= F)]
  return(pred)
}

### dummy evaluation 
if (target == "ps"){
  dummy_pred = dummy_model(train = train, test = test, target = "target")
  confusionMatrix(dummy_pred, test$target)
}
if (target == "ms1"){
  dummy_pred = dummy_model(train = train, test = test, target = "target")
  mae(dummy_pred, test$target)
}

############################### predicting classes ###########################

### load flower pc scores
spp_cogn = read.table("0_data/my_spp_cogniaux.csv", sep=",", h=T)

### load flower pc scores
center_flower_df = read.table("1_flower_analyses/center_flower_df.csv", sep=",", h=T)

### getting pore size from sections
spp_pore = spp_cogn %>% 
  mutate(pore_size = case_when(
    Cogniaux %in% large_pores    ~ 1,
    !Cogniaux %in% large_pores   ~ 0,
    TRUE                         ~ NA
    )
  ) %>% 
  select(species, pore_size)

### organizing features
my_traits = center_flower_df %>% 
  plyr::join(y = spp_pore, type = "left", by = "species") %>% 
  select(flower_size, pore_size, stamen_dim, herkogamy)

###
rf_model = randomForest(formula = target ~ . , 
                  data = dt, 
                  proximity = T) 

### orgnaizing into data frame 
my_pred = predict(rf_model, my_traits)
geo_state = center_flower_df$geo_state
species = center_flower_df$species
pred_df = data.frame(geo_state, species, my_pred)

### table name
dir_name = paste0("1_flower_analyses/pred_class_",target,".csv")
### exporting
write.table(pred_df, dir_name, sep=",", row.names = F)

############################# plotting ps predictions ############################

### load flower pc scores
pred_df = read.table("1_flower_analyses/pred_class_ps.csv", sep=",", h=T)

### testing differences
tab = table(pred_df$geo_state, pred_df$my_pred)
chisq.test(tab)

### summarizing
pred_df= pred_df %>% 
  group_by(geo_state) %>% 
  reframe(generalist = sum(my_pred == "generalist"),
         specialist = sum(my_pred == "specialist"),
         ) %>% 
  pivot_longer(cols = c(generalist, specialist), names_to = "pollination" ) %>% 
  group_by(geo_state) %>% 
  mutate(perc = 100* value/sum(value) )

### plot param
axis_title_size = 8
x_text_size = 7
y_text_size = 7
legend_text_size = 6
legend_key_size = 0.4

### plotting
pred_plot = ggplot(data = pred_df) + 
  
  geom_bar(aes(x = geo_state, 
               y = perc, 
               fill = pollination,
               color = geo_state),
           stat = "identity", 
           width = 0.9,
           linewidth = 0.8,
           alpha = 0.75) +
  
  scale_x_discrete(labels=c("AF" = "AF-endemic", 
                            "other" = "non-endemic"))+
  
  scale_fill_manual(values = c("gray", "black") )+
  
  scale_color_manual(values = c("#1E88E5","#D81B60") )+
  
  xlab("geographic distribution") +
  
  ylab("relative frequency (%)\n of pollination systems") +
  
  guides(fill=guide_legend(title="")) +
  guides(color = "none")+
  
  theme(panel.background=element_rect(fill="white"),
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=axis_title_size, face="bold"),
        axis.text.x= element_text(size= x_text_size),
        axis.text.y = element_text(size=y_text_size, angle = 0),
        legend.position = "bottom",
        legend.text = element_text(size= legend_text_size),
        legend.key = element_blank(),
        legend.key.size = unit(legend_key_size, 'cm'))

# export plot
tiff("3_graphs/pred_plot_ps.tiff", units="cm", width=7, height=6.5, res=600)
  print(pred_plot)
dev.off()

############################# plotting ms predictions ############################

### load flower pc scores
pred_df = read.table("1_flower_analyses/pred_class_ms1.csv", sep=",", h=T)

### exploring differences
pred_df %>% 
  group_by(geo_state) %>% 
  reframe(mean(my_pred), sd(my_pred), median(my_pred))

### testing differences
my_pred = pred_df$my_pred
geo_state = pred_df$geo_state
names(my_pred) = pred_df$species
names(geo_state) = pred_df$species

paov = phylANOVA(tree = mcc_phylo, 
          x = geo_state ,
          y = my_pred,
          posthoc = F)

### plot param
axis_title_size = 8
x_text_size = 7
y_text_size = 7
legend_text_size = 6
legend_key_size = 0.4

### plotting
pred_plot = ggplot(data = pred_df,
                   aes(x=geo_state, 
                       y=my_pred, 
                       fill=geo_state)) + 
  
  geom_point(aes(color=geo_state),
             position = position_jitter(width = 0.07), 
             size = 1.2, 
             alpha = 0.65) +
  
  geom_boxplot(width = 0.5, 
               outlier.shape = NA,
               alpha = 0.25)+
  
  scale_fill_manual(values=c("#1E88E5","#D81B60"))+
  scale_colour_manual(values=c("#1E88E5","#D81B60"))+
  scale_x_discrete(labels=c("AF" = "AF-endemic", 
                            "other" = "non-endemic")) +
  
  xlab("geographic distribution") +
  
  ylab("relative frequency (%)\n of self pollination") +
  
  guides(fill=guide_legend(title="")) +
  guides(color = "none")+
  
  theme(panel.background=element_rect(fill="white"),
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=axis_title_size, face="bold"),
        axis.text.x= element_text(size= x_text_size),
        axis.text.y = element_text(size=y_text_size, angle = 0),
        legend.position = "bottom",
        legend.text = element_text(size= legend_text_size),
        legend.key = element_blank(),
        legend.key.size = unit(legend_key_size, 'cm'))

# export plot
tiff("3_graphs/pred_plot_ms.tiff", units="cm", width=7, height=6.5, res=600)
  print(pred_plot)
dev.off()
