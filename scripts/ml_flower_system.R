############################## loading libraries ##############################

if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("caret")) install.packages("caret"); library("caret")
if (!require("randomForest")) install.packages("randomForest"); library("randomForest")

################################ loading data ###############################

### pollination and mating system data
pmd = read.table("0_data/pollination_mating_data.csv", sep= ",", header= T, fill= T)
full_names = unique(paste0("Miconia ", pmd$species))

### trait data
trait_path = "0_data/MelastomaTraits (Jun 30, 2023).xlsx"
trait  = readxl::read_xlsx(path= trait_path, sheet= 3 )

############################## evaluating some patterns #######################

n_ms = pmd %>% 
  filter( !is.na(hand_self_fruiting_perc) & !is.na(bagged_fruiting_perc) ) %>% 
  nrow()

n_out = pmd %>% 
  filter( hand_self_fruiting_perc ==0  & bagged_fruiting_perc==0 ) %>% 
  nrow()

n_self_apo = pmd %>% 
  filter( hand_self_fruiting_perc >  0 & bagged_fruiting_perc >  0) %>% 
  nrow()

n_self = pmd %>% 
  filter( hand_self_fruiting_perc >  0 & bagged_fruiting_perc ==  0) %>% 
  nrow()

n_apo = pmd %>% 
  filter( hand_self_fruiting_perc ==  0 & bagged_fruiting_perc >  0) %>% 
  nrow()

############################### processing data ##############################

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
      hand_self_fruiting_perc == 0 & bagged_fruiting_perc == 0 ~ "dependent",
      hand_self_fruiting_perc >  0 | bagged_fruiting_perc >  0 ~ "independent",
      TRUE                                                     ~ NA
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
  select(species, pore_size, nectar, ps, ms1, ms2, ms3)


### floral traits
flw = trait %>% 
  filter(Genus == "Miconia" & Species %in% full_names) %>% 
  mutate(
    species = Species,
    flower_size = Style_length_mean,
    dimetrism = (Stamen_length_antesepal_mean - Stamen_length_antepetal_mean)/Stamen_length_antesepal_mean,
    herkogamy = Style_length_mean - Stamen_length_antepetal_mean
  ) %>% 
  select(species, flower_size, dimetrism, herkogamy) %>% 
  filter(!is.na(flower_size))

### join floral traits and system data
data = flw %>% 
  plyr::join(y = system, type = "left", by = "species")

### data for pollination syetem
data_ps = data %>% 
  filter(!is.na(ps)) %>% 
  mutate(target = factor(ps, levels = c("generalist", "specialist") )) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(c("nectar","flower_size", "dimetrism","pore_size", "herkogamy", "target") ) )
 
### data for mating system
data_ms1 = data %>% 
  filter(!is.na(ms1)) %>% 
  mutate(target = factor(ms1, levels = c("independent", "dependent"))) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(c("nectar","flower_size", "dimetrism", "pore_size", "herkogamy", "target") ) ) 

data_ms2 = data %>% 
  filter(!is.na(ms2)) %>% 
  mutate(target = factor(ms2, levels = c("selfer", "outcrosser"))) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(c("nectar","flower_size", "dimetrism", "pore_size", "herkogamy", "target") ) ) 


################################## Random Forest ##############################

### function to set train and test
set_train_test = function(data){
  list = list()
  set.seed(40) 
  ind = sample(2, size = nrow(data), replace = T, prob = c(0.4, 0.6))
  list$train = data[ind==2,]
  list$test = data[ind==1,]
  return(list)
}

### choose which dataset
dt = data_ps
dt = data_ms2

### exclude nectar?
dt = dt %>%
  select(!nectar)

### separating train and test
list = set_train_test(data = dt)
train = list$train
test = list$test

### training model 
rf = randomForest(formula = target ~ . , 
                  data = train, 
                  proximity = T) 

### predictions
pred <- predict(rf, test)

### confusion matrix
confusionMatrix(pred, test$target)

### improtance plot
varImpPlot(rf,
           sort = F,
           main = "Variable Importance")

### importance list
importance(rf)

### setting a dummy model
dummy = function(data, target){
  target = data[[target]]
  n = length(target)
  set.seed(40) 
  pred = sample(target, size = n)
  return(pred)
}

### dummy prediction
dummy_pred = dummy(data = test, target = "target")

### confusion matrix for dummy
confusionMatrix(dummy_pred, test$target)

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
  rename(dimetrism = stamen_dim) %>% 
  select(flower_size, pore_size, dimetrism, herkogamy)

### orgnaizing into data frame 
my_pred = predict(rf, my_traits)
geo_state = center_flower_df$geo_state
species = center_flower_df$species
pred_class_df = data.frame(geo_state, species, my_pred)

### exporting
write.table(pred_class_df, "1_flower_analyses/pred_class_mating2.csv", sep=",", row.names = F)

############################# comparing classes ################################

### load flower pc scores
pred_class = read.table("1_flower_analyses/pred_class_pollin.csv", sep=",", h=T)

### testing differences
tab = table(pred_class$geo_state, pred_class$my_pred)
chisq.test(tab)

### summarizing
pred_df= pred_class %>% 
  group_by(geo_state) %>% 
  reframe(generalist = sum(my_pred == "generalist"),
         specialist = sum(my_pred == "specialist"),
         ) %>% 
  pivot_longer(cols = c(generalist, specialist), names_to = "pollination" ) %>% 
  group_by(geo_state) %>% 
  mutate(perc = 100* value/sum(value) )

### plot param
axis_title_size = 8
x_text_size = 6
y_text_size = 6
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
tiff("3_graphs/pred_plot.tiff", units="cm", width=7, height=6.5, res=600)
  print(pred_plot)
dev.off()
