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
    ms = case_when(
      hand_self_fruiting_perc == 0 & bagged_fruiting_perc == 0 ~ "dependent",
      hand_self_fruiting_perc >  0 | bagged_fruiting_perc >  0 ~ "independent",
      TRUE                                                     ~ NA
    ),
    
    anther_pore = case_when(
      Cogniaux %in% large_pores    ~ 1,
      !Cogniaux %in% large_pores   ~ 0,
      TRUE                         ~ NA
    ),
  ) %>% 
  select(species, ps, ms, anther_pore, nectar)


### características florais
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

### juntando flores e caracteríticas
data = flw %>% 
  plyr::join(y = system, type = "left", by = "species")

### data for pollination syetem
data_ps = data %>% 
  filter(!is.na(ps)) %>% 
  mutate(target = factor(ps)) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(c("nectar","flower_size", "dimetrism","anther_pore", "herkogamy", "target") ) )
 
### data for mating syetem
data_ms = data %>% 
  filter(!is.na(ms)) %>% 
  mutate(target = factor(ms)) %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(any_of(c("nectar","flower_size", "dimetrism", "anther_pore", "herkogamy", "target") ) ) 


################################## Random Forest ##############################

set.seed(42) 

set_train_test = function(data){
  list = list()
  ind = sample(2, nrow(data), replace = TRUE, prob = c(0.6, 0.4))
  list$train = data[ind==1,]
  list$test = data[ind==2,]
  return(list)
}

data1 = data_ps
data2 = data_ps %>% 
  select(!nectar)

list = set_train_test(data = data1 )
train = list$train
test = list$test

rf = randomForest(formula = target ~ . , data = train, proximity = T) 

p1 <- predict(rf, test)
confusionMatrix(p1, test$target)

varImpPlot(rf,
           sort = T,
           n.var = 4,
           main = "Variable Importance")
importance(rf)

MDSplot(rf, train$target)
