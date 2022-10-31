############################### fitting pgls models ###########################


### setting common variables
# choose phylogenetic tree
phylo = mcc_phylo

# set predictor variable
predictor = geo_state

# set dataframe with response variables
response_df = data.frame(flower_proxy[,-c(1:2)])
naming_vector = flower_proxy$species

### on the relationship
#flower_size = linear
#herkogamy = linear
#pore_size = log
#sitgma_size = log

log_these = c( "pore_size", "sitgma_size")

for (j in 1:ncol(response_df) ){
  ### picking one response variable
  response_name = colnames(response_df)[j]
  ### setting output dir
  # check if output dir exists
  dir_check = dir.exists(paths=paste("2_comparative_analyses/PGLS/", response_name, sep="") )
  # create output dir if not created yet
  if (dir_check == FALSE){
    dir.create(path= paste("2_comparative_analyses/PGLS/", response_name, sep=""), showWarnings = , recursive = FALSE, mode = "0777")
  }
  # ouput dir
  out_dir = paste("2_comparative_analyses/PGLS/", response_name, sep="")
  ### setting result objects
  best_evo_model = c("model", "aicc")
  intercept = c()
  angular = c()
  t_value = c()
  p_value = c()
  shapiro_w = c()
  shapiro_p = c()
  r2 = c()
  ### applying leave-one-out
  for (i in 1:nrow(response_df)){
    # droping n observatio in the predictor
    iter_predictor = predictor[-i]
    # selecting response variable
    iter_response = response_df[-i,j]
    iter_names = naming_vector[-i]
    names(iter_response) = iter_names
    # need transformation?
    if (response_name %in% log_these){iter_response = log(iter_response)}
    # initial values 
    init_bm = sd(iter_response)
    init_ou = mean(iter_response) - min(iter_response)
    # drop species from phylogenetic tree
    to_drop = naming_vector[i]
    iter_phylo = drop.tip(phylo, to_drop)
    # fitting models to flower traits
    fit_bm = fitContinuous(phy= iter_phylo, iter_response,  model="BM")
    fit_ou = fitContinuous(phy= iter_phylo, iter_response,  model="OU")
    # chosing best-fit model for response variable by aicc
    if (fit_bm$opt$aicc < fit_ou$opt$aicc){
      cor_mtx = corBrownian(value=init_bm, phy= iter_phylo, form=~iter_names)
      best_evo_model = rbind(best_evo_model, c("BM", fit_bm$opt$aicc))
    } else {
      if (fit_bm$opt$aicc - fit_ou$opt$aicc < 2){
        cor_mtx = corBrownian(value=init_bm, phy= iter_phylo, form=~iter_names)
        best_evo_model = rbind(best_evo_model, c("BM", fit_bm$opt$aicc))
      } else {
        cor_mtx = corMartins(value= init_ou, phy = iter_phylo, fixed = T, form=~iter_names)
        best_evo_model = rbind(best_evo_model, c("OU", fit_ou$opt$aicc))
      }
    }
    # fitting pgls
    fit_gls = gls(iter_response ~ iter_predictor, correlation=cor_mtx,  method = "REML")
    summary_gls = summary(fit_gls)
    # taking coefficients and test
    gls_test_table = summary_gls$tTable
    intercept = c(intercept, summary_gls$coefficients[1])
    angular = rbind(angular, summary_gls$coefficients[-1])
    t_value = c(t_value, gls_test_table[,3][2])
    p_value = rbind(p_value, gls_test_table[,4])
    # checking residuals
    res = resid(fit_gls)[1:length(iter_response)]
    shapiro = shapiro.test(res)
    shapiro_w = c(shapiro_w, shapiro$statistic)
    shapiro_p = c(shapiro_p, shapiro$p.value)
    # calculating R model fit
    ssr = sum(res^2)
    sst =  sum((iter_response- mean(iter_response))^2)
    r2 = c(r2, c(1 - (ssr/sst)) )
    print(paste("round:", as.character(i)))
  }
  print(paste("response done:", response_name))
  # organizing model evalution and pgls fit
  pgls_stats = data.frame(intercept, angular, t_value, p_value, r2, shapiro_w, shapiro_p)
  # exporting
  write.table(best_evo_model, paste(out_dir,"/best_evo_model.csv", sep=""), row.names=F, quote=F, sep=",")
  write.table(pgls_stats, paste(out_dir,"/pgls_stats.csv", sep=""), row.names=F, quote=F, sep=",")
}
