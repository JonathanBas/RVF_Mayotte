library("rjags")
library("bayesplot")
library("ggpubr")
load.module("dic")

source(file = "./import_data.R")

to_mcmc.list = function(mcarray_obj, to_burn, which_par = "all"){
  nchains = length(as.mcmc.list(mcarray_obj[["deviance"]]))
  
  if(which_par == "all"){
    var_mod = names(mcarray_obj)[which(!(names(mcarray_obj) %in% c("pD")))]
  }else{
    var_mod = names(mcarray_obj)[which(!grepl("pD|igg_t|inc_true|inc_obs", names(mcarray_obj)))]
  }
  
  Mch = as.list(rep(NA, nchains))
  for (chain_i in 1:nchains){
    for (na in var_mod){
      post_var_chain_i = as.mcmc.list(mcarray_obj[[na]])[[chain_i]]
      Mch[[chain_i]] = cbind(Mch[[chain_i]], post_var_chain_i)
    }
    included_in_chain = (1:nrow(Mch[[chain_i]]))[! (1:nrow(Mch[[chain_i]])) %in% to_burn]
    Mch[[chain_i]] = mcmc(Mch[[chain_i]][included_in_chain,-1])
  }
  output_Mch = as.mcmc.list(Mch)
  output_combined_Mch = do.call("rbind", Mch)
  
  list(output_var_mod = var_mod, output_Mch=output_Mch, output_combined_Mch=output_combined_Mch)
}

n_simu <- 5000
num_aggr <- "aggr_geo3"
n_zones <- list(aggr_geo1 = 1,
                aggr_geo3 = 2)
popul = list(aggr_geo1 = c("Mayotte" = 144262),
             aggr_geo3 = c("Central communes" = 68189, "Outer communes" = 76073))

dat1_filt <- filter(dat1, !is.na(dat1[,num_aggr]))
datmod <- as.data.frame(dat1_filt %>% group_by(geo = dat1_filt[,num_aggr], week) %>% summarise(nb_igg_plus = sum(type=="sero" & inf=="anc"),
                                                                                                  nb_test = sum(type=="sero"),
                                                                                                  nb_cas = sum(type=="cas")))
list_wks <- c(paste0("2018-", 41:52), paste0("2019-0", 1:9), paste0("2019-", 10:40))
list_geo <- unique(dat1_filt[,num_aggr])
reg_wks <- expand.grid(geo = list_geo, week = list_wks)

datmod <- merge(datmod, reg_wks, by=c("geo","week"), all=T)
datmod[is.na(datmod)] <- 0

inc_obs <- datmod$nb_cas
igg_obs <- datmod$nb_igg_plus
nb_test_antic <- datmod$nb_test

n_times <- length(list_wks)
tvec <- data.frame(t = (0:(nrow(datmod)-1)) %% n_times+1, wks = datmod$week, used = ifelse(datmod$week %in% paste0("2019-", 18:52), F, T))
max_lag_ig <- min(as.data.frame(datmod %>% group_by(geo) %>% summarise(val = which(nb_test != 0)[1] - 1))$val)

rm(jags)
jags <- jags.model('./mod_RVF_bug_fit_reg.bug',
                   data = list('n_zones'=n_zones[[num_aggr]], 'popul'=popul[[num_aggr]], 'n_times'=n_times, 'inc_obs'=inc_obs, 'igg_t'=igg_obs*tvec$used, 'nb_test_antic'=nb_test_antic*tvec$used, 'max_lag_ig'=max_lag_ig),
                   n.chains = 3,
                   n.adapt = 0)

rm(samples)
chrono_start <- proc.time()[3]
samples <- jags.samples(jags,
                        c('tau','serop_igg_pre','param1_log','param1','param2','param3','param4','deviance','serop_igg_pre_tot','param1_tot','serop_igg_post','serop_igg_post_tot'),
                        200000,
                        thin = 200)
print(paste("Calculation time:", proc.time()[3]-chrono_start, "secondes"))
date_time <- Sys.time() ; date_time <- gsub(":", "-", date_time)
save(samples, file=paste0("./res/samples ", date_time, ".Rdata"))

# load(file=paste0("./res/samples ", date_time, ".Rdata"))
transformed_chain = to_mcmc.list(mcarray_obj = samples, to_burn = 1:150, which_par = "not_all")
var_mod = transformed_chain[["output_var_mod"]]
Mch = transformed_chain[["output_Mch"]]
combined_Mch = transformed_chain[["output_combined_Mch"]]

plot(Mch) ; par(mfrow = c(1, 1))
vec_med <- summary(Mch)[[2]][,"50%"]
tab_inter <- HPDinterval(as.mcmc(combined_Mch))
if(!any(names(vec_med) != rownames(tab_inter))){
  tab_param <- cbind(round(vec_med,4), round(tab_inter,4))
}
print(tab_param)
write.table(tab_param, file = paste0("./res/tab_param ", date_time, ".csv"), sep=";", col.names=NA, row.names=T)

par(mfrow = c(3, 3)) ; for(i in 1:ncol(combined_Mch)){acf(combined_Mch[,i], lag.max = 30, ylim = c(-1,1), main = colnames(combined_Mch)[i])} ; par(mfrow = c(1, 1))
effectiveSize(Mch)
gelman.plot(Mch)


