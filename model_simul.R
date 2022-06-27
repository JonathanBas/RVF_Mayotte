
source(file = "./model_fit.R")

if(num_aggr == "aggr_geo1"){
  
  tau <- as.vector(samples[["tau"]])
  igg_pre <- as.vector(samples[["serop_igg_pre"]])
  param1 <- as.vector(samples[["param1"]])
  param2 <- as.vector(samples[["param2"]])
  param3 <- as.vector(samples[["param3"]])
  param4 <- as.vector(samples[["param4"]])

  sim_inc_true <- matrix(NA, n_simu, n_times)
  sim_inc_obs <- matrix(NA, n_simu, n_times)
  sim_igg <- matrix(NA, n_simu, n_times)
  serop_igg <- matrix(NA, n_simu, n_times)
  
  for(sim in 1:n_simu){
    iter_num <- sample(1:length(samples[[1]]), 1)
    tau_iter <- tau[iter_num]
    igg_pre_iter <- igg_pre[iter_num]
    param1_iter <- param1[iter_num]
    param2_iter <- param2[iter_num]
    param3_iter <- param3[iter_num]
    param4_iter <- param4[iter_num]

    for (t in 1:max_lag_ig){
      
      sim_inc_true[sim,t] <- (param1_iter/sqrt(2*3.142*param3_iter)) *exp(-(log(t)-param2_iter)^2/(2*param3_iter))/t
      sim_inc_obs[sim,t] <- rbinom(n=1, size=round(sim_inc_true[sim,t]), prob=tau_iter)
      sim_igg[sim,t] <- rbinom(n=1, size=nb_test_antic[t], prob=igg_pre_iter)
      serop_igg[sim,t] <- igg_pre_iter
    }
    for (t in (max_lag_ig+1):n_times){
      
      sim_inc_true[sim,t] <- (param1_iter/sqrt(2*3.142*param3_iter)) *exp(-(log(t)-param2_iter)^2/(2*param3_iter))/t
      sim_inc_obs[sim,t] <- rbinom(n=1, size=round(sim_inc_true[sim,t]), prob=tau_iter)
      
      dyn_igg <- rep(NA,t)
      dyn_igg[1] <- (t >= 1+param4_iter) * sim_inc_true[sim,1]
      for(i in 2:t){
        dyn_igg[i] <- dyn_igg[i-1] + (t >= i+param4_iter) * sim_inc_true[sim,i]
      }
      dyn_igg_prop <- dyn_igg /popul[[num_aggr]]
      
      sim_igg[sim,t] <- rbinom(n=1, size=nb_test_antic[t], prob=min(0.999, igg_pre_iter + dyn_igg_prop[t]))
      
      serop_igg[sim,t] <- igg_pre_iter + dyn_igg_prop[t]
    }
  }
  
  sim_inc_obs_val <- apply(X=sim_inc_obs, FUN=median, MARGIN=2, na.rm=T)
  sim_inc_obs_max <- apply(X=sim_inc_obs, FUN=quantile, MARGIN=2, probs=0.975)
  sim_inc_obs_min <- apply(X=sim_inc_obs, FUN=quantile, MARGIN=2, probs=0.025)
  
  serop_igg_val <- apply(X=serop_igg, FUN=median, MARGIN=2, na.rm=T)
  serop_igg_max <- apply(X=serop_igg, FUN=quantile, MARGIN=2, probs=0.975)
  serop_igg_min <- apply(X=serop_igg, FUN=quantile, MARGIN=2, probs=0.025)
  
  sim_igg_val <- apply(X=sim_igg, FUN=median, MARGIN=2, na.rm=T)
  sim_igg_max <- apply(X=sim_igg, FUN=quantile, MARGIN=2, probs=0.975)
  sim_igg_min <- apply(X=sim_igg, FUN=quantile, MARGIN=2, probs=0.025)
  
  sim_inc_true_val <- apply(X=sim_inc_true, FUN=median, MARGIN=2)
  sum(sim_inc_true_val)
  sim_inc_true_max <- apply(X=sim_inc_true, FUN=quantile, MARGIN=2, probs=0.975)
  sim_inc_true_min <- apply(X=sim_inc_true, FUN=quantile, MARGIN=2, probs=0.025)
  
  png(filename="./res/fig_2.png",pointsize=6,res=300,width = 24, height = 20, units = "cm")
  par(mfrow=c(2,2))
  
  igg_obs_inf <- sapply(X=1:n_times, FUN = function(z_iter){ifelse(nb_test_antic[z_iter]!=0, prop.test(igg_obs[z_iter], nb_test_antic[z_iter])$conf.int[1], NA)})
  igg_obs_sup <- sapply(X=1:n_times, FUN = function(z_iter){ifelse(nb_test_antic[z_iter]!=0, prop.test(igg_obs[z_iter], nb_test_antic[z_iter])$conf.int[2], NA)})
  grad_x_axis <- seq(1, n_times, 4)
  
  p1 = ggplot()
  p1 = p1 + geom_line(aes(x = 1:n_times, y = sim_inc_true_val), size=1)
  p1 = p1 + geom_line(aes(x = 1:n_times, y = sim_inc_true_max), linetype="dashed")
  p1 = p1 + geom_line(aes(x = 1:n_times, y = sim_inc_true_min), linetype="dashed")
  p1 = p1 + ggtitle(paste0(names(popul[[num_aggr]]), "\nIncident infections (reported and unreported)"))
  p1 = p1 + theme_bw()
  p1 = p1 + scale_x_continuous(labels = tvec$wks[grad_x_axis], breaks = grad_x_axis)
  p1 = p1 + theme(axis.text.x = element_text(angle = 90))
  p1 = p1 + xlab("Weeks") + ylab("Number of new infections")
  print(sum(sim_inc_true_val))
  print(sum(inc_obs))
  
  p2 = ggplot()
  p2 = p2 + geom_point(aes(x = 1:n_times, y = inc_obs))
  p2 = p2 + geom_line(aes(x = 1:n_times, y = sim_inc_obs_val), size=1, col="dodgerblue4")
  p2 = p2 + geom_line(aes(x = 1:n_times, y = sim_inc_obs_max), linetype="dashed", col="dodgerblue4")
  p2 = p2 + geom_line(aes(x = 1:n_times, y = sim_inc_obs_min), linetype="dashed", col="dodgerblue4")
  p2 = p2 + theme_bw()
  p2 = p2 + scale_x_continuous(labels = tvec$wks[grad_x_axis], breaks = grad_x_axis)
  p2 = p2 + theme(axis.text.x = element_text(angle = 90))
  p2 = p2 + ggtitle("\nIncident reported infections") + xlab("Weeks") + ylab("Number of new reported infections")
  
  p3 = ggplot() + theme_minimal()
  
  p4 = ggplot()
  p4 = p4 + geom_linerange(aes(x = 1:n_times, ymin = 100*igg_obs_inf, ymax = 100*igg_obs_sup))
  p4 = p4 + geom_point(aes(x = 1:n_times, y = 100*igg_obs/nb_test_antic, fill = tvec$used), shape=21)
  p4 = p4 + scale_fill_manual(values = c("TRUE"="black", "FALSE"="white"), labels = c("TRUE"="Fitting data", "FALSE"="Data not used\nfor fitting"), name=NULL)
  p4 = p4 + geom_line(aes(x = 1:n_times, y = 100*serop_igg_val), size=1, col="dodgerblue4")
  p4 = p4 + geom_line(aes(x = 1:n_times, y = 100*serop_igg_max), linetype="dashed", col="dodgerblue4")
  p4 = p4 + geom_line(aes(x = 1:n_times, y = 100*serop_igg_min), linetype="dashed", col="dodgerblue4")
  p4 = p4 + theme_bw() + ylim(0,40)
  p4 = p4 + scale_x_continuous(labels = tvec$wks[grad_x_axis], breaks = grad_x_axis)
  p4 = p4 + theme(axis.text.x = element_text(angle = 90),
                  legend.position = c(0.8, 0.8),
                  legend.direction = "vertical",
                  legend.background = element_rect(fill = "white", color = "black"))
  p4 = p4 + ggtitle("\nIgG seroprevalence") + xlab("Weeks") + ylab("Percentage of IgG+ people (%)")
  
  p = ggarrange(p1,p2,p3,p4, nrow=2, ncol=2, labels=c("A","B","","C"))
  plot(p)
  
  dev.off()
  par(mfrow=c(1,1))
  
}else{
  
  sim_inc_true <- list() ; for(z in 1:as.numeric(n_zones[num_aggr])){sim_inc_true[[z]] <- matrix(NA, n_simu, n_times)}
  sim_inc_obs <- list() ; for(z in 1:as.numeric(n_zones[num_aggr])){sim_inc_obs[[z]] <- matrix(NA, n_simu, n_times)}
  sim_igg <- list() ; for(z in 1:as.numeric(n_zones[num_aggr])){sim_igg[[z]] <- matrix(NA, n_simu, n_times)}
  serop_igg <- list() ; for(z in 1:as.numeric(n_zones[num_aggr])){serop_igg[[z]] <- matrix(NA, n_simu, n_times)}
  
  for(sim in 1:n_simu){
    iter_num <- sample(1:nrow(combined_Mch), 1)
    tau_iter <- combined_Mch[iter_num,"tau"]
    param4_iter <- combined_Mch[iter_num,"param4"]

    igg_pre_iter = param1_iter = param2_iter = param3_iter = c()
    for(z in 1:as.numeric(n_zones[num_aggr])){
      igg_pre_iter <- c(igg_pre_iter, combined_Mch[iter_num,paste0("serop_igg_pre[",z,"]")])
      param1_iter <- c(param1_iter, combined_Mch[iter_num,paste0("param1[",z,"]")])
      param2_iter <- c(param2_iter, combined_Mch[iter_num,paste0("param2[",z,"]")])
      param3_iter <- c(param3_iter, combined_Mch[iter_num,paste0("param3[",z,"]")])
    }
    
    for (z in 1:as.numeric(n_zones[num_aggr])){
      for (t in 1:max_lag_ig){
        
        sim_inc_true[[z]][sim,t] <- (param1_iter[z]/sqrt(2*3.142*param3_iter[z]))*exp(-(log(t)-param2_iter[z])^2/(2*param3_iter[z]))/t
        sim_inc_obs[[z]][sim,t] <- rbinom(n=1, size=round(sim_inc_true[[z]][sim,t]), prob=tau_iter)
        sim_igg[[z]][sim,t] <- rbinom(n=1, size=nb_test_antic[(z - 1)*n_times + t], prob=igg_pre_iter[z])
        serop_igg[[z]][sim,t] <- igg_pre_iter[z]
      }
      for (t in (max_lag_ig+1):n_times){
        
        sim_inc_true[[z]][sim,t] <- (param1_iter[z]/sqrt(2*3.142*param3_iter[z]))*exp(-(log(t)-param2_iter[z])^2/(2*param3_iter[z]))/t
        sim_inc_obs[[z]][sim,t] <- rbinom(n=1, size=round(sim_inc_true[[z]][sim,t]), prob=tau_iter)
        
        dyn_igg <- rep(NA,t)
        dyn_igg[1] <- (t >= 1+param4_iter) * sim_inc_true[[z]][sim,1]
        for(i in 2:t){
          dyn_igg[i] <- dyn_igg[i-1] + (t >= i+param4_iter) * sim_inc_true[[z]][sim,i]
        }
        dyn_igg_prop <- dyn_igg /popul[[num_aggr]][z]
        
        sim_igg[[z]][sim,t] <- rbinom(n=1, size=nb_test_antic[(z - 1)*n_times + t], prob=min(0.999, igg_pre_iter[z] + dyn_igg_prop[t]))
        
        serop_igg[[z]][sim,t] <- igg_pre_iter[z] + dyn_igg_prop[t]
      }
    }
  }
  
  sim_inc_obs_val <- lapply(X=sim_inc_obs, FUN=function(tab){return(apply(X=tab, FUN=median, MARGIN=2, na.rm=T))})
  sim_inc_obs_max <- lapply(X=sim_inc_obs, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.975))})
  sim_inc_obs_min <- lapply(X=sim_inc_obs, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.025))})
  
  serop_igg_val <- lapply(X=serop_igg, FUN=function(tab){return(apply(X=tab, FUN=median, MARGIN=2, na.rm=T))})
  serop_igg_max <- lapply(X=serop_igg, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.975))})
  serop_igg_min <- lapply(X=serop_igg, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.025))})
  
  sim_igg_val <- lapply(X=sim_igg, FUN=function(tab){return(apply(X=tab, FUN=median, MARGIN=2, na.rm=T))})
  sim_igg_max <- lapply(X=sim_igg, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.975))})
  sim_igg_min <- lapply(X=sim_igg, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.025))})
  
  sim_inc_true_val <- lapply(X=sim_inc_true, FUN=function(tab){return(apply(X=tab, FUN=median, MARGIN=2, na.rm=T))})
  lapply(X=sim_inc_true_val, FUN=sum)
  sim_inc_true_max <- lapply(X=sim_inc_true, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.975))})
  sim_inc_true_min <- lapply(X=sim_inc_true, FUN=function(tab){return(apply(X=tab, FUN=quantile, MARGIN=2, probs=0.025))})
  
  p = list()
  
  for(z in 1:as.numeric(n_zones[num_aggr])){

    time_iter_z <- ((z - 1)*n_times+1):((z - 1)*n_times+n_times)
    igg_obs_inf <- sapply(X=time_iter_z, FUN = function(z_iter){ifelse(nb_test_antic[z_iter]!=0, prop.test(igg_obs[z_iter], nb_test_antic[z_iter])$conf.int[1], NA)})
    igg_obs_sup <- sapply(X=time_iter_z, FUN = function(z_iter){ifelse(nb_test_antic[z_iter]!=0, prop.test(igg_obs[z_iter], nb_test_antic[z_iter])$conf.int[2], NA)})
    grad_x_axis <- seq(1, n_times, 4)
    
    p1 = ggplot()
    p1 = p1 + geom_line(aes(x = 1:n_times, y = sim_inc_true_val[[z]]), size=1)
    p1 = p1 + geom_line(aes(x = 1:n_times, y = sim_inc_true_max[[z]]), linetype="dashed")
    p1 = p1 + geom_line(aes(x = 1:n_times, y = sim_inc_true_min[[z]]), linetype="dashed")
    p1 = p1 + ggtitle(paste0(names(popul[[num_aggr]])[z], "\nIncident infections (reported and unreported)"))
    p1 = p1 + theme_bw() + ylim(0, max(unlist(sim_inc_true_max)))
    p1 = p1 + scale_x_continuous(labels = tvec$wks[time_iter_z][grad_x_axis], breaks = grad_x_axis)
    p1 = p1 + theme(axis.text.x = element_text(angle = 90))
    p1 = p1 + xlab(NULL) + ylab("Number of new infections")
    print(sum(sim_inc_true_val[[z]]))
    print(sum(inc_obs[time_iter_z]))
    
    p2 = ggplot()
    p2 = p2 + geom_point(aes(x = 1:n_times, y = inc_obs[time_iter_z]))
    p2 = p2 + geom_line(aes(x = 1:n_times, y = sim_inc_obs_val[[z]]), size=1, col="dodgerblue4")
    p2 = p2 + geom_line(aes(x = 1:n_times, y = sim_inc_obs_max[[z]]), linetype="dashed", col="dodgerblue4")
    p2 = p2 + geom_line(aes(x = 1:n_times, y = sim_inc_obs_min[[z]]), linetype="dashed", col="dodgerblue4")
    p2 = p2 + theme_bw()
    p2 = p2 + scale_x_continuous(labels = tvec$wks[time_iter_z][grad_x_axis], breaks = grad_x_axis)
    p2 = p2 + scale_y_continuous(breaks=seq(0,1000,2), limits=c(0, max(inc_obs, unlist(sim_inc_obs_max))))
    p2 = p2 + theme(axis.text.x = element_text(angle = 90))
    p2 = p2 + ggtitle("Incident reported infections") + xlab(NULL) + ylab("Number of new reported infections")
    
    p3 = ggplot() + theme_minimal()
    
    p4 = ggplot()
    p4 = p4 + geom_linerange(aes(x = 1:n_times, ymin = 100*igg_obs_inf, ymax = 100*igg_obs_sup))
    p4 = p4 + geom_point(aes(x = 1:n_times, y = 100*igg_obs[time_iter_z]/nb_test_antic[time_iter_z], fill = tvec$used[time_iter_z]), shape=21)
    p4 = p4 + scale_fill_manual(values = c("TRUE"="black", "FALSE"="white"), labels = c("TRUE"="Fitting data", "FALSE"="Data not used\nfor fitting"), name=NULL)
    p4 = p4 + geom_line(aes(x = 1:n_times, y = 100*serop_igg_val[[z]]), size=1, col="dodgerblue4")
    p4 = p4 + geom_line(aes(x = 1:n_times, y = 100*serop_igg_max[[z]]), linetype="dashed", col="dodgerblue4")
    p4 = p4 + geom_line(aes(x = 1:n_times, y = 100*serop_igg_min[[z]]), linetype="dashed", col="dodgerblue4")
    p4 = p4 + theme_bw() + ylim(0,40)
    p4 = p4 + scale_x_continuous(labels = tvec$wks[time_iter_z][grad_x_axis], breaks = grad_x_axis)
    p4 = p4 + theme(axis.text.x = element_text(angle = 90),
                    legend.position = c(0.8, 0.8),
                    legend.direction = "vertical",
                    legend.background = element_rect(fill = "white", color = "black"))
    p4 = p4 + ggtitle("IgG seroprevalence") + xlab("Weeks") + ylab("Percentage of IgG+ people (%)")
    
    p[[z]] = ggarrange(p1,p2,p4,ncol=1, align="v", labels=LETTERS[c(1,2,3)+(z-1)*3], label.x = 0.93)
    
  }
  
  p_glob = ggarrange(plotlist = p, ncol = as.numeric(n_zones[num_aggr]))
  png(filename="./res/fig_2.png",pointsize=6,res=300,width = 12*as.numeric(n_zones[num_aggr]), height = 28, units = "cm")
  plot(p_glob)
  dev.off()
  par(mfrow=c(1,1))
  
  pdf(file="./res/fig_2.pdf", width=4.5*as.numeric(n_zones[num_aggr]), height=10)
  print(p_glob)
  dev.off()
  
}
par(mfrow=c(1,1))

