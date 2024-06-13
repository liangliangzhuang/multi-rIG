# ==========================================================
# ============ Major functions used in this paper ==========
# ==========================================================

# 1. Reparameterized Inverse Gaussian Distribution ====
# 1.1 Generate random numbers ====
rrIG = function(n1, lambda1, gamma1){
  re = rinvGauss(n1, nu = lambda1/gamma1, lambda = lambda1^2)
  return(re)
}
# 1.2 Generate density function ====
drIG = function(n1, lambda1, gamma1){
  re = dinvGauss(n1, nu = lambda1/gamma1, lambda = lambda1^2)
  return(re)
}
# 1.3 Generates distribution function ====
prIG = function(p, lambda1, gamma1){
  re = pinvGauss(p, nu = lambda1/gamma1, lambda = lambda1^2)
  return(re)
}


# 2 Lambda related functions ====
# Scen.I -> "pp", Scen.II -> "pl", Scen.III -> "lp", Scen.IV -> "ll"
# "p" means "power-law" form $\Lambda_k(t; \alpha_k, \beta_k) = \beta_k t^{\alpha_k}$
# "l" means "log-linear" form $\Lambda_k(t; \alpha_k, \beta_k) = \beta_k \left[ \exp (\alpha_k t) - 1 \right]$
# 2.1 Lambda functions in different scenarios ====
Lambda_fun = function(type, par0, m=m, bin = 1){
  alpha0 = par0[1]; beta0 = par0[2]
  alpha = par0[2+1:K]; beta = par0[2+K+1:K]
  t = seq(0,m,1/bin)
  Lambda = matrix(NA, length(t)-1, K)
  Lambda0 = numeric()
  if(type == "pp"){
    for(j in 1:(length(t)-1)){
      Lambda0[j] = beta0 * (t[j+1]^alpha0 - t[j]^alpha0)
      Lambda[j,] = beta * (t[j+1]^alpha - t[j]^alpha)
    }
  } else if(type == "pl"){
    for(j in 1:(length(t)-1)){
      Lambda0[j] = beta0 * (t[j+1]^alpha0 - t[j]^alpha0)
      Lambda[j,] = beta * (exp(t[j+1]*alpha) - exp(t[j]*alpha))
    }
  } else if(type == "lp"){
    for(j in 1:(length(t)-1)){
      Lambda0[j] = beta0 * (exp(t[j+1]*alpha0) - exp(t[j]*alpha0))
      Lambda[j,] = beta * (t[j+1]^alpha - t[j]^alpha)
    }
  } else if(type == "ll"){
    for(j in 1:(length(t)-1)){
      Lambda0[j] = beta0 * (exp(t[j+1]*alpha0) - exp(t[j]*alpha0))
      Lambda[j,] = beta * (exp(t[j+1]*alpha) - exp(t[j]*alpha))
    }
  }
  return(list("Lambda0" = Lambda0,"Lambda" = Lambda, "Time" = t[-1]))
}
# 2.2 Cumulative sum ====
Lambda_cum = function(par, type){
  Lambda_re <- Lambda_fun(type = type, par0 = par[-length(par)],m=m)
  delta_Lambda0_hat <- Lambda_re$Lambda0
  delta_Lambda_hat <- Lambda_re$Lambda
  gamma_hat = par[length(par)] 
  Lambda0_hat = cumsum(delta_Lambda0_hat)
  Lambda_hat = apply(delta_Lambda_hat,2,cumsum)
  par_hat = list(Lambda0_hat,Lambda_hat,gamma_hat)
  return(par_hat) 
}
# 2.3 Partial derivative  ====
Lambda_fun_der = function(type, par0, time){
  alpha0 = par0[1]; beta0 = par0[2]
  alpha = par0[2+1:K]; beta = par0[2+K+1:K]
  der_alpha = der_beta = matrix(NA,K,m)
  der_alpha0 = der_beta0 = numeric()
  if(type == "pp"){
    for(j in 1:m){
      der_alpha0[j] = beta0 * (time[j+1]^(alpha0) * log(time[j+1]) - time[j]^(alpha0) * log(time[j]))
      der_beta0[j] = time[j+1]^alpha0 - time[j]^alpha0
      for(k in 1:K){ # beta * (tij^alpha - tij-1^alpha)
        der_alpha[k,j] = beta[k] * (time[j+1]^(alpha[k]) * log(time[j+1]) - time[j]^(alpha[k]) * log(time[j]))
        der_beta[k,j] = time[j+1]^alpha[k] - time[j]^alpha[k]
      }
    }
  } else if(type == "pl"){
    for(j in 1:m){
      # power
      der_alpha0[j] = beta0 * (time[j+1]^(alpha0) * log(time[j+1]) - time[j]^(alpha0) * log(time[j]))
      der_beta0[j] = time[j+1]^alpha0 - time[j]^alpha0
      for(k in 1:K){
        # log-linear
        der_alpha[k,j] = beta[k] * (time[j+1] * exp(alpha[k]*time[j+1]) - time[j] * exp(alpha[k]*time[j]))
        der_beta[k,j] =  exp(alpha[k]*time[j+1]) - exp(alpha[k]*time[j])
      }
    }
  }else if(type == "lp"){
    for(j in 1:m){
      # log-linear
      der_alpha0[j] = beta0 * (time[j+1] * exp(alpha0*time[j+1]) - time[j] * exp(alpha0*time[j]))
      der_beta0[j] = exp(alpha0 * time[j+1]) - exp(alpha0 * time[j])
      for(k in 1:K){
        # power
        der_alpha[k,j] = beta[k] * (time[j+1]^(alpha[k]) * log(time[j+1]) - time[j]^(alpha[k]) * log(time[j]))
        der_beta[k,j] = time[j+1]^alpha[k] - time[j]^alpha[k]
      }
    }
  }else if(type == "ll"){
    for(j in 1:m){
      der_alpha0[j] = beta0 * (time[j+1] * exp(alpha0*time[j+1]) - time[j] * exp(alpha0*time[j]))
      der_beta0[j] = exp(alpha0 * time[j+1]) - exp(alpha0 * time[j])
      for(k in 1:K){
        der_alpha[k,j] = beta[k] * (time[j+1] * exp(alpha[k]*time[j+1]) - time[j] * exp(alpha[k]*time[j]))
        der_beta[k,j] =  exp(alpha[k]*time[j+1]) - exp(alpha[k]*time[j])
      }
    }
  }
  return(list("der_alpha" = der_alpha, "der_beta" = der_beta,
              "der_alpha0" = der_alpha0, "der_beta0" = der_beta0))
}
# 3. Generate simulated data ====
sim.dat.path = function(type, m, K, n, t, par){
  # Denote Parameters
  alpha0 = par[1];beta0 = par[2]; 
  alpha = par[2+1:K]; beta = par[2+K+1:K]
  gamma = par[length(par)]
  # Compute lambda
  Lambda_re = Lambda_fun(type = type, par0 = c(alpha0,beta0,alpha,beta), m = m)
  Lambda0 = Lambda_re$Lambda0; Lambda = Lambda_re$Lambda
  # Generate simulated data
  diff_X_t = diff_Y_t = diff_Z_t = list()
  X_t = Y_t = Z_t = list()
  for(i in 1:n){
    diff_X_t[[i]] = diff_Y_t[[i]] = matrix(NA,m,K) 
    diff_Z_t[[i]] = numeric()
    for(j in 1:m){
      diff_Z_t[[i]][j] = rrIG(1, Lambda0[j], gamma) # Z
      if(is.na(diff_Z_t[[i]][j])) diff_Z_t[[i]][j] = 0
      for(k in 1:K){
        diff_X_t[[i]][j,k] = rrIG(1, Lambda[j,k], gamma) # X
        diff_Y_t[[i]][j,k] = diff_X_t[[i]][j,k] + diff_Z_t[[i]][j] # Y
      }
    }
    diff_Y_t[[i]] = data.frame(diff_Y_t[[i]])
    colnames(diff_Y_t[[i]]) = paste("PC",1:K,sep='')
    Y_t[[i]] = apply(diff_Y_t[[i]],2,cumsum) 
    Y_t[[i]] = data.frame(Y_t[[i]])
    colnames(Y_t[[i]]) = paste("PC",1:K,sep='')
  }
  return(list("Y_t" = Y_t, "diff_Y_t" = diff_Y_t,
              "diff_X_t" = diff_X_t, "diff_Z_t" = diff_Z_t,
              "Lambda" = Lambda, "Lambda0" = Lambda0))
}
# 4. Relevant graph plotting ====
# 4.1 Plot degenerate Data Paths ====
degradation.path.plot = function(data, leg.pos = "none",ech = 5, scale1 = 'fixed'){
  # Adding initial values
  sim_dat_plot <- lapply(data, function(df) {
    new_row <- rep(0,K) # add zero
    rbind(new_row, df)
  })
  data1 = map(sim_dat_plot, ~ mutate(.x, Time = 0:(n()-1)))
  merged_df <- bind_rows(data1, .id = "Unit")
  p1 = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value") %>% 
    ggplot(aes(Time,Value,color = Unit)) + 
    geom_line(alpha=0.8) + geom_point(size=0.6) +
    facet_wrap(vars(PC),nrow = 1, scales = scale1) + 
    theme_bw() +
    scale_color_viridis(discrete = T) + 
    scale_x_continuous(breaks = seq(0, m, by = ech),limits = c(0, m)) + 
    theme(legend.position = leg.pos) #panel.grid = element_blank()
  return(p1)
}
# 4.2 Plot EM iteration ========
EM_iter_plot = function(em_par, par, ncol=4, ture_value = TRUE, 
                        f_names = list(expression(hat(alpha)[0]), expression(hat(alpha)[1]), 
                                       expression(hat(alpha)[2]), expression(hat(alpha)[3]),
                                       expression(hat(beta)[1]), expression(hat(beta)[2]), 
                                       expression(hat(beta)[3]), expression(hat(gamma)))){
  # em_para$par_re[1:em_para$iter,] -> Iteration results of EM algorithm.
  orders = c(paste("alpha",0:K,sep=""),paste("beta",1:K,sep=""),"gamma")
  # Adding mathematical formulas
  f_labeller <- function(variable, value){return(f_names[value])}
  
  d1 = em_par %>% data.frame() %>% 
    mutate("index" = 1:dim(.)[1]) %>% 
    pivot_longer(cols= !index, names_to = "para", values_to = "value") 
  d1$para = factor(d1$para, levels = orders, ordered = TRUE)
  
  if(ture_value == TRUE){
    dummy <- data.frame("para" = colnames(em_par), Z = par)
    dummy$para = factor(dummy$para, levels = orders, ordered = TRUE)
    p1 = ggplot(d1,aes(index,value,color = para)) + geom_line() + 
      facet_wrap(vars(para), ncol = ncol, scales = "free",labeller = f_labeller) +
      scale_color_aaas(name = "Parameters") + 
      geom_hline(data = dummy, aes(yintercept = Z),linetype = "dashed") +
      theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none') +
      xlab("iteration")
  } else{
    p1 = ggplot(d1,aes(index,value,color = para)) + geom_line() + 
      facet_wrap(vars(para), ncol = ncol, scales = "free",labeller = f_labeller) +
      scale_color_aaas(name = "Parameters") + 
      theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none') + 
      xlab("iteration")
  }
  return(p1)
}

# 4.3 Plot path fitting =======
mean.path.fit.plot = function(data, mean_y, leg.pos = "none",
                              ech = 5, ci = TRUE){
  # Data processing
  ## Add the first row of data to be 0
  sim_dat_plot <- lapply(data, function(df) {
    new_row <- rep(0,K)
    rbind(new_row, df)
  })
  Yhat_dat_plot <- lapply(mean_y, function(df) {
    new_row <- rep(0,K)
    rbind(new_row, df)
  })
  # Merge two dataset
  data1 = map(sim_dat_plot, ~ mutate(.x, Time = 0:(n()-1)))
  data_mean = map(Yhat_dat_plot, ~ mutate(.x, Time = 0:(n()-1)))
  merged_df1 <- bind_rows(data1, .id = "Unit")
  merged_df2 <- bind_rows(data_mean, .id = "Unit") # Conflicts with the above "Unit"
  if(ci == TRUE){
    merged_df2$Unit = rep(c("Low","Mean","Up"),each = length(0:m))
  }else{
    merged_df2$Unit = rep("Mean",each = length(0:m))
  }
  colnames(merged_df2)= colnames(merged_df1)
  merged_df = rbind(merged_df1,merged_df2)
  # Plot
  mer_dat = merged_df %>% pivot_longer(cols = !c(Time,Unit), names_to = "PC", values_to = "Value") 
  p1 = mer_dat %>% ggplot(aes(Time,Value,color = factor(Unit), linetype= factor(Unit), size= factor(Unit))) + 
    geom_line(alpha=0.8) + 
    facet_wrap(vars(PC),nrow = 1) + 
    theme_bw() +
    scale_color_manual(name= "", values = c(rep("gray60",n),"#21908C","#440154","#21908C"))+ #"#21908C","#440154","#21908C""blue","red","blue"
    scale_linetype_manual(values = c(rep(1,n),2,1,2)) +
    scale_size_manual(values = c(rep(0.5,n),0.8,0.9,0.8)*0.8) +
    scale_x_continuous(breaks = seq(0, m, by = ech), limits = c(0, m)) +
    scale_y_reverse() +
    theme(legend.position = 'none') #panel.grid = element_blank()
  return(p1)
}
