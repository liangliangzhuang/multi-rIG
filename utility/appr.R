# ======== CDF approximation of GL approximate failure time ============

# 1. Data augmentation =====
pivot_expand_dat = function(bin = 2, data = sim_dat$Y_t){
  # bin -> scale
  seq1 = seq(0, m, 1/bin)
  Y_new_dat = list()
  for(i in 1:n){
    Y_new_dat[[i]] = matrix(NA,length(seq1)-1,K)
    for(k in 1:K){
      Y_dat = c(0,data[[i]][,k])
      Y_new = Y_dat[1]
      for(j in 1:(m)){
        Y2 = seq(Y_dat[j], Y_dat[j+1], (Y_dat[j+1] - Y_dat[j])/bin)[-1]
        Y_new = c(Y_new,Y2)
      }
      Y_new_dat[[i]][,k] = Y_new[-1]
    }
  }
  return(list("Time" = seq1[-1], "Y_new" = Y_new_dat))
}
pivot_expand_lambda = function(bin = 2, par = em_para_re){
  # bin -> scale
  Lambda_re <- Lambda_fun(type = types, par0 = par[-length(par)],bin = bin,m=m)
  delta_Lambda0_hat <- Lambda_re$Lambda0
  delta_Lambda_hat <- Lambda_re$Lambda
  gamma_hat = par[length(par)] 
  Lambda0_hat = cumsum(delta_Lambda0_hat) 
  Lambda_hat = apply(delta_Lambda_hat,2,cumsum) 
  par_hat = list(Lambda0_hat,Lambda_hat,gamma_hat)
  return("par_hat" = par_hat)
}

# 2. Calculate the distribution of failure time ====
p_fz <- function(z, j, para2, D = D) {
  # Function within integral Eq.(10)
  re <- numeric()
  lambda0 = para2[[1]];  lambda = para2[[2]]; gamma2 = para2[[3]]
  for (k in 1:K) {
    re[k] <- prIG(D[k] - z, lambda1 = lambda[j, k], gamma1 = gamma2)
  }
  ru <- (1 - prod(re)) * drIG(z, lambda1 = lambda0[j], gamma1 = gamma2)
  return(ru)
}

# 3. CDF approximation method (empirical distribution vs GL method) ====
Approx_CDF_failure_time <- function(i = i, n1 = 30, bins = 5, div = 10, mc_n = 1000,
                                    D, par, data2) {
  # 3.1 Data augmentation ====
  par_hat = pivot_expand_lambda(bin = bins, par = par)
  expand_dat = pivot_expand_dat(bin = bins, data = data2)
  time1 = expand_dat$Time; data = expand_dat$Y_new
  # 3.2 GL method ====
  SYS_t = Sys.time() 
  rules <- legendre.quadrature.rules(n1)[[n1]] # obtain weights and nodes provided in the R package {gaussquad}
  Y_min <- apply(data[[i]], 1, min)
  weight <- rules$w
  der_pr <- (rules$x + 1) %o% Y_min / 2 # n1*m
  uu <- rr <- matrix(NA, length(time1), n1)
  for (p in 1:n1) {
    for (j in 1:length(Y_min)) {
      uu[j, p] <- tryCatch(p_fz(der_pr[p, j], j = j, para2 = par_hat, D = D),
                           error = function(e) return(NA))
    }
  }
  for (h in 1:length(time1)) rr[h, ] <- uu[h, ] * weight
  re_gh <- Y_min / 2 * apply(rr, 1, sum, na.rm = T) # 矩阵
  gh_time = Sys.time() - SYS_t
  # 3.3 Empirical approximation method ====
  SYS_t = Sys.time() 
  epoch <- Y_min / div
  der_pr2 <- t(epoch / 2 + epoch %o% ((1:div) - 1))
  ree <- matrix(NA, length(Y_min), div)
  for (p in 1:div) {
    for (j in 1:length(Y_min)) {
      ree[j, p] <- tryCatch(p_fz(z = der_pr2[p, j], j = j, para2 = par_hat, D = D),
                            error = function(e) return(NA)) * epoch[j]
    }
  }
  re_em <- apply(ree, 1, sum, na.rm = T)
  em_time = Sys.time() - SYS_t
  # 3.4 MCMC method ====
  SYS_t = Sys.time() 
  mcm <- matrix(NA, length(Y_min), mc_n)
  for (p in 1:mc_n) {
    for (j in 1:length(Y_min)) {
      x1 <- runif(1,0,Y_min[j])
      mcm[j, p] = tryCatch(p_fz(z = x1, j = j, para2 = par_hat, D = D),error = function(e) return(NA))  * Y_min[j]
    }
  }
  re_mcmc <- apply(mcm,1, mean, na.rm = T)
  mcmc_time = Sys.time() - SYS_t
  # 3.5 Compute approximate results ====
  ft_dat <- data.frame("Time" = time1[-1], "GH" = re_gh[-1], "Empirical" = re_em[-1], "MCMC" = re_mcmc[-1])
  p1 = ggplot(data = ft_dat) + 
    geom_bar(aes(x = Time,y = Empirical), stat = "identity", fill = "gray") + 
    geom_line(aes(x = Time,y = GH)) +
    geom_line(aes(x = Time,y = MCMC)) 
  
  return(list(ft_dat,p1,"cost_time" = c(gh_time,em_time,mcmc_time)))
}
# 4. Plot the integral approximation results under different scenarios ====
scen_appr_CDF = function(types, D, sim_dat, para, limits=c(5, 10), leg.pos = "none"){
  # Results in different order (l)
  arr_re = Approx_CDF_failure_time(i = 1, n1 = 5, bins = 10, div = 30, D=D, par=para, data2 = sim_dat$Y_t)
  arr_re2 = Approx_CDF_failure_time(i = 1, n1 = 7, bins = 10, div = 30, D=D, par=para, data2 = sim_dat$Y_t)
  arr_re3 = Approx_CDF_failure_time(i = 1, n1 = 10, bins = 10, div = 30, D=D, par=para, data2 = sim_dat$Y_t)
  arr_re4 = Approx_CDF_failure_time(i = 1, n1 = 15, bins = 10, div = 30, D=D, par=para, data2 = sim_dat$Y_t)
  # Results summary
  all_appr = data.frame(cbind(arr_re[[1]][,c(1,3)],arr_re[[1]][,2],arr_re2[[1]][,2],arr_re3[[1]][,2],arr_re4[[1]][,2],arr_re[[1]][,4]))
  colnames(all_appr) = c("Time","Empirical","GL-5","GL-7","GL-10","GL-15","MCMC") 
  # Plot
  pp_com = all_appr %>% ggplot() + 
    geom_bar(aes(x = Time, y = Empirical,fill = "Empirical"), stat = "identity", alpha = 0.5) + 
    geom_line(aes(x = Time, y = `GL-5`, color = "GL (l=5)"), linetype=6, size = 1) +
    geom_line(aes(x = Time, y = `GL-7`, color = "GL (l=7)"), linetype=8, size = 1) +
    geom_line(aes(x = Time, y = `GL-10`, color = "GL (l=10)"), linetype=1, size = 1) +
    geom_line(aes(x = Time, y = `GL-15`, color = "GL (l=15)"), linetype=4, size = 1) +
    geom_line(aes(x = Time, y = `MCMC`, color = "MCMC"), linetype=10, size = 1) +
    scale_x_continuous(limits =limits) +
    labs(x = "Time", y = "CDF") +
    scale_color_manual(values = c("Empirical" = "gray60", "GL (l=5)" = "#7AD151", 
                                  "GL (l=7)" = "#FDE724", "GL (l=10)" = "#693476",
                                  "GL (l=15)" = "#48A29F", "MCMC" = "gray20"
    )) +
    scale_fill_manual("",values="gray60")+
    theme_bw() + theme(panel.grid = element_blank(), legend.position = leg.pos,
                       legend.key=element_blank(),
                       legend.title=element_blank()) 
  # Time comparison =====
  pp_time = c(arr_re$cost_time[1],arr_re2$cost_time[1],arr_re3$cost_time[1],arr_re4$cost_time[1],arr_re$cost_time[3])
  
  return(list("plot" = pp_com, "cost_time" = pp_time))
}




