# ========================================================================
# ===== Multivariate inverse Gaussian process with common effects ========
# Tested in R 4.2.0., implemented on 2024.06.12.  
# The latest code updates can be found at https://github.com/liangliangzhuang/multi-rIG 
# The fatigue crack-size data are openly available in Meeker, W.Q., et al., 2022. Statistical Methods for Reliability Data. John Wiley & Sons. 
# ========================================================================
# load packages
library(SuppDists)
library(gaussquad)
library(tidyverse)
library(viridis)
library(ggsci)
library(latex2exp)
library(cowplot)
library(openxlsx)
library(ggplot2)
# load functions
source("utility/fct.R")
source("utility/em.R")

#===============================================================
#======= Data Analysis for the Fatigue Crack-size Data =========
#===============================================================
# 1. Load the data from .xlsx file =======
fc_dat = fc_dat1 = fc_diff_dat = list()
n = 6 # Number of units
for (i in 1:n) {
  fc_dat[[i]] = read.xlsx("case/Fatigue-crack-size.xlsx",sheet = i) 
  fc_dat1[[i]] =  (fc_dat[[i]][-1,] - 0.9)*10 
  fc_diff_dat[[i]] = apply(fc_dat[[i]],2,diff) *10 
}

m = dim(fc_diff_dat[[1]])[1] # Number of measurements
K = dim(fc_diff_dat[[1]])[2] # Number of PCs
thr = c(0.315,0.245,0.14) *10 # Degradation Threshold

# Figure 10, plot of the data
degradation.path.plot(data=fc_dat1, leg.pos = "right", ech = 2) +
  xlab(TeX(r'(Millions of cycles)')) +
  ylab(TeX(r'(Crack Length (inches$\times 10$))') ) 
# ggsave(paste("case/result/crack_dat.pdf",sep=''), width = 7, height = 3)

# 2. Model fitting and selection ========
# 2.1 Four scenarios ========
# Scen.I -> "pp", Scen.II -> "pl", Scen.III -> "lp", Scen.IV -> "ll"
max_iter = 1000; tol1 = 10^-4
type3 = c("pp","pl","lp","ll") 
em_para_re2 = para_ini2 = matrix(NA,4,9)
loglik = aic = numeric()
em_plot = em_iter = Yhat_dat = list()
ttt = Sys.time()
for(k in 1:length(type3)){
  # 2.1.1 Determine initial parameter estimation values (see Section 4.2)
  if(type3[k] == "pp"){co1 = 1;co2 = 1}else if(type3[k] == "pl"){co1 = 1;co2 = 0.1}else if(type3[k] == "lp"){co1 = 0.1;co2 = 1}else if(type3[k] == "ll"){co1 = 0.1;co2 = 0.1}
  para_ini3 = init_para_guess(m = m, K = K, time = 0:m, type = type3[k], data1 = fc_diff_dat, data2 = fc_dat1, 
                              init_param = c(1*co1,1*co2,1*co2,1*co2,1*1,1*1,1*1))
  para_ini2[k,] = c(para_ini3[1],1,para_ini3[2:length(para_ini3)])
  # 2.1.2 EM algorithm for point estimation (see Section 4.1)
  em_para2 = EM(data = fc_diff_dat, par0 = para_ini2[k,], types = type3[k], tol1 = tol1, max_iter = max_iter, time = 0:m, loglik = TRUE)
  em_iter[[k]] = em_para2$par_re[1:em_para2$iter,-2]
  # 2.1.3 Trace plots
  em_plot[[k]] = EM_iter_plot(em_par = em_iter[[k]], 
                              par = rep(NA,2*K+2), ncol=4, ture_value = F)
  em_para_re2[k,] = em_para2$par_re[em_para2$iter,]
  loglik[k] = em_para2$log_lik
  aic[k] = -2*em_para2$log_lik + 2*(length(em_para_re2) - 1) # -1 is to remove beta0
  # Calculate the mean degradation value of the estimated results
  Yhat_dat[[k]] = matrix(NA,m,K)
  lam_dat = Lambda_cum(par = em_para_re2[k,], type = type3[k])
  for(j in 1:K){Yhat_dat[[k]][,j] = (lam_dat[[1]] + lam_dat[[2]][,j])/lam_dat[[3]]}
  Yhat_dat[[k]] = data.frame(Yhat_dat[[k]])
  colnames(Yhat_dat[[k]]) = paste("PC",1:K,sep='')
}

# EM iteration —— Figure S2
cowplot::plot_grid(plotlist = em_plot, nrow = 4,labels = c("I","II","III","IV")) + xlab("")
# ggsave(paste("case/result/em-iteration.pdf",sep=''), width = 9, height = 12)

# Estimated mean degradation path —— Figure 11 (b)
mean_path_fit =list()
for(k in 1:4){
  mean_path_fit[[k]] = mean.path.fit.plot(data=fc_dat1, mean_y = list(Yhat_dat[[k]]), leg.pos = "none", ech = 5, ci = F)
}
cowplot::plot_grid(plotlist = mean_path_fit)
# ggsave(paste("case/result/em-all-crack.pdf",sep=''), width = 9, height = 12)

# Table 4 =======
tabel4 = round(cbind(em_para_re2,aic),3)  
tabel4 = data.frame(tabel4[,-2])
colnames(tabel4) = c("alpha0", paste("alpha",1:K,sep=''),paste("beta",1:K,sep=''),"gamma","AIC")
tabel4
# write.csv(tabel4, paste("case/result/Para_summary.csv",sep=''))
# save.image(file = paste("case/result/final-crack.RData", sep=''))






