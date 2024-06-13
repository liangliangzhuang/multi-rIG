#=============== EM algorithm =============
# 1. Initial value estimation ====
init_para_guess = function(m, K, time=0:m, type = "pl", data1, data2, init_param = c(1,1,1,1,3,4,5)){
  unit_mean = unit_mean_diff = unit_var_diff = matrix(NA,m,K)
  # utility funciton
  list_select = function(j = 1, k = 1, data_list = data1){
    re = lapply(data_list, function(df) df[j,k])
    return(unlist(re))
  }
  
  for(k in 1:K){
    for(j in 1:m){
      unit_mean_diff[j,k] = mean(list_select(j = j, k = k, data_list = data1)) # Calculate the mean value of $\delta Y$ at each time $t_j$
      unit_mean[j,k] = mean(list_select(j = j, k = k, data_list = data2)) # Calculate the mean value of $Y$ at each time $t_j$
      unit_var_diff[j,k] = var(list_select(j = j, k = k, data_list = data1)) * n /(n-1) # Calculate the variance value of $Y$
    }
  }
  lambda0_ini_hat = numeric()
  lambda_ini_hat = matrix(NA,m,K)
  gamma_ini_hat = mean(sqrt(unit_mean_diff / unit_var_diff))
  
  # Estimate alpha and beta =====
  # Define the loss function
  loss_function <- function(params, type1 = type, gamma1 = gamma_ini_hat, time = time, unit_mean) {
    # params: alpha 0:K, beta 0:K
    a0 <- params[1];
    ak <- params[1+1:K]; bk <- params[1+K+1:K]
    re = re2 = uu = matrix(NA,m,K)
    if(type1 == "pp") {
      for(k in 1:K){
        uu[,k] = (time[-1]^a0 + bk[k] * time[-1]^ak[k])
        re[,k] = gamma1^3 * unit_mean[,k]^2 / uu[,k] + gamma1 * uu[,k]  
      }
    } else if(type1 == "pl") {
      for(k in 1:K) {
        re[,k] = ((time[-1]^a0 + bk[k] * (exp(time[-1]*ak[k])-1))/gamma1 - unit_mean[,k])^2 
      }
    } else if(type1 == "lp") {
      for(k in 1:K) {
        re[,k] = (((exp(time[-1]*a0)-1) + bk[k] * time[-1]^ak[k])/gamma1 - unit_mean[,k])^2 
      }
    } else if(type1 == "ll") {
      for(k in 1:K) {
        re[,k] = (((exp(time[-1]*a0)-1) + bk[k] * (exp(time[-1]*ak[k])-1))/gamma1 - unit_mean[,k])^2 
      }
    }
    return(sum(re)) 
  }
  
  result0 <- optim(init_param, loss_function, type1=type, unit_mean = unit_mean,
                   time = time, method = "L-BFGS-B",control = list(maxit = 3000),
                   lower = rep(0,2*K+1), upper = c(rep(5,K+1),rep(2000,K)))
  
  para_ini = c(result0$par,gamma_ini_hat) 
  return(para_ini)
}

# 2. Calculate Expectation ====
E_z <- function(method = 5, par1, data, n1 = 10) {
  # Compute the posterior distribution of latent variables p(Z|Y)
  # par1 = c(alpha0, beta0, alpha, beta, gamma)
  # n1 = 10 means order (l) in GL algorithm 
  n = length(data)
  m <- dim(data[[1]])[1]
  K <- dim(data[[1]])[2]
  # GL algorithm 
  rules <- legendre.quadrature.rules(n1)[[n1]] # nodes and weights
  weight = rules$w # weights
  p_fz = function(z,i){
    exp( 
      -1/2 * (apply(par1$delta_Lambda^2 * (data[[i]] - z)^-1 + par1$gamma^2 * (data[[i]] - z),1,sum) + 
                (par1$delta_Lambda0^2 * z^-1 + par1$gamma^2 * z)) +
        -3/2 * log(z) - apply(3/2 * log(data[[i]] - z),1,sum))
  }
  GH_appr = function(rules, data){
    weight = rules$w
    der_pr = map(1:n, ~ (rules$x + 1) %o% apply(data[[.]], 1, min) / 2)
    re <- map(1:n, ~ {
      uu <- sapply(1:n1, function(p) p_fz(der_pr[[.x]][p,], .x))
      rr <- uu %*% diag(weight)
      apply(data[[.x]], 1, min) / 2 * rowSums(rr)
    })
    return(do.call(rbind, re)) # 导出所有n的结果，n*n1
  }
  Pro_z.y <- function(z, data) {
    p_Y.re <- GH_appr(rules = rules, data = data) #n*n1
    p_z.y <- map(1:n, ~ p_fz(z[,.x], i = .x) / p_Y.re[.x,])
    p_z.y <- do.call(rbind, p_z.y) # Compute p(Z|Y)  #n*n1
    return(list("p_z.y" = p_z.y, "p_Y.re" = p_Y.re)) 
  }
  der_pr = map(1:n, ~ (rules$x + 1) %o% apply(data[[.]], 1, min) / 2) # list
  
  uu = rr = list()
  if (method == 1) {
    # log(y-z) 
    re_z1 = function(z, k = k, data = data) Pro_z.y(z, data=data)$p_z.y * t(log((sapply(data, function(x) x[, k]) - z))) 
    E_z.re <- map(1:n, ~ matrix(NA, m, K))
    # Approximation
    for (k1 in 1:K) {
      uu <- lapply(1:n1, function(p) re_z1(sapply(der_pr, function(x) x[p, ]), k = k1, data = data)) #m*n1 
      rr <- lapply(1:n, function(i) sapply(uu, function(x) x[i, ]) %*% diag(weight)) #i-th unit
      for(i in 1:n){
        E_z.re[[i]][, k1] = apply(data[[i]],1,min) /2 * apply(rr[[i]],1,sum)
      }
    }
  } else if (method == 2) {
    # (y - z)^-1 
    p_Y.re <- GH_appr(rules = rules,data = data)
    re_z2 <- function(z, i = i, k = k, data = data) exp( 
      -1/2 * (apply(par1$delta_Lambda^2 * (data[[i]] - z[,i])^-1 + par1$gamma^2 * (data[[i]] - z[,i]),1,sum) + 
                (par1$delta_Lambda0^2 * z[,i]^-1 + par1$gamma^2 * z[,i])) +
        -3/2 * log(z[,i]) - apply(3/2 * log(data[[i]] - z[,i]),1,sum) -log(data[[i]][, k] - z[,i]))  / p_Y.re[i,] 
    
    # Approximation
    E_z.re <- map(1:n, ~ matrix(NA, m, K))
    for (k1 in 1:K) {
      uu <- lapply(1:n1, function(p) {
        t(sapply(1:n, function(i) re_z2(sapply(der_pr, function(x) x[p, ]), i = i, k = k1, data = data)))
      })
      rr <- lapply(1:n, function(i) sapply(uu, function(x) x[i, ]) %*% diag(weight))
      for(i in 1:n){
        E_z.re[[i]][, k1] = apply(data[[i]],1,min) /2 * apply(rr[[i]],1,sum)
      }
    }
    
  } else if (method == 3) {
    # z 
    p_Y.re <- GH_appr(rules = rules, data = data)
    E_z.re <- list()
    re_z3 <- function(z,  data = data){t(sapply(1:n, function(i) p_fz(z[, i], i = i) * z[, i] / p_Y.re[i, ]))}
    # Approximation
    uu <- lapply(1:n1, function(p) re_z3(sapply(der_pr, function(x) x[p, ]), data = data))
    rr <- lapply(1:n, function(i) sapply(uu, function(x) x[i, ]) %*% diag(weight))
    E_z.re <- lapply(1:n, function(i) {apply(data[[i]], 1, min) / 2 * apply(rr[[i]], 1, sum)})
  } else if(method == 4){
    # log(z) 
    re_z4 = function(z, data = data) Pro_z.y(z, data=data)$p_z.y  * t(log(z))
    E_z.re <- list()
    # Approximation
    uu <- lapply(1:n1, function(p) re_z4(sapply(der_pr, function(x) x[p, ]), data = data))
    rr <- lapply(1:n, function(i) sapply(uu, function(x) x[i, ]) %*% diag(weight))
    E_z.re <- lapply(1:n, function(i) {apply(data[[i]], 1, min) / 2 * apply(rr[[i]], 1, sum)})
  }else if(method == 5){
    # z^-1
    p_Y.re <- GH_appr(rules = rules, data = data)
    E_z.re <- list()
    re_z5 <- function(z, i = i, data = data){t(sapply(1:n, function(i) p_fz(z[, i], i = i) * exp(-log(z[, i])) / p_Y.re[i, ])) }
    # Approximation
    uu <- lapply(1:n1, function(p) re_z5(sapply(der_pr, function(x) x[p, ]), data = data))
    rr <- lapply(1:n, function(i) sapply(uu, function(x) x[i, ]) %*% diag(weight))
    E_z.re <- lapply(1:n, function(i) {apply(data[[i]], 1, min) / 2 * apply(rr[[i]], 1, sum)})
  }
  return(E_z.re)
}

# 3. Main EM algorithm ====
EM = function(data, par0, types, tol1 = 0.001, max_iter, time, loglik = F){
  pb <- txtProgressBar(min = 1, max = max_iter, style = 3) # View Progress: Specify minimum, maximum, and style
  n <- length(data) # Number of units
  m <- dim(data[[1]])[1] # Number of measurements
  K <- dim(data[[1]])[2] # Number of PCs
  iter <- 1
  # Set initial parameters
  par_re <- matrix(NA, max_iter, length(par0))
  par_re[1, ] <- par0
  
  while (iter < max_iter) {
    # Iterations
    iter <- iter + 1
    # parameter settings
    init_para = par_re[iter - 1,]
    Lambda_re <- Lambda_fun(type = types, par0 = init_para[-length(init_para)],m=m)
    par <- list("delta_Lambda" = Lambda_re$Lambda, "delta_Lambda0" = Lambda_re$Lambda0, "gamma" = init_para[length(par0)])
    # E-step: Calculate the posterior distribution of the latent variable p(Z|Y) + corresponding expectation ============
    E2 = E_z(method = 2, par1 = par, data = data, n1 = 10) # Equation S10 (1)
    E3 = E_z(method = 3, par1 = par, data = data, n1 = 10) # Equation S10 (2)
    E5 = E_z(method = 5, par1 = par, data = data, n1 = 10) # Equation S10 (3)
    
    # M-step: ===============
    # gamma 
    der_gamma = function(data, raw_par, E = E3){
      gamma_hat = raw_par[length(raw_par)]
      Lambda_re <- Lambda_fun(type = types, par0 = raw_par[-length(raw_par)],m=m)
      delta_Lambda0_hat <- Lambda_re$Lambda0
      delta_Lambda_hat <- Lambda_re$Lambda
      par1 <- list("delta_Lambda" = delta_Lambda_hat, "delta_Lambda0" = delta_Lambda0_hat, "gamma" = gamma_hat)
      E_z.re = matrix(NA, n, m)
      part1 = part2 = numeric()
      for(i in 1:n){
        part1[i] = sum(par1$delta_Lambda0 + apply(par1$delta_Lambda,1,sum))
        part2[i] = sum(-(K-1) * E[[i]] + apply(data[[i]][,],1,sum))
      }
      re = sum(part1, na.rm = T)/sum(part2, na.rm = T) 
      return(re)
    }
    # alpha0 
    alpha0_prime <- function(x, raw_par, data, E = E5) {
      raw_par[1] = x; gamma_hat = raw_par[length(raw_par)]
      # Recalculate Lambda functions
      Lambda_re <- Lambda_fun(type = types, par0 = raw_par,m=m)
      Lambda_hat <- Lambda_re$Lambda
      Lambda0_hat <- Lambda_re$Lambda0 
      part2_re = ree = numeric()
      for(i in 1:n){
        lam_der = Lambda_fun_der(type = types, par0 = raw_par, time = 0:m)$der_alpha0
        part2_re = lam_der * (1/Lambda0_hat + gamma_hat - Lambda0_hat * E5[[i]])
        ree[i] = sum(part2_re, na.rm = T)
      }
      return(sum(ree, na.rm = T))
    }
    # alphak
    alphak_prime <- function(x, raw_par, data, k=1, E = E2) {
      raw_par[2+k] = x; gamma_hat = raw_par[length(raw_par)]
      # Recalculate Lambda functions
      Lambda_re <- Lambda_fun(type = types, par0 = raw_par,m=m)
      Lambda_hat <- Lambda_re$Lambda
      Lambda0_hat <- Lambda_re$Lambda0 
      part2_re = ree = numeric()
      for(i in 1:n){
        lam_der = Lambda_fun_der(type = types, par0 = raw_par, time = 0:m)$der_alpha[k,]
        part2_re = lam_der * (1/Lambda_hat[,k] + gamma_hat - Lambda_hat[,k] * E[[i]][,k])
        ree[i] = sum(part2_re, na.rm = T)
      }
      return(sum(ree, na.rm = T))
    }
    # betak 
    betak_prime <- function(x, raw_par, data, k=1, E = E2) {
      raw_par[2+K+k] = x; gamma_hat = raw_par[length(raw_par)]
      # Recalculate Lambda functions
      Lambda_re <- Lambda_fun(type = types, par0 = raw_par,m=m)
      Lambda_hat <- Lambda_re$Lambda
      Lambda0_hat <- Lambda_re$Lambda0 
      part2_re = ree = numeric()
      for(i in 1:n){
        lam_der = Lambda_fun_der(type = types, par0 = raw_par, time = 0:m)$der_beta[k,]
        part2_re = lam_der * (1/Lambda_hat[,k] + gamma_hat - Lambda_hat[,k] * E[[i]][,k])
        ree[i] = sum(part2_re, na.rm = T)
      }
      return(sum(ree, na.rm = T))
    }
    
    init_para_new = init_para  # Replace estimated parameters
    
    # Derivative of gamma parameter ====
    hat_gamma = der_gamma(data = data, raw_par = init_para_new)
    init_para_new[length(init_para_new)] = hat_gamma 
    
    # Derivative of alpha0 parameter ====
    if(types == "pp" | types == "pl"){
      hat_alpha0 = tryCatch(uniroot(alpha0_prime, raw_par = init_para_new, data = data, c(0.001,1.4))$root,
                            error = function(e) return(NA))
    } else if(types == "lp" | types == "ll"){
      hat_alpha0 = tryCatch(uniroot(alpha0_prime, raw_par = init_para_new, data = data, c(0.001,2)*0.2)$root,
                            error = function(e) return(NA))
    }
    init_para_new[1] = hat_alpha0 
    # Derivative of beta0 parameter (Since the model is not identifiable, let beta0 be 1) ====
    hat_beta0 = 1 
    init_para_new[2] = hat_beta0
    # Derivative of alphak parameter ====
    hat_alphak = numeric()
    for(k in 1:K){
      if(types == "pp" | types == "lp"){
        hat_alphak[k] = tryCatch(uniroot(alphak_prime, k=k, raw_par = init_para_new, data = data, c(0.001,5))$root,
                                 error = function(e) return(NA))
      } else if(types == "ll" | types == "pl"){
        hat_alphak[k] = tryCatch(uniroot(alphak_prime, k=k, raw_par = init_para_new, data = data, c(0.001,5)*0.1)$root,
                                 error = function(e) return(NA))
      }
      init_para_new[2+k] = hat_alphak[k] 
    }
    # Derivative of betak parameter ====
    hat_betak = numeric()
    for(k in 1:K){
      hat_betak[k] = tryCatch(uniroot(betak_prime, k=k, raw_par = init_para_new, data = data, c(0.01,100))$root,
                              error = function(e) return(NA))
      init_para_new[2+K+k] = hat_betak[k] 
    }
    # =======================
    setTxtProgressBar(pb, iter) # Display progress bar
    par_re[iter, ] = init_para_new
    # Stopping criteria
    if(sqrt(sum((par_re[iter, ] - par_re[iter-1, ])^2,na.rm=T)) < tol1) break
  }
  
  colnames(par_re) = c("alpha0","beta0",paste0("alpha",1:K, sep = ""), paste0("beta",1:K, sep = ""),"gamma")
  close(pb) # Close progress bar
  
  # Calculate the log-likelihood function
  if(loglik == TRUE){
    Lambda_re_final <- Lambda_fun(type = types, par0 = init_para_new[-length(init_para_new)],m=m)
    Lambda0_final_hat <- Lambda_re_final$Lambda0
    Lambda_final_hat <- Lambda_re_final$Lambda
    gamma_final_hat = init_para_new[length(init_para_new)]
    par_final <- list("delta_Lambda" = Lambda_final_hat, "delta_Lambda0" = Lambda0_final_hat, "gamma" = gamma_final_hat)
    E1 = E_z(method = 1, par1 = par_final, data = data, n1 = 10)
    E4 = E_z(method = 4, par1 = par_final, data = data, n1 = 10)
    part1_re = matrix(NA,m,K); part2_re = ree = numeric()
    for(i in 1:n){
      for(k in 1:K){
        part1_re[,k] = log(Lambda_final_hat[,k]) + gamma_final_hat * Lambda_final_hat[,k] - 3/2 * E1[[i]][,k] - (Lambda_final_hat[,k]^2)/2 * E2[[i]][,k] - gamma_final_hat^2/2 * (data[[i]][,k] - E3[[i]]) - log(sqrt(2*pi))
      }
      part2_re = log(Lambda0_final_hat) + gamma_final_hat * Lambda0_final_hat - 3/2 * E4[[i]] - (Lambda0_final_hat^2)/2 * E5[[i]] - gamma_final_hat^2/2 * E3[[i]] - log(sqrt(2*pi))
      ree[i] = sum(part1_re + part2_re, na.rm = T)
    }
    log_lik = sum(ree)
  }else{log_lik = NULL}
  
  return(list("par_re" = par_re, 
              "iter" = iter, 
              "mse_all" = sqrt(sum((par_re[iter, ] - par_re[iter-1, ])^2,na.rm=T)),
              "log_lik" = log_lik))
}


