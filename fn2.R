# This .R includes the code of MCEM to approximate the mle in the correlated random individual effect model (D3).
# This approximated mle will be used as the initialization of MCMC to reduce the computational cost for convergence.
# Instead of running this .R, you can run run1.R that will use the following functions.

#The MCEM algorithm
EM_SIR = function(Y1, X1, Z1, Y2, X2, Z2, Y3, X3, Z3, 
                  rho, sigma_u1, sigma_u2, param1, sigma_g, param2,
                  noI, noI2, size, df, max_iter) { 
  # first initialize the parameter at (approximate) mle
  # at each iteration:
  #   update joint random effects
  #   update reproduction random effects
  #   update growth fixed effect parameters
  #   update reproduction fixed effect parameters
  #   update growth variance parameter
  #   update random effect distribution, including variance and correlation parameters
  param_store      = matrix(0, nrow = max_iter, ncol = 4 + length(param1) + length(param2))
  param_store[1, ] = c(rho, sigma_u1, sigma_u2, param1, sigma_g, param2)
  D = matrix(c(sigma_u1**2, rep(rho*sigma_u1*sigma_u2, 2), sigma_u2**2), nrow = 2)
  sig = c(chol(D)[1,1], chol(D)[1,2], chol(D)[2,2])
  pred1 = X1 %*% param1
  pred2 = X2 %*% param2
  pred3 = X3 %*% param2
  for (iter in 2:max_iter) {
    posterior1 = matrix(0, nrow = size, ncol = noI )
    posterior2 = matrix(0, nrow = size, ncol = noI )
    posterior3 = matrix(0, nrow = size, ncol = noI2)
    D_inv = solve(D)
    u_store = sir_mean12(D_inv, sigma_g, Y1, Z1, Y2, Z2, pred1, pred2, noI)
    v_store = sir_mean3(sigma_u2, Y3, Z3, pred3, noI2)
    for (i in 1:noI) {
      temp = SIR_u12(size, df, u_store[i,], D_inv, D, sigma_g,
                     Y1[Z1[,i]==1], pred1[Z1[,i]==1], Y2[Z2[,i]==1], pred2[Z2[,i]==1])
      posterior1[,i] = temp[,1]
      posterior2[,i] = temp[,2]
    }
    for (i in 1:noI2) posterior3[,i] = SIR_u3(size, v_store[i], sigma_u2, Y3[Z3[,i]==1], pred3[Z3[,i]==1])
    
    lm_Y = Y1 - Z1 %*% colMeans(posterior1)
    param1 = lm(lm_Y ~ X1 - 1)$coefficients
    param2 = optim(param2, logistic, large_Y = c(Y2,Y3), large_X = rbind(X2,X3),
                   Z_u = rbind(Z2 %*% t(posterior2), Z3 %*% t(posterior3)), size = size, method="BFGS")$par
    pred1 = X1 %*% param1
    pred2 = X2 %*% param2
    pred3 = X3 %*% param2
    sigma_g = sqrt(sum((lm_Y - pred1)**2)/length(Y1))
    
    sig = optim(sig, rn_update,posterior1=posterior1, posterior2=posterior2,noI=noI,method="BFGS")$par
    D = matrix(c(sig[1], sig[2],0,sig[3]), nrow = 2) %*% matrix(c(sig[1], 0,sig[2],sig[3]), nrow = 2)
    sigma_u1 = sqrt(D[1,1])
    rho = D[1,2] / sigma_u1 / sigma_u2
    D[2,2] = (sum(colMeans(posterior2**2)) + sum(colMeans(posterior3**2))) / (noI + noI2)
    sigma_u2 = sqrt(D[2,2])
    D[1,2] = rho * sigma_u1 *sigma_u2; D[2,1] = rho * sigma_u1 *sigma_u2
    sig = c(chol(D)[1,1], chol(D)[1,2], chol(D)[2,2])
    param_store[iter, ] = c(rho, sigma_u1, sigma_u2, param1, sigma_g, param2)
  }
  return(param_store)
}

###update joint random effects###
#posterior mean of the joint random effects
sir_mean12 = function(D_inv, sigma_g, Y1, Z1, Y2, Z2, pred1, pred2, noI) { 
  u_store = matrix(0, nrow = noI, ncol = 2)
  for (i in 1:noI) {
    u_store[i,2] = optim(0, mean_fn12, D_inv=D_inv,sigma_g=sigma_g,
                         part_Y1=Y1[Z1[,i]==1],
                         part_Y2=Y2[Z2[,i]==1],
                         part_pred1=pred1[Z1[,i]==1],
                         part_pred2=pred2[Z2[,i]==1],
                         n=sum(Z1[,i]==1),method = "BFGS")$par
    u_store[i,1] = (sum(Y1[Z1[,i]==1] - pred1[Z1[,i]==1])/sigma_g**2 -
                      u_store[i,2] * D_inv[1,2])/
      (D_inv[1,1] + sum(Z1[,i]==1)/sigma_g**2)
  }
  return(u_store)
}

#cost fn for estimating the posterior mean of the joint random effects 
mean_fn12 = function(u2, D_inv, sigma_g, part_Y1, part_pred1, part_Y2, part_pred2, n) {
  u1 =(sum(part_Y1-part_pred1)/sigma_g**2-u2 * D_inv[1,2]) / (D_inv[1,1] + n/sigma_g**2)
  p = 1/(1+exp(-part_pred2-u2))
  (sum(part_Y2 - p) - D_inv[2,2] * u2 - D_inv[1,2] * u1)**2
}

#sample importance resampling for the joint random effects
SIR_u12 = function(size, df, u, D_inv, D, sigma_g, part_Y1, part_pred1, part_Y2, part_pred2) { 
  p = 1/(1+exp(- part_pred2 - u[2]))
  fisher = matrix(c(D_inv[1,1] + length(part_Y1)/sigma_g**2, D_inv[1,2],
                    D_inv[1,2], D_inv[2,2] + sum(p*(1-p))), nrow = 2)
  prop = rmvt(size, solve(fisher), df = df, delta = u)
  l = apply(prop, 1, function(u) likeli_con12(u, part_Y1, part_pred1, part_Y2, part_pred2, sigma_g)) #likelihood
  q = dmvt(prop, u, solve(fisher), df = df, log = TRUE) #proposal likelihood
  l_2 = dmvnorm(prop, c(0,0), D, log = TRUE) #prior likelihood
  
  prob = exp(l+l_2-q-min(l+l_2-q))
  prob = prob / sum(prob)
  idx = sample(1:size, size, prob=prob,replace=TRUE)
  return(prop[idx,])
}

#loglikelihood of the joint random effects
likeli_con12 = function(u, part_Y1, part_pred1, part_Y2, part_pred2, sigma_g) { 
  l1 = sum(dnorm(part_Y1, part_pred1 + u[1], sigma_g, log = TRUE))
  p  = 1/(1+exp(-part_pred2 - u[2]))
  l2 = sum(part_Y2 * log(p) + (1-part_Y2) * log(1-p))
  return(l1 + l2)
}
###update joint random effects###

###update reproduction random effects###
#posterior mean of the reproduction random effects
sir_mean3 = function(sigma_u2, Y3, Z3, pred3, noI2) { 
  v_store = numeric(noI2)
  for (i in 1:noI2) {
    v_store[i] = optim(0,mean_fn3,sigma_u2=sigma_u2,
                       part_Y3=Y3[Z3[,i]==1],part_pred3=pred3[Z3[,i]==1],method = "BFGS")$par
  }
  return(v_store)
}

#cost fn for estimating the posterior mean of the reproduction random effects 
mean_fn3  = function(v, sigma_u2, part_Y3, part_pred3) {
  p = 1/(1+exp(-part_pred3-v))
  (sum(part_Y3 - p) - v / sigma_u2**2)**2
}

#sample importance resampling for the reproduction random effects
SIR_u3 = function(size, v, sigma_u2, part_Y3, part_pred3) { 
  p = 1/(1+exp(- part_pred3 - v))
  sim_sd = 1 / sqrt(1 / sigma_u2**2 + sum(p*(1-p)))
  prop = rnorm(size, v, sim_sd)
  l = unlist(lapply(prop, function(u) likeli_con3(v, part_Y3, part_pred3))) #likelihood
  q = dnorm(prop, v, sim_sd, log = TRUE) #proposal likelihood
  l_2 = dnorm(prop, 0, sigma_u2, log = TRUE) #prior likelihood
  
  prob = exp(l+l_2-q-min(l+l_2-q))
  prob = prob / sum(prob)
  idx = sample(1:size, size, prob=prob,replace=TRUE)
  return(prop[idx])
}

#loglikelihood of the reproduction random effects
likeli_con3  = function(v, part_Y3, part_pred3) { 
  p  = 1/(1+exp(- part_pred3 - v))
  l3 = sum(part_Y3 * log(p) + (1-part_Y3) * log(1-p))
  return(l3)
}
###update reproduction random effects###




###update reproduction fixed parameters###
#loglikelihood of reproduction given random effects
logistic = function(param2, large_Y, large_X, Z_u, size) { 
  pred = large_X %*% param2
  temp = matrix(rep(pred,size), ncol = size)
  p = rowSums(1/(1+exp(- temp - Z_u))) / size
  return(-sum(large_Y * log(p) + (1-large_Y) * log(1-p)))
}
###update reproduction fixed parameters###

###update random effect distribution ###
#loglikelihood of the random effect distribution
rn_update = function(sig, posterior1, posterior2, noI) { 
  D = matrix(c(sig[1],sig[2],0,sig[3]), nrow = 2) %*% matrix(c(sig[1],0,sig[2],sig[3]), nrow = 2)
  D_inv = solve(D)
  l12 = 0
  for (i in 1:noI) {
    l12 = l12 + mean(apply(cbind(posterior1[,i],posterior2[,i]), 1, function(x) 
      x[1]**2 * D_inv[1,1] + 2*x[1]*x[2]*D_inv[1,2] + x[2]**2 * D_inv[2,2] ))
  }
  return(noI*log(det(D)) + l12)
}
###update random effect distribution ###
