# This .R includes the nimble code and required libraries for estimating vital rate parameters.
# Instead of running this .R, you can run run1.R that will use the following functions.
library(MASS)
library(nimble)
library(lme4)
library(parallel)
library(mvtnorm)
indString <- nimbleCode({ # Model I1: vanilla model.
  pred1[1:N1,1] <- X1[1:N1,1:2] %*% param1[1:2]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2]
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) #pred1[i]
  for (i in 1:N2) Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
})

iitString <- nimbleCode({ # Model I2: correlated random year effect model.
  pred1[1:N1,1] <- X1[1:N1,1:2] %*% param1[1:2] + Z1[1:N1,1:noI] %*% u[1:noI,1]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2] + Z2[1:N2,1:noI] %*% u[1:noI,2]
  for (i in 1:N1) {
    Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) 
  }
  for (i in 1:N2) {
    Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  }
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
  sigma_u1~dunif(0,5)
  sigma_u2~dunif(0,5)
  for (i in 1:noI) { #W is a precision matrix
    u[i,1] ~ dnorm(0,sd = sigma_u1)
    u[i,2] ~ dnorm(0,sd = sigma_u2)
  }
})

iiiString <- nimbleCode({ # Model I3: correlated random individual effect model.
  pred1[1:N1,1] <- X1[1:N1,1:2] %*% param1[1:2] + Z1[1:N1,1:noI ] %*% u[1:noI,1]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2] + Z2[1:N2,1:noI ] %*% u[1:noI,2]
  pred3[1:N3,1] <- X3[1:N3,1:2] %*% param2[1:2] + Z3[1:N3,1:noI2] %*% v[1:noI2 ]
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) 
  for (i in 1:N2) Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  for (i in 1:N3) Y3[i] ~ dbinom(1/(1+exp(-pred3[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
  sigma_u1~dunif(0,5)
  sigma_u2~dunif(0,5)
  for (i in 1:noI) { #W is a precision matrix
    u[i,1] ~ dnorm(0,sd = sigma_u1)
    u[i,2] ~ dnorm(0,sd = sigma_u2)
  }
  for (i in 1:noI2) v[i] ~ dnorm(0,sd = sigma_u2)
})

repString <- nimbleCode({ # Model D1a: reproduction conditional model.
  pred1[1:N1,1] <- X1[1:N1,1:3] %*% param1[1:3]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2]
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) #pred1[i]
  for (i in 1:N2) Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param1[3]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
})

# The closed-form likelihood is quite complicated. Instead of running it on nimble, we hand-coded the function,
# which you will find after the nimble code for Model D3.

driString <- nimbleCode({ # Model D2a: shared drivers model.
  pred1[1:N1,1] <- X1[1:N1,1:3] %*% param1[1:3]
  pred2[1:N2,1] <- X2[1:N2,1:3] %*% param2[1:3]
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) #pred1[i]
  for (i in 1:N2) Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param1[3]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  param2[3]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
})

mitString <- nimbleCode({ # Model D2b: correlated random year effect model.
  pred1[1:N1,1] <- X1[1:N1,1:2] %*% param1[1:2] + Z1[1:N1,1:noI] %*% u[1:noI,1]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2] + Z2[1:N2,1:noI] %*% u[1:noI,2]
  for (i in 1:N1) {
    Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) 
  }
  for (i in 1:N2) {
    Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  }
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
  W[1:2,1:2] ~ dwish(Omega[1:2,1:2],3)
  for (i in 1:noI) { #W is a precision matrix
    u[i,1:2] ~ dmnorm(mu[1:2],W[1:2,1:2])
  }
  W_inv[1:2,1:2] <- inverse(W[1:2,1:2])
  sigma_u1 <- sqrt(W_inv[1,1])
  sigma_u2 <- sqrt(W_inv[2,2])
  rho <- W_inv[1,2] / sigma_u1 / sigma_u2
})

miiString <- nimbleCode({ # Model D3: correlated random individual effect model.
  pred1[1:N1,1] <- X1[1:N1,1:2] %*% param1[1:2] + Z1[1:N1,1:noI ] %*% u[1:noI,1]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2] + Z2[1:N2,1:noI ] %*% u[1:noI,2]
  pred3[1:N3,1] <- X3[1:N3,1:2] %*% param2[1:2] + Z3[1:N3,1:noI2] %*% v[1:noI2 ]
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) 
  for (i in 1:N2) Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  for (i in 1:N3) Y3[i] ~ dbinom(1/(1+exp(-pred3[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
  W[1:2,1:2] ~ dwish(Omega[1:2,1:2],3)
  for (i in 1:noI ) u[i,1:2] ~ dmnorm(mu[1:2],W[1:2,1:2]) #W is a precision matrix
  for (i in 1:noI2) v[i]     ~ dnorm(0,sd = sigma_u2)
  
  W_inv[1:2,1:2] <- inverse(W[1:2,1:2])
  sigma_u1 <- sqrt(W_inv[1,1])
  sigma_u2 <- sqrt(W_inv[2,2])
  rho <- W_inv[1,2] / sigma_u1 / sigma_u2
})

###Model D1b: copula model###
# note that we use MH with symmetric proposal here
F_3 <- function(y_3, q) ifelse(y_3 > 1, q + (y_3 - 1) * (1 - q), y_3 * q)
f_3 <- function(y_3, q) ifelse(y_3 > 1, 1 - q, q)

den = function(y3,q,y1,sd_g,D) { # density (or likelihood) of copula model 
  sum(dmvnorm(cbind(y1,pnorm(F_3(y3,q))),c(0,0),D,log=TRUE) + log(f_3(y3,q)) - dnorm(pnorm(F_3(y3,q)),log=TRUE)) - 
    length(y1)*log(sd_g) 
}

MCMC_cop = function(Y1, Y2, X1, X2, Y0, X0, propo, inits, burn, iteration) {
  n1 = 2 # no of fixed effect, variance parameter w.r.t. continuous parameter
  n2 = 2 # no of fixed effect parameter w.r.t. binary parameter
  n3 = 2 # no of random effect distribution parameter, i.e. 3
  
  param = matrix(0, nrow = iteration, ncol = n1+n2+n3) # all parameters of interest
  param[1,] = inits
  D = matrix(c(1,param[1,6],param[1,6],1), nrow = 2)
  D = matrix(c(1,0,0,1), nrow = 2)
  u = runif(length(Y1))
  Y3= Y2 + u
  p1= param[1,1:2]
  p2= param[1,3:4]
  sd_g = param[1,5]
  z = (Y1 - X1 %*% p1) / sd_g
  q = 1 - 1/(1+exp(-X2 %*% p2))
  for (iter in 2:iteration) {
    u = MH_u(Y2,Y3,propo[7],z,q,D,sd_g)
    D = MH_D(Y3,propo[6],z,q,D,sd_g)
    sd_g = MH_sig(Y3,propo[5],z,q,D,sd_g)
    z = (Y1 - X1 %*% p1) / sd_g
    p2 = MH_p2(Y3,X2,propo[3:4],p2,z,D,sd_g,Y0,X0)
    q = 1 - 1/(1+exp(-X2 %*% p2))
    p1 = MH_p1(Y1,X1,propo[1:2],p1,q,D,sd_g,Y3)
    z = (Y1 - X1 %*% p1) / sd_g
    param[iter,] = c(p1,p2,sd_g,D[1,2])
    if ((iter %% 50 == 0) & (iter < burn)) {
      freq = apply(param[1:iter,1:6], 2, function(x) length(unique(x))) / iter
      propo[1:6] = ifelse(freq > 0.4, propo[1:6] * 1.1, propo[1:6] * 0.9) 
    }
  }
  return(param)
}

MH_p1 = function(Y1,X1,propo_1,p1,q,D,sd_g,Y3) { # update fixed effect, variance w.r.t. continuous response by MH
  for (j in 1:2) { # update fixed effect parameters
    temp_p1 = p1
    temp_p1[j] = temp_p1[j] + rnorm(1,0,propo_1[j]) #propose parameter
    z_old = (Y1 - X1 %*% p1) / sd_g
    z_new = (Y1 - X1 %*% temp_p1) / sd_g
    a = den(Y3,q,z_new,sd_g,D) - den(Y3,q,z_old,sd_g,D)
    if (runif(1,0,1) < exp(a)) p1[j] = temp_p1[j]
  }
  return(p1)
} 

MH_sig = function(Y3,propo_s,z,q,D,sd_g) {
  sd_new = rlnorm(1,log(sd_g), propo_s) #propose parameter
  z_new = z * sd_g / sd_new
  # a = den(Y3,q,z_new,sd_new,D) - den(Y3,q,z,sd_g,D)
  a = sum(dnorm(z_new,log=TRUE)-dnorm(z,log=TRUE))-length(z)*(log(sd_new)-log(sd_g))
  if (runif(1) < exp(a)) sd_g = sd_new
  return(sd_g)
}

MH_p2 = function(Y3,X2,propo_2,p2,z,D,sd_g,Y0,X0) { # update fixed effect parameters w.r.t. binary response by MH
  for (j in 1:2) {
    temp_p2 = p2
    temp_p2[j] = temp_p2[j] + rnorm(1,0,propo_2[j]) # propose parameter
    q_old = 1 - 1/(1+exp(-X2 %*% p2))
    q_new = 1 - 1/(1+exp(-X2 %*% temp_p2))
    a = den(Y3,q_new,z,sd_g,D) - den(Y3,q_old,z,sd_g,D)
    # a: data likelihood, symmetric proposal, diffuse prior
    p_old = 1/(1+exp(-X0 %*% p2))
    p_new = 1/(1+exp(-X0 %*% temp_p2))
    a = a + sum(Y0 * (log(p_new) - log(p_old)) + (1 - Y0) * (log(1 - p_new) - log(1 - p_old)))
    if (runif(1,0,1) < exp(a)) p2[j] = temp_p2[j]
  }
  return(p2)
} 

MH_D = function(Y3,propo_a,z,q,D,sd_g) {
  alpha_old = D[1,2]
  if (abs(alpha_old) < 1 - propo_a) {alpha_new = runif(1, alpha_old - propo_a, alpha_old + propo_a)}
  else if (alpha_old > 1 - propo_a) {alpha_new = runif(1, 1 - 2*propo_a, 1)}
  else {alpha_new = runif(1,-1,-1 + 2*propo_a)}
  D_new = matrix(c(1, alpha_new, alpha_new, 1), nrow = 2)
  a = den(Y3,q,z,sd_g,D_new) - den(Y3,q,z,sd_g,D)
  # a: data likelihood, symmetric proposal, diffuse prior
  if (runif(1) < exp(a)) D = D_new
  return(D)
}

MH_u = function(Y2,Y3,propo_u,z,q,D,sd_g) {
  u = Y3 - Y2
  for (i in 1:length(u)) {
    u_old = u[i]
    if ((u_old < 1 - propo_u) & (u_old > 0 + propo_u)) {u_new = runif(1, u_old - propo_u, u_old + propo_u)}
    else if (u_old > 1 - propo_u) {u_new = runif(1, 1 - 2*propo_u, 1)}
    else {u_new = runif(1, 0, 2*propo_u)}
    a = den(Y2[i] + u_new, q[i], z[i], sd_g, D) - den(Y2[i] + u_old, q[i], z[i], sd_g, D)
    # a: data likelihood, symmetric proposal, diffuse prior
    if (runif(1) < exp(a)) u[i] = u_new
  }
  return(u)
}
###Model D1b: copula model###
