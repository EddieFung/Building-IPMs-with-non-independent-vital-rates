# This .R includes the nimble code and required libraries for estimating vital rate parameters.
# Instead of running this .R, you can run run1.R that will use the following functions.
library(MASS)
library(nimble)
library(lme4)
library(parallel)
library(mvtnorm)
indString <- nimbleCode({ # Ind(a): independent model.
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

iitString <- nimbleCode({ # Ind(b): independent model with uncorrelated random year effects.
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

iiiString <- nimbleCode({ # Ind(c): independent model with uncorrelated random individual effects.
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

driString <- nimbleCode({ # M1: shared drivers model
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

mitString <- nimbleCode({ # M2: correlated random year effect model
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

miiString <- nimbleCode({ # M3: correlated random individual effect model
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

repString <- nimbleCode({ # M4: reproduction conditional model
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

copString <- nimbleCode({  # M5:copula model
  pred1[1:N1,1] <- X1[1:N1,1:2] %*% param1[1:2] 
  pred2[1:N1,1] <- X1[1:N1,1:2] %*% param2[1:2]
  pred3[1:N3,1] <- X3[1:N3,1:2] %*% param2[1:2] 
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) 
  for (i in 1:N1) Y2[i]~dbinom(iprobit((probit(1-expit(pred2[i,1]))-alpha/sigma_g*(Y1[i]-pred1[i,1]))/sqrt(1 - alpha**2)),1) 
  for (i in 1:N3) Y3[i] ~ dbinom(1/(1+exp(-pred3[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
  alpha~dunif(-1,1)
})

mtrString <- nimbleCode({ # M6 (in appendix): correlated random year effect model + reproduction conditional model
  pred1[1:N1,1] <- X1[1:N1,1:3] %*% param1[1:3] + Z1[1:N1,1:noI] %*% u[1:noI,1]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2] + Z2[1:N2,1:noI] %*% u[1:noI,2]
  for (i in 1:N1) {
    Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) 
  }
  for (i in 1:N2) {
    Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  }
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param1[3]~dnorm(0, sd = 100)
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

drrString <- nimbleCode({ # M7 (in appendix): shared drivers model + reproduction conditional model
  pred1[1:N1,1] <- X1[1:N1,1:4] %*% param1[1:4]
  pred2[1:N2,1] <- X2[1:N2,1:3] %*% param2[1:3]
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) #pred1[i]
  for (i in 1:N2) Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param1[3]~dnorm(0, sd = 100)
  param1[4]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  param2[3]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
})

dmiString <- nimbleCode({ # M8 (in appendix): correlated random individual effect model + shared drivers model
  pred1[1:N1,1] <- X1[1:N1,1:3] %*% param1[1:3] + Z1[1:N1,1:noI ] %*% u[1:noI,1]
  pred2[1:N2,1] <- X2[1:N2,1:3] %*% param2[1:3] + Z2[1:N2,1:noI ] %*% u[1:noI,2]
  pred3[1:N3,1] <- X3[1:N3,1:3] %*% param2[1:3] + Z3[1:N3,1:noI2] %*% v[1:noI2 ]
  for (i in 1:N1) Y1[i] ~ dnorm(pred1[i,1], sd = sigma_g) 
  for (i in 1:N2) Y2[i] ~ dbinom(1/(1+exp(-pred2[i,1])),1)
  for (i in 1:N3) Y3[i] ~ dbinom(1/(1+exp(-pred3[i,1])),1)
  
  param1[1]~dnorm(0, sd = 100)
  param1[2]~dnorm(0, sd = 100)
  param1[3]~dnorm(0, sd = 100)
  param2[1]~dnorm(0, sd = 100)
  param2[2]~dnorm(0, sd = 100)
  param2[3]~dnorm(0, sd = 100)
  sigma_g~dunif(0,5)
  W[1:2,1:2] ~ dwish(Omega[1:2,1:2],3)
  for (i in 1:noI ) u[i,1:2] ~ dmnorm(mu[1:2],W[1:2,1:2]) #W is a precision matrix
  for (i in 1:noI2) v[i]     ~ dnorm(0,sd = sigma_u2)
  
  W_inv[1:2,1:2] <- inverse(W[1:2,1:2])
  sigma_u1 <- sqrt(W_inv[1,1])
  sigma_u2 <- sqrt(W_inv[2,2])
  rho <- W_inv[1,2] / sigma_u1 / sigma_u2
})

tmiString <- nimbleCode({ # M9 (in appendix): correlated random individual effect model + correlated random year effect model
  pred1[1:N1,1] <- X1[1:N1,1:2] %*% param1[1:2] + Z1[1:N1,1:noI ] %*% u[1:noI,1] + P1[1:N1,1:noT] %*% p[1:noT,1]
  pred2[1:N2,1] <- X2[1:N2,1:2] %*% param2[1:2] + Z2[1:N2,1:noI ] %*% u[1:noI,2] + P2[1:N2,1:noT] %*% p[1:noT,2]
  pred3[1:N3,1] <- X3[1:N3,1:2] %*% param2[1:2] + Z3[1:N3,1:noI2] %*% v[1:noI2 ] + P3[1:N3,1:noT] %*% p[1:noT,2]
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
  
  V[1:2,1:2] ~ dwish(Omega[1:2,1:2],3)
  for (i in 1:noT) p[i,1:2] ~ dmnorm(mu[1:2],V[1:2,1:2]) 
  V_inv[1:2,1:2] <- inverse(V[1:2,1:2])
  sigma_v1 <- sqrt(V_inv[1,1])
  sigma_v2 <- sqrt(V_inv[2,2])
  rh2 <- V_inv[1,2] / sigma_v1 / sigma_v2
})
