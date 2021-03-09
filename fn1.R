library(MASS)
library(nimble)
library(lme4)
library(parallel)
library(mvtnorm)
indString <- nimbleCode({ #independent model
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
repString <- nimbleCode({ #reproduction conditional model
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
copString <- nimbleCode({  #copula model
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
driString <- nimbleCode({ #shared drivers model
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
iitString <- nimbleCode({ #uncorrelated random year effect model
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

mitString <- nimbleCode({ #correlated random year effect model
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

iiiString <- nimbleCode({ #uncorrelated random individual effect model
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

miiString <- nimbleCode({ #correlated random individual effect model
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
