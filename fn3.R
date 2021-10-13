# This .R includes the code and required libraries to approximate log lambda
# given that you have run run1.R to obtain the posterior samples of parameters.
# Instead of running this .R, you can run run2.R that will use the following functions.

library(MASS)
library(nimble)
library(lme4)
library(parallel)
library(mvtnorm)
###IPM (midpoint rule)###
# the midpoint rule is applicable for Model I1, Model D1a, Model D1b
bin_ker=function(x,param) { #kernel for bernoulli vital rate
  u=1/(1+exp(-param[1]-param[2]*x))
  return(u)
}
gau_ker=function(x_new,x,param) { #kernel for gaussian vital rate
  u=dnorm(x_new,param[1]+param[2]*x,param[3])
  return(u)
}
gau_adjust = function(G, mesh) { 
  #in case the probability of growth/ inheritance over the bound is positive
  idx_ev = which(colSums(G) < 1)
  for (id in idx_ev) {
    if (id <= mesh / 2) {
      G[1, id] = G[1, id] + 1 - sum(G[,id])
    } else {
      G[mesh, id] = G[mesh, id] + 1 - sum(G[,id])
    }
  }
  idx_over = which(colSums(G) > 1)
  for (id in idx_over) {G[, id] = G[, id] / sum(G[,id])}
  return(G)
}
gb_ker=function(x_new,x,param_g,param_r,param_c) { # growth kernel for breeder
  z = (x_new-param_g[1]-param_g[2]*x)*param_c[1]/param_g[3]
  p = 1/(1+exp(-param_r[1]-param_r[2]*x))
  delta = pnorm((qnorm(1-p)-z)/sqrt(1-param_c[1]**2))
  u=(1-delta)*dnorm(x_new,param_g[1]+param_g[2]*x,param_g[3])
  return(u)
}
gn_ker=function(x_new,x,param_g,param_r,param_c) { # growth kernel for non-breeder
  z = (x_new-param_g[1]-param_g[2]*x)*param_c[1]/param_g[3]
  p = 1/(1+exp(-param_r[1]-param_r[2]*x))
  delta = pnorm((qnorm(1-p)-z)/sqrt(1-param_c[1]**2))
  u=delta*dnorm(x_new,param_g[1]+param_g[2]*x,param_g[3])
  return(u)
}
lam_ind=function(mesh,min_size,max_size,param_s,param_i,param_g,param_r) { # Model I1: vanilla models
  h=(max_size-min_size)/mesh # step size
  b=min_size+c(0:mesh)*h # boundary points
  y=0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  S=bin_ker(y,param_s) #survival kernel
  R=bin_ker(y,param_r) #reproduction kernel
  G=h*outer(y,y,gau_ker,param=param_g) #growth kernel
  I=h*outer(y,y,gau_ker,param=param_i) #inheritance kernel
  G=gau_adjust(G, mesh)
  I=gau_adjust(I, mesh)
  
  for(i in 1:mesh) I[,i]=I[,i]*R[i]
  K=G+I # placeholder
  for(i in 1:mesh) K[,i]=K[,i]*S[i] 
  lam <- Re(eigen(K)$values[1])
  u <- Re(eigen(K)$vectors[,1])
  u <- u / sum(u)
  mu = sum(y*u)
  sig2 = sum(u*(y**2)) - mu**2
  return(c(log(lam),mu,sig2))
}

#Model D1a: reproduction conditional models
lam_rep=function(mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_e) { 
  h  =(max_size-min_size)/mesh # step size
  b  =min_size+c(0:mesh)*h # boundary points
  y  =0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  S  =bin_ker(y,param_s) #survival kernel
  R  =bin_ker(y,param_r) #reproduction kernel
  G_0=h*outer(y,y,gau_ker,param=param_g) #growth kernel for non-breeder
  G_1=h*outer(y,y,gau_ker,param=c(param_g[1]+param_e[1],param_g[2:3])) #growth kernel for breeder
  I  =h*outer(y,y,gau_ker,param=param_i) #inheritance kernel
  G_0=gau_adjust(G_0, mesh)
  G_1=gau_adjust(G_1, mesh)
  I  =gau_adjust(I, mesh)
  K  =matrix(0, nrow = mesh*2, ncol = mesh*2)
  
  for(i in 1:mesh) I[,i]=I[,i]*R[i]*S[i] 
  for(i in 1:2) K[1:mesh, 1:mesh + (i-1)*mesh] = I
  for(i in 1:mesh) G_0[,i]=G_0[,i]*S[i] 
  K[1:mesh + mesh, 1:mesh + mesh] = G_0
  
  for(i in 1:mesh) {
    G_1[,i]=G_1[,i]*R[i]*S[i] 
    G_0[,i]=G_0[,i]*(1-R[i]) 
  }
  K[1:mesh + mesh, 1:mesh] = G_0+G_1
  lam <- Re(eigen(K)$values[1])
  u <- Re(eigen(K)$vectors[,1])
  u <- u / sum(u)
  mu = sum(y*u)
  sig2 = sum(u*(y**2)) - mu**2
  return(c(log(lam),mu,sig2))
}

lam_cop=function(mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_c) { #Model D1b: copula models
  h =(max_size-min_size)/mesh # step size
  b =min_size+c(0:mesh)*h # boundary points
  y =0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  S =bin_ker(y,param_s) #survival kernel
  R =bin_ker(y,param_r) #reproduction kernel
  G1=h*outer(y,y,gb_ker,param_g=param_g,param_r=param_r,param_c=param_c) #growth kernel for breeder
  G0=h*outer(y,y,gn_ker,param_g=param_g,param_r=param_r,param_c=param_c) #growth kernel for non-breeder
  I =h*outer(y,y,gau_ker,param=param_i) #inheritance kernel
  I =gau_adjust(I , mesh)
  
  for(i in 1:mesh) I[,i]=I[,i]*R[i]
  K=G1+G0+I # placeholder
  for(i in 1:mesh) K[,i]=K[,i]*S[i] 
  lam <- Re(eigen(K)$values[1])
  u <- Re(eigen(K)$vectors[,1])
  u <- u / sum(u)
  mu = sum(y*u)
  sig2 = sum(u*(y**2)) - mu**2
  return(c(log(lam),mu,sig2))
}
###IPM (midpoint rule)###


###IPM (element-selection simulation)###
# the element-selection simulation is applicable for IPMs with temporal heterogeneity, Model I2, Model D2a, Model D2b

#matrix construction for shared drivers models (D2a)
mat_dri=function(mesh,min_size,max_size,S,I,param_g,param_r,param_d,ran) {
  h =(max_size-min_size)/mesh # step size
  b =min_size+c(0:mesh)*h # boundary points
  y =0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  R =bin_ker(y,c(param_r[1] + ran * param_d[4],param_r[2])) #reproduction kernel
  G =h*outer(y,y,gau_ker,param=c(param_g[1]+ran * param_d[3],param_g[2:3]))
  G =gau_adjust(G , mesh)
  
  for(i in 1:mesh) I[,i]=I[,i]*R[i]
  K=G+I # placeholder
  for(i in 1:mesh) K[,i]=K[,i]*S[i] 
  return(K)
}

#simulation for shared drivers models (D2a) with gaussian predictive distribution on NAO
SIM_dri1 = function(T0,T,mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_d) {
  h =(max_size-min_size)/mesh # step size
  b =min_size+c(0:mesh)*h # boundary points
  y =0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  
  S =bin_ker(y,param_s) #no need to rebuild survival and inheritance kernel at every iteration
  I =h*outer(y,y,gau_ker,param=param_i) #inheritance kernel
  I =gau_adjust(I , mesh)
  
  store = matrix(0, nrow = T, ncol = 3) #abundance, population mass mean, variance
  init = runif(mesh,0,1)
  init = init / sum(init)
  ran = rnorm(T, param_d[1], param_d[2])
  for (t in 1:T) {
    K=mat_dri(mesh,min_size,max_size,S,I,param_g,param_r,param_d,ran[t])
    init = K %*% init
    store[t,1] = sum(init)
    init = init / store[t,1]
    
    store[t,2] = sum(init*y)
    store[t,3] = sum(init*(y**2)) - store[t,2]**2
  }
  return(c(mean(log(tail(store[,1],-T0))), colMeans(store[,2:3])))
}

#simulation for shared drivers models (D2a) with bootstrapping NAO
SIM_dri2 = function(T0,T,mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_d,NAO) {
  h =(max_size-min_size)/mesh # step size
  b =min_size+c(0:mesh)*h # boundary points
  y =0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  
  S =bin_ker(y,param_s) #no need to rebuild survival and inheritance kernel at every iteration
  I =h*outer(y,y,gau_ker,param=param_i) #inheritance kernel
  I =gau_adjust(I , mesh)
  
  store = matrix(0, nrow = T, ncol = 3) #abundance, population mass mean, variance
  init = runif(mesh,0,1)
  init = init / sum(init)
  idx  = sample(1:dim(NAO)[1], T, replace = TRUE)
  for (t in 1:T) {
    K=mat_dri(mesh,min_size,max_size,S,I,param_g,param_r,param_d,NAO[idx[t],2])
    init = K %*% init
    store[t,1] = sum(init)
    init = init / store[t,1]
    
    store[t,2] = sum(init*y)
    store[t,3] = sum(init*(y**2)) - store[t,2]**2
  }
  return(c(mean(log(tail(store[,1],-T0))), colMeans(store[,2:3])))
}

#matrix construction for (correlated) random year effect models (I2/ I2b)
mat_mit=function(mesh,min_size,max_size,S,I,param_g,param_r,param_m) {
  h =(max_size-min_size)/mesh # step size
  b =min_size+c(0:mesh)*h # boundary points
  y =0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  ran = rmvnorm(1, c(0,0), matrix(c(param_m[2]**2, rep(param_m[1]*param_m[2]*param_m[3],2), param_m[3]**2), nrow=2))
  R =bin_ker(y,c(param_r[1]+ran[2],param_r[2])) #reproduction kernel
  G =h*outer(y,y,gau_ker,param=c(param_g[1]+ran[1],param_g[2:3]))
  G =gau_adjust(G , mesh)
  
  for(i in 1:mesh) I[,i]=I[,i]*R[i]
  K=G+I # placeholder
  for(i in 1:mesh) K[,i]=K[,i]*S[i] 
  return(K)
}

#simulation for (correlated) random year effect models (I2/ D2b)
SIM_mit = function(T0,T,mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_m) {
  h =(max_size-min_size)/mesh # step size
  b =min_size+c(0:mesh)*h # boundary points
  y =0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  
  S =bin_ker(y,param_s) #no need to rebuild survival and inheritance kernel at every iteration
  I =h*outer(y,y,gau_ker,param=param_i) #inheritance kernel
  I =gau_adjust(I , mesh)
  
  store = matrix(0, nrow = T, ncol = 3) #abundance, population mass mean, variance
  init = runif(mesh,0,1)
  init = init / sum(init)
  for (t in 1:T) {
    K=mat_mit(mesh,min_size,max_size,S,I,param_g,param_r,param_m)
    init = K %*% init
    store[t,1] = sum(init)
    init = init / store[t,1]
    
    store[t,2] = sum(init*y)
    store[t,3] = sum(init*(y**2)) - store[t,2]**2
  }
  return(c(mean(log(tail(store[,1],-T0))), colMeans(store[,2:3])))
}
###IPM (element-selection simulation)###

###IPM (intermediate method)###
# the intermediate method is applicable for IPMs with persistent individual heterogeneity, Model I3, Model D3
# Here because we know the distribution of z (random individual effects), 
# instead of using mid-point rule on z,
# we use mesh point that are of uniform quantile

#matrix construction for (correlated) random individual effect models (I3/ D3)
mat_mii=function(meshi,mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_m) {
  D = matrix(c(param_m[2]**2,rep(param_m[1]*param_m[2]*param_m[3],2),param_m[3]**2),nrow=2)
  yi = matrix(0, nrow = meshi**2, ncol = 2)
  for (i in 1:meshi) {
    for (j in 1:meshi) {
      yi[(i-1)*meshi+j,1] = qnorm((i*2-1)/(2*meshi)) #qnorm(runif(1,0+(i-1)/meshi,i/meshi))
      yi[(i-1)*meshi+j,2] = qnorm((j*2-1)/(2*meshi)) #qnorm(runif(1,0+(j-1)/meshi,j/meshi))
    }
  }
  yi = t(apply(yi, 1, function(x) t(chol(D)) %*% x))
  
  h=(max_size-min_size)/mesh # step size
  b=min_size+c(0:mesh)*h # boundary points
  y=0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  S=bin_ker(y,param_s) #survival kernel
  I=h*outer(y,y,gau_ker,param=param_i) #inheritance kernel
  I=gau_adjust(I, mesh)
  K = list()
  for (i in 1:meshi**2) {
    I2=I
    R=bin_ker(y,c(param_r[1]+yi[i,2],param_r[2])) #reproduction kernel 
    G=h*outer(y,y,gau_ker,param=c(param_g[1]+yi[i,1],param_g[2:3])) #growth kernel
    G=gau_adjust(G, mesh)
    for(k in 1:mesh) {
      I2[,k]=I[,k]*R[k]*S[k] 
      G[,k] =G[,k]*S[k]
    }
    K[[(i-1)*2+1]] = I2
    K[[(i-1)*2+2]] = G
  }
  return(K)
}

#simulation for (correlated) random individual effect models (I3/ D3)
SIM_mii = function(T,meshi,mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_m) {
  K=mat_mii(meshi,mesh,min_size,max_size,param_s,param_i,param_g,param_r,param_m)
  
  h=(max_size-min_size)/mesh # step size
  b=min_size+c(0:mesh)*h # boundary points
  y=0.5*(b[1:mesh]+b[2:(mesh+1)]) # mesh points, mid_pt rule
  
  store_N = numeric(T)
  ind = runif(mesh*meshi**2,0,1)
  ind = ind / sum(ind)
  for (t in 1:T) {
    new_ind = numeric(mesh*meshi**2)
    for (i in 1:meshi**2) {
      new_ind[1:mesh + (i-1)*mesh] = new_ind[1:mesh + (i-1)*mesh] + K[[(i-1)*2+2]] %*% ind[1:mesh + (i-1)*mesh]
      new_ind = new_ind + rep((K[[(i-1)*2+1]] %*% ind[1:mesh + (i-1)*mesh]) / meshi**2, meshi**2)
    }
    store_N[t] = sum(new_ind)
    ind = new_ind / store_N[t]
  }
  lam = store_N[t]
  u = ind
  mu = sum(y*u)
  sig2 = sum(u*(y**2)) - mu**2
  return(c(log(lam),mu,sig2))
}
###IPM (intermediate method)###
