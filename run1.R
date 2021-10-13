# This is the first .R file you should run.
# This .R include code for estimating vital rate parameters by MCMC on nimble. 
# The posterior samples of vital rate parameters will be stored in a new directory name "sample".
# The computational time is a few minutes per model, except model D1b and D3,
# which may run around a day
source("fn1.R") #nimble/ MCMC string
source("fn2.R") #MCEM

###NAO preparation###
# NAO data is obtained from https://crudata.uea.ac.uk/cru/data/nao/nao.dat
# then we compute the mean in Dec, Jan, Feb, Mar
NAO = read.table("NAO")
NAO = cbind(NAO[1:11,1], (rowSums(NAO[2:12,2:4]) + NAO[1:11,13]) / 4)
colnames(NAO) = c("sheep.yr", "NAO")
###NAO preparation###

###Raw data preparation###
df = read.csv("SHEEP.csv")
df[,5] = log(df[,5]); df[,6] = log(df[,6]) #log body mass
df[df[,7]>1,7]=1 #ignore twin

df_s = df[!is.na(df[,4 ])&!is.na(df[,5]), c(5,4 ,3)] #survival only depends on log body mass
df_i = df[!is.na(df[,10])&!is.na(df[,5]), c(5,10,3)] #inheritance only depends on log body mass
df_i[,2] = log(df_i[,2]) #log body mass for offspring

df_g = df[!is.na(df[,6 ])&!is.na(df[,5]), c(5,6 ,1,2,3)] #drop unrelated column
df_r = df[!is.na(df[,7 ])&!is.na(df[,5]), c(5,7 ,1,2,3)]
df_g = merge(df_g, NAO, by = "sheep.yr", all.x = T)
df_r = merge(df_r, NAO, by = "sheep.yr", all.x = T)
df_g = df_g[,c(2,3,4,1,5,6)] #reorder the df
df_r = df_r[,c(2,3,4,1,5,6)]
###Raw data preparation###

###mle for initialization of MCMC###
fit_s = glm(df_s[,2]~df_s[,1], family = "binomial") #mle for survival 
fit_i =  lm(df_i[,2]~df_i[,1]) #mle for inheritance
fit_f_g =  lm(wtt1~wt, data = df_g) #mle for growth
fit_f_r = glm( rec~wt, data = df_r, family = "binomial") #mle for reproduction
fit_d_g =  lm(wtt1~wt+NAO, data = df_g) #mle for growth with NAO records
fit_d_r = glm( rec~wt+NAO, data = df_r, family = "binomial") #mle for reproduction with NAO records
fit_t_g =  lmer(wtt1~wt+(1|sheep.yr), data = df_g) #mle for growth with random year effects
fit_t_r = glmer( rec~wt+(1|sheep.yr), data = df_r, family = "binomial") #mle for reproduction with random year effects
###mle for initialization of MCMC###

###estimating parameters of inheritance and survival###
N1 = dim(df_i)[1]
N2 = dim(df_s)[1]
X1 = cbind(rep(1, N1), df_i[,1])
X2 = cbind(rep(1, N2), df_s[,1])
Y1 = df_i[,2]
Y2 = df_s[,2]

data  <- list(Y1=Y1,Y2=Y2)
const <- list(N1=N1,N2=N2,X1=X1[1:N1,1:2],X2=X2[1:N2,1:2])
inits <- list(param1=fit_i$coefficient,param2=fit_s$coefficient,
              pred1=X1[1:N1,1:2] %*% summary(fit_i)$coefficients[,1],
              pred2=X2[1:N2,1:2] %*% summary(fit_s)$coefficients[,1],
              sigma_g=sqrt(sum((Y1[1:N1]-fit_i$fitted.values)**2)/(N1-2))
)
model <- nimbleModel(code = indString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
ISSample <- runMCMC(modelMCMCC, nburnin = 20000, niter = 100000)
###estimating parameters of inheritance and survival###

###estimating parameters in the vanilla model (I1)###
N1 = dim(df_g)[1]
N2 = dim(df_r)[1]
X1 = cbind(rep(1, N1), df_g[,1])
X2 = cbind(rep(1, N2), df_r[,1])
Y1 = df_g[,2]
Y2 = df_r[,2]

data  <- list(Y1=Y1,Y2=Y2)
const <- list(N1=N1,N2=N2,X1=X1[1:N1,1:2],X2=X2[1:N2,1:2])
inits <- list(param1=fit_f_g$coefficient,param2=fit_f_r$coefficient,
              pred1=X1[1:N1,1:2] %*% summary(fit_f_g)$coefficients[,1],
              pred2=X2[1:N2,1:2] %*% summary(fit_f_r)$coefficients[,1],
              sigma_g=sqrt(sum((Y1[1:N1]-fit_f_g$fitted.values)**2)/(N1-2))
)
model <- nimbleModel(code = indString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
indSample <- runMCMC(modelMCMCC, nburnin = 20000, niter = 100000)
###estimating parameters in the vanilla model (I1)###

###estimating parameters in the reproduction conditional model (D1a)###
N2 = dim(df_r)[1]
X2 = cbind(rep(1, N2), df_r[,1])
Y2 = df_r[,2]

df_gr = merge(df_g, df_r, by = c("id", "sheep.yr", "wt", "age", "NAO"))
N1 = dim(df_gr)[1]
X1 = cbind(rep(1, N1), df_gr[,3], (df_gr[,4] <= 1) * df_gr[,7])
Y1 = df_gr[,6]; Y1 = Y1[!is.na(X1[,3])]
X1 = X1[!is.na(X1[,3]),]
fit_e_g = lm(Y1~X1[,2]+X1[,3])
N1 = dim(X1)[1]

data  <- list(Y1=Y1,Y2=Y2)
const <- list(N1=N1,N2=N2,X1=X1,X2=X2)
inits <- list(param1=fit_e_g$coefficient,param2=fit_f_r$coefficient,
              pred1=X1 %*% fit_e_g$coefficient,
              pred2=X2 %*% fit_f_r$coefficient,
              sigma_g=sqrt(sum((Y1-fit_e_g$fitted.values)**2/(N1-2)))
)
model <- nimbleModel(code = repString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
repSample <- runMCMC(modelMCMCC, nburnin = 20000, niter = 100000)
###estimating parameters in the reproduction conditional model (D1a)###

###estimating parameters in the copula model (D1b)###
df_gra = merge(df_g, df_r, by = c("id", "sheep.yr", "wt", "age", "NAO"), all = TRUE)

N1 = sum(!is.na(df_gra[,6]))
N3 = sum( is.na(df_gra[,6]))
X1 = cbind(numeric(N1)+1, as.matrix(df_gra)[!is.na(df_gra[,6]),3])
X3 = cbind(numeric(N3)+1, as.matrix(df_gra)[ is.na(df_gra[,6]),3])
Y1 = df_gra[!is.na(df_gra[,6]),6]
Y2 = 1 - df_gra[!is.na(df_gra[,6]),7]
Y3 = df_gra[ is.na(df_gra[,6]),7]

inits = c(fit_f_g$coefficients,fit_f_r$coefficients,sigma(fit_f_g),0)
propo = c(0.1,0.1,0.2,0.2,0.1,0.15,0.1) # sd for the proposal distributions
copSample = MCMC_cop(Y1, Y2, X1, X1, Y3, X3, propo, inits, 20000, 100000)
###estimating parameters in the copula model (D1b)###

###estimating parameters in the shared drivers model (D2a)###
N1 = dim(df_g)[1]
N2 = dim(df_r)[1]
X1 = cbind(rep(1, N1), df_g[,1], df_g[,6])
X2 = cbind(rep(1, N2), df_r[,1], df_r[,6])
Y1 = df_g[,2]
Y2 = df_r[,2]

data  <- list(Y1=Y1,Y2=Y2)
const <- list(N1=N1,N2=N2,X1=X1[1:N1,1:3],X2=X2[1:N2,1:3])
inits <- list(param1=fit_d_g$coefficient,param2=fit_d_r$coefficient,
              pred1=X1[1:N1,1:3] %*% summary(fit_d_g)$coefficients[,1],
              pred2=X2[1:N2,1:3] %*% summary(fit_d_r)$coefficients[,1],
              sigma_g=sqrt(sum((Y1[1:N1]-fit_d_g$fitted.values)**2)/(N1-2))
)

model <- nimbleModel(code = driString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
driSample <- runMCMC(modelMCMCC, nburnin = 20000, niter = 100000)
###estimating parameters in the shared drivers model (D2a)###

###estimating parameters in the uncorrelated random year effect model (I2)###
X1 = cbind(numeric(dim(df_g)[1])+1, as.matrix(df_g))
X2 = cbind(numeric(dim(df_r)[1])+1, as.matrix(df_r))
N1 = dim(X1)[1]
N2 = dim(X2)[1]
Z1 = model.matrix(~0+as.factor(X1[,'sheep.yr']))
Z2 = model.matrix(~0+as.factor(X2[,'sheep.yr']))
noI = dim(Z1)[2]
colnames(Z1) = 1:noI
colnames(Z2) = 1:noI
const <- list(N1=N1,N2=N2,noI=noI,Omega=matrix(c(0.1,0,0,0.1),nrow=2),mu=c(0,0),
              X1=X1[1:N1,1:2],X2=X2[1:N2,1:2],Z1=Z1[1:N1,1:noI],Z2=Z2[1:N2,1:noI])
data  <- list(Y1=X1[,3],Y2=X2[,3])
inits <- list(param1=summary(fit_t_g)$coefficients[,1],param2=summary(fit_t_r)$coefficients[,1],
              pred1=X1[1:N1,1:2] %*% summary(fit_t_g)$coefficients[,1],
              pred2=X2[1:N2,1:2] %*% summary(fit_t_r)$coefficients[,1],
              u=cbind(data.frame(ranef(fit_t_g))[,4], data.frame(ranef(fit_t_r))[,4]),
              sigma_g=as.data.frame(VarCorr(fit_t_g))[2,5],
              sigma_u1=as.data.frame(VarCorr(fit_t_g))[1,5],sigma_u2=as.data.frame(VarCorr(fit_t_r))[1,5]
)
model <- nimbleModel(code = iitString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g","sigma_u1","sigma_u2"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
iitSample <- runMCMC(modelMCMCC, nburnin = 20000, niter = 100000)
###estimating parameters in the uncorrelated random year effect model (I2)###

###estimating parameters in the correlated random year effect model (D2b)###
const <- list(N1=N1,N2=N2,noI=noI,Omega=matrix(c(0.001,0,0,0.001),nrow=2),mu=c(0,0),
              X1=X1[1:N1,1:2],X2=X2[1:N2,1:2],Z1=Z1[1:N1,1:noI],Z2=Z2[1:N2,1:noI])
data  <- list(Y1=X1[,3],Y2=X2[,3])
inits <- list(param1=summary(fit_t_g)$coefficients[,1],param2=summary(fit_t_r)$coefficients[,1],
              pred1=X1[1:N1,1:2] %*% summary(fit_t_g)$coefficients[,1],
              pred2=X2[1:N2,1:2] %*% summary(fit_t_r)$coefficients[,1],
              u=cbind(data.frame(ranef(fit_t_g))[,4], data.frame(ranef(fit_t_r))[,4]),
              sigma_g=as.data.frame(VarCorr(fit_t_g))[2,5],
              W=matrix(c(1/as.data.frame(VarCorr(fit_t_g))[1,4],0,0,1/as.data.frame(VarCorr(fit_t_r))[1,4]),
                       nrow=2),
              W_inv=matrix(c(as.data.frame(VarCorr(fit_t_g))[1,4],0,0,as.data.frame(VarCorr(fit_t_r))[1,4]),
                           nrow=2),
              sigma_u1=as.data.frame(VarCorr(fit_t_g))[1,5],sigma_u2=as.data.frame(VarCorr(fit_t_r))[1,5],
              rho=0
)
model <- nimbleModel(code = mitString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g","sigma_u1","sigma_u2","rho"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
mitSample <- runMCMC(modelMCMCC, nburnin = 20000, niter = 100000)
###estimating parameters in the correlated random year effect model (D2b)###

###estimating parameters in the uncorrelated random individual effect model (I3)###
idx = numeric(dim(df_g)[1])
for (i in 1:dim(df_r)[1]) {
  if (sum(df_g["id"] == df_r[i,3])) {
    idx[i] = TRUE
  } else {
    idx[i] = FALSE
  }
}
idx1 = which(idx == 1)
idx0 = which(idx == 0)

fit_i_g =  lmer(wtt1~wt+(1|id), data = df_g)
fit_i_r = glmer( rec~wt+(1|id), data = df_r[idx1,], family = "binomial")
fit_i_r2= glmer( rec~wt+(1|id), data = df_r[idx0,], family = "binomial")

X1 = cbind(numeric(dim(df_g)[1])+1, as.matrix(df_g))
X2 = cbind(numeric(length(idx1))+1, as.matrix(df_r[idx1,]))
X3 = cbind(numeric(length(idx0))+1, as.matrix(df_r[idx0,]))
N1 = dim(X1)[1]
N2 = dim(X2)[1]
N3 = dim(X3)[1]
Z1 = model.matrix(~0+as.factor(X1[,'id']))
Z2 = model.matrix(~0+as.factor(X2[,'id']))
Z3 = model.matrix(~0+as.factor(X3[,'id']))
noI = dim(Z1)[2]
noI2= dim(Z3)[2]
colnames(Z1) = 1:noI
colnames(Z2) = 1:noI
colnames(Z3) = 1:noI2
Y1=X1[,3]
Y2=X2[,3]
Y3=X3[,3]
X1=X1[,1:2]
X2=X2[,1:2]
X3=X3[,1:2]

const <- list(N1=N1,N2=N2,N3=N3,noI=noI,noI2=noI2,
              X1=X1,X2=X2,X3=X3,
              Z1=Z1,Z2=Z2,Z3=Z3,
              Omega=matrix(c(0.001,0,0,0.001),nrow=2),mu=c(0,0)
)
data  <- list(Y1=Y1,Y2=Y2,Y3=Y3)
inits <- list(param1=summary(fit_i_g)$coefficients[,1],param2=summary(fit_i_r)$coefficients[,1],
              pred1=X1[1:N1,1:2] %*% summary(fit_i_g)$coefficients[,1],
              pred2=X2[1:N2,1:2] %*% summary(fit_i_r)$coefficients[,1],
              u=cbind(data.frame(ranef(fit_i_g))[,4], data.frame(ranef(fit_i_r))[,4]),
              v=data.frame(ranef(fit_i_r2))[,4],
              sigma_g=as.data.frame(VarCorr(fit_i_g))[2,5],
              sigma_u1=as.data.frame(VarCorr(fit_i_g))[1,5],sigma_u2=as.data.frame(VarCorr(fit_i_r))[1,5]
)
model <- nimbleModel(code = iiiString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g","sigma_u1","sigma_u2"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
iiiSample <- runMCMC(modelMCMCC, nburnin = 40000, niter = 200000)
###estimating parameters in the uncorrelated random individual effect model (I3)###

###estimating parameters in the correlated random individual effect model (D3)###
# To increase the speed of convergence, we initialized the chain around the mle,
# which is approximated by MCEM.
# You can use the below code to run the MCEM to get the approximated mle, but it may take around 2-3 hours
# Alternatively, you can simply use the numeric results we obtained from MCEM

###initialization from MCEM###
# param1=summary(fit_i_g)$coefficients[,1]
# param2=fit_f_r$coefficients
# sigma_g=as.data.frame(VarCorr(fit_i_g))[2,5]
# rho = 0
# sigma_u1 = as.data.frame(VarCorr(fit_i_g))[1,5]
# sigma_u2 = as.data.frame(VarCorr(fit_i_r))[1,5]
# size = 20000; df = 25; max_iter = 50
# p_store = EM_SIR(Y1, X1, Z1, Y2, X2, Z2, Y3, X3, Z3,
#                  rho, sigma_u1, sigma_u2, param1, sigma_g, param2,
#                  noI, noI2, size, df, max_iter)
# param1=c(tail(p_store[,4],1),tail(p_store[,5],1))
# param2=c(tail(p_store[,7],1),tail(p_store[,8],1))
# sigma_g=tail(p_store[,6],1)
# rho = tail(p_store[,1],1)
# sigma_u1 = tail(p_store[,2],1)
# sigma_u2 = tail(p_store[,3],1)
# D = matrix(c(sigma_u1**2, rep(rho*sigma_u1*sigma_u2, 2),sigma_u2**2), nrow = 2)
# D_inv = solve(D)

# Initialization for MCMC, which is the approximated mle obtained by MCEM
rho = -0.8122 
sigma_u1 = 0.0474; sigma_u2 = 0.643
param1 = c(1.619, 0.487); sigma_g = 0.0694
param2 = c(-9.167, 2.753)
D = matrix(c(sigma_u1**2, rep(rho*sigma_u1*sigma_u2, 2),sigma_u2**2), nrow = 2)
D_inv = solve(D)
###initialization from MCEM###

const <- list(N1=N1,N2=N2,N3=N3,noI=noI,noI2=noI2,
              X1=X1,X2=X2,X3=X3,
              Z1=Z1,Z2=Z2,Z3=Z3,
              Omega=matrix(c(0.001,0,0,0.001),nrow=2),mu=c(0,0)
)
data  <- list(Y1=Y1,Y2=Y2,Y3=Y3)
inits <- list(param1=param1,param2=param2,
              pred1=X1 %*% param1,
              pred2=X2 %*% param2,
              pred3=X3 %*% param2,
              u=sir_mean12(D_inv, sigma_g, Y1, Z1, Y2, Z2,X1 %*% param1, X2 %*% param2, noI),
              v=sir_mean3(sigma_u2, Y3, Z3, X3 %*% param2, noI2),
              sigma_g=sigma_g, W=D_inv, W_inv=D,
              sigma_u1=sigma_u1,sigma_u2=sigma_u2,rho=rho
)
model <- nimbleModel(code = miiString, name = "model", constants = const,
                     data = data, inits = inits)
modelMCMC <- buildMCMC(model,monitors = c("param1","param2","sigma_g","sigma_u1","sigma_u2","rho"))
modelC    <- compileNimble(model)
modelMCMCC<- compileNimble(modelMCMC, project = model)
miiSample <- runMCMC(modelMCMCC, nburnin = 200000, niter = 1000000) 
###estimating parameters in the correlated random individual effect model (D3)###


dir.create("sample") #create a new directory to store poterior samples of parameters
setwd(paste(getwd(), "/sample", sep=""))
write.table(ISSample,  "ISSample", row.names = FALSE)
write.table(indSample, "indSample", row.names = FALSE)
write.table(driSample, "driSample", row.names = FALSE)
write.table(repSample, "repSample", row.names = FALSE)
write.table(copSample, "copSample", row.names = FALSE)
write.table(iitSample, "iitSample", row.names = FALSE)
write.table(mitSample, "mitSample", row.names = FALSE)
write.table(iiiSample, "iiiSample", row.names = FALSE)
write.table(miiSample, "miiSample", row.names = FALSE)
