# This is the second .R file you should run, given that you already ran run1.R to obtain the posterior samples.
# This .R include code for approximate the log lambda.
# The estimates of log lambda will be stored in a new directory name "result".
# The computational time depends on the number of available core. 
# It took us less than 1 hour to finish with 5 cores

source("fn3.R")
replication = 100
core = 5 #number of cores for used, it took us less than 1 hours to finish with 5 cores
store1 = matrix(0, nrow = replication, ncol = 8) #store lambda
store2 = matrix(0, nrow = replication, ncol = 8) #store population mean mass
store3 = matrix(0, nrow = replication, ncol = 8) #store population variance mass
NAO  = read.table("NAOc", header = FALSE)
NAO1 = read.table("NAO" , header = FALSE)
NAO1 = cbind(NAO1[2:12,1], (rowSums(NAO1[2:12,2:4]) + NAO1[1:11,13]) / 4) #NAO data in 1986-1996
NAO3 = cbind(NAO[22:51,1], (rowSums(NAO[22:51,2:4]) + NAO[21:50,13]) / 4) #NAO data in 1990-2019
NAO5 = cbind(NAO[ 2:51,1], (rowSums(NAO[ 2:51,2:4]) + NAO[ 1:50,13]) / 4) #NAO data in 1970-2019

# posterior sample of vital rate parameters
ISSample  = as.matrix(read.table(paste(getwd(), "/sample/ISSample" , sep=""), header = TRUE))
indSample = as.matrix(read.table(paste(getwd(), "/sample/indSample", sep=""), header = TRUE))
iitSample = as.matrix(read.table(paste(getwd(), "/sample/iitSample", sep=""), header = TRUE))
iiiSample = as.matrix(read.table(paste(getwd(), "/sample/iiiSample", sep=""), header = TRUE))
repSample = as.matrix(read.table(paste(getwd(), "/sample/repSample", sep=""), header = TRUE))
copSample = as.matrix(read.table(paste(getwd(), "/sample/copSample", sep=""), header = TRUE))
driSample = as.matrix(read.table(paste(getwd(), "/sample/driSample", sep=""), header = TRUE))
mitSample = as.matrix(read.table(paste(getwd(), "/sample/mitSample", sep=""), header = TRUE))
miiSample = as.matrix(read.table(paste(getwd(), "/sample/miiSample", sep=""), header = TRUE))
ISSample  = ISSample[sample(1:dim(ISSample)[1], replication, replace = TRUE),]
indSample = indSample[sample(1:dim(indSample)[1], replication, replace = TRUE),]
iitSample = iitSample[sample(1:dim(iitSample)[1], replication, replace = TRUE),]
iiiSample = iiiSample[sample(1:dim(iiiSample)[1], replication, replace = TRUE),]
repSample = repSample[sample(1:dim(repSample)[1], replication, replace = TRUE),]
copSample = copSample[sample(1:dim(copSample)[1], replication, replace = TRUE),]
driSample = driSample[sample(1:dim(driSample)[1], replication, replace = TRUE),]
mitSample = mitSample[sample(1:dim(mitSample)[1], replication, replace = TRUE),]
miiSample = miiSample[sample(1:dim(miiSample)[1], replication, replace = TRUE),]

time = Sys.time()
temp = mclapply(1:replication, function(idx) 
  lam_ind(50,1,4,ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
          indSample[idx,c(1,2,5)],indSample[idx,3:4]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,1] = temp[,1]; store2[,1] = temp[,2]; store3[,1] = temp[,3] # I1: vanilla IPMs

temp = mclapply(1:replication, function(idx) 
  SIM_mit(1000,10000,50,1,4,
          ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
          iitSample[idx,c(1,2,5)],iitSample[idx,3:4],c(0, iitSample[idx,c(6,7)])),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,2] = temp[,1]; store2[,2] = temp[,2]; store3[,2] = temp[,3] # I2: uncorrelated random year effect models

temp = mclapply(1:replication, function(idx) 
  SIM_mii(50,30,50,1,4,
          ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
          iiiSample[idx,c(1,2,5)],iiiSample[idx,3:4],c(0, iiiSample[idx,c(6,7)])),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,3] = temp[,1]; store2[,3] = temp[,2]; store3[,3] = temp[,3]  # I3: uncorrelated random individual effect models

temp = mclapply(1:replication, function(idx) 
  lam_rep(50,1,4,ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
          repSample[idx,c(1,2,6)],repSample[idx,4:5],repSample[idx,3]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,4] = temp[,1]; store2[,4] = temp[,2]; store3[,4] = temp[,3] # D1a: reproduction conditional IPMs

temp = mclapply(1:replication, function(idx) 
  lam_ind(50,1,4,ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
          copSample[idx,c(1,2,5)],copSample[idx,3:4]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,5] = temp[,1]; store2[,5] = temp[,2]; store3[,5] = temp[,3] # D1b: copula models

temp = mclapply(1:replication, function(idx) # NAO follows a normal distribution with mean -0.019, sd 1.09
  SIM_dri1(1000,10000,50,1,4,
           ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
           driSample[idx,c(1,2,7)],driSample[idx,4:5],c(-0.019,1.09,driSample[idx,c(3,6)])),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,6] = temp[,1]; store2[,6] = temp[,2]; store3[,6] = temp[,3] # D2a: shared driver models

#shared driver models, with bootstrapping NAO
# temp = mclapply(1:replication, function(idx) 
#   SIM_dri2(1000,10000,50,1,4,
#            ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
#            driSample[idx,c(1,2,7)],driSample[idx,4:5],c(0,0,driSample[idx,c(3,6)]),NAO5),
#   mc.cores = core)

temp = mclapply(1:replication, function(idx) 
  SIM_mit(1000,10000,50,1,4,
          ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
          mitSample[idx,c(1,2,6)],mitSample[idx,3:4],mitSample[idx,c(5,7,8)]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,7] = temp[,1]; store2[,7] = temp[,2]; store3[,7] = temp[,3]  # D2b: correlated random year effect models

temp = mclapply(1:replication, function(idx) 
  SIM_mii(50,30,50,1,4,
          ISSample[idx,3:4],ISSample[idx,c(1,2,5)],
          miiSample[idx,c(1,2,6)],miiSample[idx,3:4],miiSample[idx,c(5,7,8)]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,8] = temp[,1]; store2[,8] = temp[,2]; store3[,8] = temp[,3]  # D3: correlated random individual effect models

print(Sys.time() - time)
dir.create("result") #create a new directory to store log lambda
setwd(paste(getwd(), "/result", sep=""))
write.table(store1, "lambda", row.names = FALSE)
write.table(store2, "mean", row.names = FALSE)
write.table(store3, "variance", row.names = FALSE)
