#sheep lambda approximation
replication = 500
core = 5 #number of cores for used, it took us around 2 hours to finish with 5 cores
source("fn3.R")
store1 = matrix(0, nrow = replication, ncol = 8) #store lambda
store2 = matrix(0, nrow = replication, ncol = 8) #store population mean mass
store3 = matrix(0, nrow = replication, ncol = 8) #store population variance mass
NAO  = read.table("NAOc", header = FALSE)
NAO1 = read.table("NAO" , header = FALSE)
NAO1 = cbind(NAO1[2:12,1], (rowSums(NAO1[2:12,2:4]) + NAO1[1:11,13]) / 4) #NAO data in 1986-1996
NAO3 = cbind(NAO[22:51,1], (rowSums(NAO[22:51,2:4]) + NAO[21:50,13]) / 4) #NAO data in 1990-2019
NAO5 = cbind(NAO[ 2:51,1], (rowSums(NAO[ 2:51,2:4]) + NAO[ 1:50,13]) / 4) #NAO data in 1970-2019

ISSample  = as.matrix(read.table(paste(getwd(), "/sample/ISSample" , sep=""), header = TRUE))
indSample = as.matrix(read.table(paste(getwd(), "/sample/indSample", sep=""), header = TRUE))
repSample = as.matrix(read.table(paste(getwd(), "/sample/repSample", sep=""), header = TRUE))
copSample = as.matrix(read.table(paste(getwd(), "/sample/copSample", sep=""), header = TRUE))
driSample = as.matrix(read.table(paste(getwd(), "/sample/driSample", sep=""), header = TRUE))
iitSample = as.matrix(read.table(paste(getwd(), "/sample/iitSample", sep=""), header = TRUE))
mitSample = as.matrix(read.table(paste(getwd(), "/sample/mitSample", sep=""), header = TRUE))
iiiSample = as.matrix(read.table(paste(getwd(), "/sample/iiiSample", sep=""), header = TRUE))
miiSample = as.matrix(read.table(paste(getwd(), "/sample/miiSample", sep=""), header = TRUE))
id0 = sample(1:dim(miiSample)[1], replication, replace = TRUE)
miiSample = miiSample[id1,] #we do sampling on model VIII here to reduce the memory
n1 = dim(ISSample)[1]

time = Sys.time()
id0 = sample(1:n1, replication, replace = TRUE)
id1 = sample(1:n1, replication, replace = TRUE)
temp = mclapply(1:replication, function(x) 
  lam_ind(50,1,4,ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
          indSample[id1[x],c(1,2,5)],indSample[id1[x],3:4]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,1] = temp[,1]; store2[,1] = temp[,2]; store3[,1] = temp[,3] #independent IPMs

id0 = sample(1:n1, replication, replace = TRUE)
id1 = sample(1:n1, replication, replace = TRUE)
temp = mclapply(1:replication, function(x) 
  SIM_mit(1000,5000,50,1,4,
          ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
          iitSample[id1[x],c(1,2,5)],iitSample[id1[x],3:4],c(0, iitSample[id1[x],c(6,7)])),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,2] = temp[,1]; store2[,2] = temp[,2]; store3[,2] = temp[,3] #uncorrelated random year effect models

id0 = sample(1:n1, replication, replace = TRUE)
id1 = sample(1:dim(iiiSample)[1], replication, replace = TRUE)
iiiSample = iiiSample[id1,]
temp = mclapply(1:replication, function(x) 
  SIM_mii(50,30,50,1,4,
          ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
          iiiSample[x,c(1,2,5)],iiiSample[x,3:4],c(0, iiiSample[x,c(6,7)])),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,3] = temp[,1]; store2[,3] = temp[,2]; store3[,3] = temp[,3]  #uncorrelated random individual effect models

id0 = sample(1:n1, replication, replace = TRUE)
id1 = sample(1:n1, replication, replace = TRUE)
temp = mclapply(1:replication, function(x) 
  lam_rep(50,1,4,ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
          repSample[id1[x],c(1,2,6)],repSample[id1[x],4:5],repSample[id1[x],3]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,4] = temp[,1]; store2[,4] = temp[,2]; store3[,4] = temp[,3] #reproduction conditional IPMs


id0 = sample(1:n1, replication, replace = TRUE)
id1 = sample(1:n1, replication, replace = TRUE)
temp = mclapply(1:replication, function(x) 
  lam_cop(50,1,4,ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
          copSample[id1[x],c(1,2,6)],copSample[id1[x],3:4],copSample[id1[x],5]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,5] = temp[,1]; store2[,5] = temp[,2]; store3[,5] = temp[,3] #copula models


id0 = sample(1:n1, replication, replace = TRUE)
id1 = sample(1:n1, replication, replace = TRUE)
temp = mclapply(1:replication, function(x) 
  SIM_dri1(1000,5000,50,1,4,
           ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
           driSample[id1[x],c(1,2,7)],driSample[id1[x],4:5],c(-0.019,1.09,driSample[id1[x],c(3,6)])),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,6] = temp[,1]; store2[,6] = temp[,2]; store3[,6] = temp[,3] #shared driver models

#shared driver models, with bootstrapping NAO
# temp = mclapply(1:replication, function(x) 
#   SIM_dri2(1000,5000,50,1,4,
#            ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
#            driSample[id1[x],c(1,2,7)],driSample[id1[x],4:5],c(0,0,driSample[id1[x],c(3,6)]),NAO5),
#   mc.cores = core)

id0 = sample(1:n1, replication, replace = TRUE)
id1 = sample(1:n1, replication, replace = TRUE)
temp = mclapply(1:replication, function(x) 
  SIM_mit(1000,5000,50,1,4,
          ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
          mitSample[id1[x],c(1,2,6)],mitSample[id1[x],3:4],mitSample[id1[x],c(5,7,8)]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,7] = temp[,1]; store2[,7] = temp[,2]; store3[,7] = temp[,3]  #correlated random year effect models

id0 = sample(1:n1, replication, replace = TRUE)
temp = mclapply(1:replication, function(x) 
  SIM_mii(50,30,50,1,4,
          ISSample[id0[x],3:4],ISSample[id0[x],c(1,2,5)],
          miiSample[x,c(1,2,6)],miiSample[x,3:4],miiSample[x,c(5,7,8)]),
  mc.cores = core)
temp = matrix(unlist(temp), ncol = 3, byrow = TRUE)
store1[,8] = temp[,1]; store2[,8] = temp[,2]; store3[,8] = temp[,3]  #uncorrelated random individual effect models

print(Sys.time() - time)
dir.create("result") #create a new directory to store log lambda
setwd(paste(getwd(), "/result", sep=""))
write.table(store1, "lambda", row.names = FALSE)
write.table(store2, "mean", row.names = FALSE)
write.table(store3, "variance", row.names = FALSE)
