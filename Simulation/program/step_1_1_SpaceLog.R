args=commandArgs(TRUE)

if(length(args)==0){
  stop("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

imethod=as.numeric(args[[1]])
k=as.numeric(args[[2]])
n=as.numeric(args[[3]])
p=as.numeric(args[[4]])
e=as.numeric(args[[5]])
model=as.numeric(args[[6]])
#--set environment----------------------------------
vicky.dir = file.path('program')
source(file.path('_space_joint_lib.R'))

library(spacelog)
library(plyr)
library(dplyr)
library(data.table)
library(graphics)

setwd(vicky.dir)
imodel=ifelse(model==1, 'BA', 'ER')

if (imodel=='ER')
{
  if (n==400 & p==100) {
    ie=ifelse(e==1, 'pE_2.00e-02', 'pE_4.00e-02')
  } else if (n==400 & p==200){
    ie=ifelse(e==1, 'pE_1.00e-02', 'pE_2.00e-02')
  } else if (n==400 & p==300){
    ie=ifelse(e==1, 'pE_6.67e-03','pE_1.33e-02')
  }
load(file.path(vicky.dir,'simu_data', paste0(imodel, '_simu_n_',n, '_p_',p,'_',ie,'_min_beta_0.1_100A.RData')))
}

if (imodel=='BA'){
  load(file.path(vicky.dir,'simu_data', paste0(imodel, '_simu_n_',n, '_p_',p,'_e_',e,'_min_beta_0.1_100A.RData')))
}

ls()


Data = X[[k]]
TrueData=A[[k]]
nlambda=100
ntau=10
nDat = apply(Data, 2,scale)
dat=nDat
meanx = apply(dat, 2, mean)
normx = sqrt(rowSums((t(dat) - meanx)^2)/n)
nDat = scale(dat, meanx, normx)
corr = abs(crossprod(nDat, nDat))
diag(corr) = 0
lamMax = max(corr)
thresh = 2 * exp(seq(log(lamMax), log(1/n), len = nlambda))
lambda = thresh/10
tau = 10^(seq(-6, 0, length.out = ntau))


lamb1.vec=lambda
lamb3.vec=tau
seed=777

iter=3
nn=n

r1 = space.joint.extBIC(nDat, imethod, lamb1.vec, lam2=0, lamb3.vec, seed)
r2 = space.joint.BIC(nDat, imethod, lamb1.vec, lam2=0, lamb3.vec, seed)
r3 = space.joint.min.num(nDat, TrueData, imethod, lamb1.vec, lam2=0, lamb3.vec, seed)
r4 = space.joint.F1(nDat, TrueData, imethod, lamb1.vec, lam2=0, lamb3.vec, seed)

save(r1,r2,r3,r4, file=file.path(vicky.dir, 'result', paste0("space_joint_log_imethod_", imethod,"_iData_",k,"_imodel_",imodel,"_n_",n, "_p_",p,"_e_",e,"_100A.RData")))

