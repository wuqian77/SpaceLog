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
library(data.table)
library(graphics)
library(devtools)
install_github("Sun-lab/penalized_estimation/PEN")
install_github("Sun-lab/penalized_estimation/PenPC")
install.packages("pcalg")
library(PEN)
library(pcalg)
library(PenPC)

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



set.seed(777)
Data = X[[k]]
TrueData=A[[k]]
nDat = apply(Data, 2,scale)

coefs_PEN = ne.PEN.givenLambda(dat=nDat,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=FALSE, Model.selection='None')
ind_F1 = PEN.min.num(coefs=coefs_PEN$coefNI[-1,,], Data=nDat, TrueData = TrueData,option='F1')
ind_num = PEN.min.num(coefs=coefs_PEN$coefNI[-1,,], Data=nDat, TrueData = TrueData,option='FP+FN')
coef_log_F1_all = coefs_PEN$coefNI[-1,,ind_F1$w2use]
coef_log_num_all = coefs_PEN$coefNI[-1,,ind_num$w2use]
#coef_log_extBIC  = ne.PEN(dat=nDat,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=FALSE)
#coef_log_BIC  = ne.PEN(dat=nDat,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=FALSE, Model.selection='BIC')
#coef_las_BIC = ne.PEN.lasso(dat=nDat, V=1:p, order=TRUE, verbose=FALSE, 
#                            Model.selection="BIC")
#coef_las_extBIC = ne.PEN.lasso(dat=nDat, V=1:p, order=TRUE, verbose=FALSE, 
#                               Model.selection="ExtendedBIC")

coef_log_BIC_new=NULL
coef_log_extBIC_new=NULL
coef_log_qiu=NULL
coef_log_F1=NULL
coef_log_num=NULL
lambda=NULL
for (j in 1:p){
BIC.vec=NULL
extBIC.vec=NULL
b.vec=NULL
nNon0.vec=NULL
nn=nrow(nDat)
espilon=0.01
s=sqrt(n/(log(p)^3))
qiu.lambda = sqrt(2*(2+espilon)*log(p/s)/nn)
num = PEN.min.num.byrow(coefs=coefs_PEN$coefNI[-1,,],ind=j,Data=nDat, TrueData = TrueData,option='FP+FN')
F1 = PEN.min.num.byrow(coefs=coefs_PEN$coefNI[-1,,],ind=j,Data=nDat, TrueData = TrueData,option='F1')

for (i in 1:1000){
coefs = coefs_PEN$coefNI[,j,i]
coefs[j+1]=1
nNon0  = sum(coefs!=0)-2
Pred = cbind(rep(1,n),nDat)%*%(coefs) - nDat[,j]
a=(Pred-nDat[,j])
b=mean(a^2)
b.vec[i]=b
nNon0.vec[i]=nNon0
BIC.vec[i]=nn*log(b)+log(nn)*nNon0
gamma  =  1-1/(2*log(p)/log(nn))
logtau =  lchoose((p-1), nNon0)
extBIC.vec[i]=nn*log(b)+log(nn)*nNon0+2*gamma*logtau
}
coef_log_BIC_new = cbind(coef_log_BIC_new, coefs_PEN$coefNI[-1,j,which.min(BIC.vec)])
coef_log_extBIC_new = cbind(coef_log_extBIC_new, coefs_PEN$coefNI[-1,j,which.min(extBIC.vec)])
coef_log_qiu = cbind(coef_log_qiu,coefs_PEN$coefNI[-1,j,which.min(abs(coefs_PEN$lambda.vec-qiu.lambda))])
coef_log_F1 = cbind(coef_log_F1,coefs_PEN$coefNI[-1,j,F1$w2use])
coef_log_num = cbind(coef_log_num,coefs_PEN$coefNI[-1,j,num$w2use])
lambda = rbind(lambda, c(BIC=which.min(BIC.vec), extBIC=which.min(extBIC.vec), 
           qiu=which.min(abs(coefs_PEN$lambda.vec-qiu.lambda)), 
           num=num$w2use, F1=F1$w2use))

}
coef_log_BIC_new[is.na(coef_log_BIC_new)]=0
coef_log_extBIC_new[is.na(coef_log_extBIC_new)]=0
coef_log_qiu[is.na(coef_log_qiu)]=0
coef_log_F1[is.na(coef_log_F1)]=0
coef_log_num[is.na(coef_log_num)]=0
coef_log_F1_all[is.na(coef_log_F1_all)]=0
coef_log_num_all[is.na(coef_log_num_all)]=0


save(coefs_PEN,coef_log_F1_all, coef_log_num_all, coef_log_BIC_new,coef_log_extBIC_new,coef_log_qiu,coef_log_F1,coef_log_num, lambda, 
    file=file.path(vicky.dir,'result',paste0("PEN_bylambda_iData_",k,"_imodel_",imodel, "_n_",n, "_ps_",p, "_max.edge_",e,"_100A.RData")))






# lasso
set.seed(77)
coef_lasso = array(list(),p)
lambda.mat = matrix(0,4*p,4*p)
qiu.ind=NULL
coef_lasso_BIC=matrix(0,p-1,p)
coef_lasso_extBIC=matrix(0,p-1,p)
coef_lasso_qiu=matrix(0,p-1,p)
coef_lasso_F1=matrix(0,p-1,p)
coef_lasso_num=matrix(0,p-1,p)
ind.lasso.list=NULL
for (i in 1:p){
  X=nDat[,-i]
  y=nDat[,i]
  betaHat = glmnet(x=X, y=y)
  nn=nrow(X)
  pp=ncol(X)
  coef_lasso[[i]]=betaHat$beta
  lambda.mat[1:length(betaHat$lambda),i]=betaHat$lambda
  
  num = LASSO.min.num.byrow(coefs=betaHat$beta,ind=i,Data=nDat, TrueData = TrueData,option='FP+FN')
  F1 = LASSO.min.num.byrow(coefs=betaHat$beta,ind=i,Data=nDat, TrueData = TrueData,option='F1')

  espilon=0.01
  s=sqrt(n/(log(pp)^3))
  qiu.lambda = sqrt(2*(2+espilon)*log(pp/s)/nn)
  qiu.ind[i] = which.min(abs(betaHat$lambda-qiu.lambda))
  
  BIC.gamma = 1 - 1/(2*log(pp)/log(nn))
  BICs = rep(NA, ncol(betaHat$beta))
  extBICs = rep(NA, ncol(betaHat$beta))
  
  for(j in 1:ncol(betaHat$beta)){
    beta.j  = betaHat$beta[,j]
    nNon0   = length(which(beta.j != 0))
    resid.j = y - X %*% beta.j
    sigma.j = var(resid.j)*(nn-1)/nn
    BICs[j] = nn*log(sigma.j)+log(nn)*nNon0
    tmp = 2*BIC.gamma*lchoose((pp-1), nNon0)
    extBICs[j] = nn*log(sigma.j)+log(nn)*nNon0 + tmp
  }
  
  BICs[which(BICs< -10^6)]=NA
  extBICs[which(extBICs< -10^6)]=NA

  coef_lasso_BIC[,i] = betaHat$beta[,which.min(BICs)]
  coef_lasso_extBIC[,i] = betaHat$beta[,which.min(extBICs)]
  coef_lasso_qiu[,i] = betaHat$beta[,which.min(abs(betaHat$lambda-qiu.lambda))]
  coef_lasso_F1[,i] = ifelse(sum(F1$w2use)<1,0,betaHat$beta[,F1$w2use])
  coef_lasso_num[,i] = betaHat$beta[,num$w2use]
  ind.lasso.list = rbind(ind.lasso.list, c(nlambda=length(betaHat$lambda), BIC=which.min(BICs), extBIC=which.min(extBICs), 
                           qiu=qiu.ind[i], 
                           num=num$w2use, F1=F1$w2use))
}

A=B=C=D=E=matrix(0,p,p)
A[,1]=c(0,coef_lasso_BIC[,1])
B[,1]=c(0,coef_lasso_extBIC[,1])
C[,1]=c(0,coef_lasso_qiu[,1])
D[,1]=c(0,coef_lasso_F1[,1])
E[,1]=c(0,coef_lasso_num[,1])

for (i in 2:(p-1)){
  A[,i]=c(coef_lasso_BIC[1:(i-1),i],0,coef_lasso_BIC[(i):(p-1),i])
  B[,i]=c(coef_lasso_extBIC[1:(i-1),i],0,coef_lasso_extBIC[(i):(p-1),i])
  C[,i]=c(coef_lasso_qiu[1:(i-1),i],0,coef_lasso_qiu[(i):(p-1),i])
  D[,i]=c(coef_lasso_F1[1:(i-1),i],0,coef_lasso_F1[(i):(p-1),i])
  E[,i]=c(coef_lasso_num[1:(i-1),i],0,coef_lasso_num[(i):(p-1),i])
}
A[,p]=c(coef_lasso_BIC[,p],0)
B[,p]=c(coef_lasso_extBIC[,p],0)
C[,p]=c(coef_lasso_qiu[,p],0)
D[,p]=c(coef_lasso_F1[,p],0)
E[,p]=c(coef_lasso_num[,p],0)

coef_lasso_BIC=A
coef_lasso_extBIC=B
coef_lasso_qiu=C
coef_lasso_F1=D
coef_lasso_num=E

coefs_lasso = ne.Lasso.givenLambda(dat=nDat,nlambda=100,V=1:p,order=TRUE,verbose=FALSE, Model.selection='None')


lasso_F1 = LASSO.min.num(coefs=coefs_lasso$coefNI, Data=nDat, TrueData = TrueData,option='F1')
lasso_num = LASSO.min.num(coefs=coefs_lasso$coefNI, Data=nDat, TrueData = TrueData,option='FP+FN')
coef_lasso_F1_all = coefs_lasso$coefNI[,,lasso_F1$w2use]
coef_lasso_num_all = coefs_lasso$coefNI[,,lasso_num$w2use]

#load(file=file.path(out.dir,'Simulation_new', paste0("NS_LASSO_bylambda_imethod_", imethod,"_iData_",k,"_imodel_",imodel, "_n_",n, "_p_",p, "_e_",e,"_100A.RData")))
save(coefs_lasso, coef_lasso_BIC,coef_lasso_F1_all, coef_lasso_num_all, coef_lasso_extBIC,coef_lasso_qiu, coef_lasso_F1,coef_lasso_num, ind.lasso.list, 
    file=file.path(vicky.dir,'result', paste0("NS_LASSO_bylambda_imethod_iData_",k,"_imodel_",imodel, "_n_",n, "_ps_",p, "_max.edge_",e,"_100A.RData")))






