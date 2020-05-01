##################################
######: R package "spacelog", which is based on R package "space" and add "log" penalty
######: to keep the consistence, we didn't change the function from space, except "space.joint" and "space.neighbor" and "jsrm"
###### R functions for fitting JSRM and MB:


##################################################################
############################ function for users
#################################################################



space.joint<-function(Y.m, lam1, lam2=0, lam3=0, sig=NULL, weight=NULL,iter=2)
{
########################parameters:
## Y.m: n (sample size) by p (number of variables) data matrix;
##            (In this function, Y will be first standardized to mean zero and l_2 norm 1)
## lam1:l_1 penalty paramter; corresponding to Y scaled as norm being one,
##     suggested range is: O(n^{3/2}\Phi^{-1}(1-\alpha/(2p^2))) for small $\alpha$ such as 0.1
## lam2: l_2 penalty paramter
## lam3: l_0 penalty paramter for log penalty
## sig: vector of {sigma^{ii}}; If sig==NULL, then set the initial value as sig=rep(1,p) and then iter>=2
##     remark: the scale of sig does not matter
## weight: if weight==NULL: equal weight; if weight==1, set weight=sig and iter>=2;
##        if weight==2. set weight==degree and iter>=2; otherwise, weight=user specified vector of length p
##        remark: the scale of weight does not matter. In this function, weight will be rescaled to have mean=1
## iter: number of iterations used for estimating \rho and \sigma; default==2; or user specified
#       note:  in some cases iter is forced to be at least 2.
#######################

####################### return value
## A list: the estimated \{\rho^{ij}\}: p by p matrix, and $\{\sigma^{ii}\}$: p by 1 vector
#######################

n=nrow(Y.m)
p=ncol(Y.m)
ITER=iter



################### preparation
if(!is.null(sig))
  { #### if the user specify "sig", do not need to update sig.
     SIG.update=F
     SIG=sig
  } else { #### otherwise, need to update sig in each iteration
     SIG.update=T
     SIG=rep(1, p)
  }

if(length(weight)==0 | (length(weight)>1 & length(weight)<p)) ### weight==NULL
  {
    WEIGHT.tag=0 ### equal weight
    WEIGHT=rep(1, p)
    WEIGHT.update=F
  }
if(length(weight)==1)
  {
    if(weight==1)
     {
      WEIGHT.tag=1 ### sig based weight
      WEIGHT=SIG
      WEIGHT.update=T
      ITER=max(2, iter)
     } else {
      WEIGHT.tag=2 ### degree based weight
      WEIGHT=rep(1,p)
      WEIGHT.update=T
      ITER=max(2, iter)
     }
  }
if(length(weight)==p)
  {
     WEIGHT.tag=3 ## prespecified weight
     WEIGHT=weight
     WEIGHT.update=F
  }
################## begin to iterate

for(i in 1:iter)
  {
    print(paste("iter=", i, sep=""))

    Y.u<-Y.m*matrix(sqrt(WEIGHT),n,p,byrow=TRUE)
    sig.u<-SIG/WEIGHT

    try(jsrm.result<-jsrm(Y.u,sig.u,n,p,lam1,lam2,lam3),silent=TRUE)
    ParCor.fit<-jsrm.result$beta.m
    diag(ParCor.fit)<-1

    coef<-ParCor.fit[upper.tri(ParCor.fit)]
    beta.cur<-Beta.coef(coef,SIG)


    if(!WEIGHT.update & !SIG.update)
     {
       break
     } else { ## update sig or weight
        if(SIG.update)
         {
             SIG<-InvSig.diag.new(Y.m,beta.cur)
         }
        if(WEIGHT.update)
         {
             if(WEIGHT.tag==1)
             {        #### sig based
               WEIGHT=SIG
             } else { #### degree based
               temp.w<-apply(abs(ParCor.fit)>1e-6,1,sum)
               temp.w<-temp.w+max(temp.w)
               WEIGHT<-temp.w/sum(temp.w)*p
             }
          }
     } ### end updating WEIGHT and SIG
  } ### end iteration

  result<-list(ParCor=ParCor.fit,sig.fit=SIG, jsrm=jsrm.result)
  return(result)
}



##########################
# the same function from space package (Peng et al., 2007)

space.neighbor<-function(Y.m, lam1, lam2=0)
{
########################parameters:
## Y.m: n (sample size) by p (number of variables) data matrix;
##            (In this function, Y will be first standardized to mean zero and l_2 norm 1)
## lam1:l_1 penalty paramter; corresponding to Y scaled as norm being one,
##     suggested range is: O(n^{3/2}\Phi^{-1}(1-\alpha/(2p^2))) for small $\alpha$ such as 0.1
## lam2: l_2 penalty paramter
#########################

####################### return value
## A list: the estimated \{\rho^{ij}\}: p by p matrix, and $\{\sigma^{ii}\}$: p by 1 vector
#######################

   #dyn.load("MBshoot.so")


   n=nrow(Y.m)
   p=ncol(Y.m)
   Y_data<-as.vector(t(Y.m))
   n_iter=500

   Rho_output<-rep(0, p*p)
   sigma_output<-rep(0, p)

   junk<-.C("MBshoot",
         as.integer(n),
         as.integer(p),
         as.single(lam1),
         as.single(lam2),
         as.single(Y_data),
         sigma_output=as.single(sigma_output),
         as.integer(n_iter),
         Rho_output=as.single(Rho_output)
         )

   rho.v<-junk$Rho_output
   rho.m<-matrix(rho.v, p, p, byrow=T)

   sigma.new=junk$sigma_output

   result<-list(ParCor=rho.m,sig.fit=sigma.new)
   return(result)
}



############################## function for example
######## generate starr based network: three hubs
# the same function from space package (Peng et al., 2007)

Hub.net<-function(p,hub.no=3,hub.size=16,umin=0.5,umax=1)
{

degree.gen<-sample(1:3,p,replace=T)
degree.gen[1:hub.no]<-hub.size          ##hub size

if(sum(degree.gen)/2 != round(sum(degree.gen)/2))
degree.gen[p]<-degree.gen[p]+1

net.adj<-RanGra(degree.gen,p)
diag(net.adj)<-0

ParCor.gen<-GenPar(net.adj,umin,umax,flip=TRUE,factor=1.5)
##generate original "paratial correlation"

##check whether the generatee partial correlation is p.d?
all(eigen(ParCor.gen)$values>0)                         ##p.d.?

##truncation the partial correlations to make small ones larger
thre<-0.1

ParCor.trun<-ParCor.gen
ParCor.trun[abs(ParCor.gen)<thre&abs(ParCor.gen)>0]<-thre
all(eigen(ParCor.trun)$values>0)                        ##still p.d.?

##get covariance matrix and save in Sig.Rdata
Sig<-GenCov(ParCor.trun)

##
return(Sig)
}


##################################################################
############################ internal functions
#################################################################
### use cross validation to select lambda
space.joint.cv <- function(Data, imethod, lamb1.vec=exp(seq(-5, 0.5, 0.1)),lam2=0, lam3=lam3,
                           nfold=10, seed=777){

  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  
  l2   = 0
  iter = 3
  
  # split datasets into 10 fold
  set.seed(seed)
  MSE.mat = NULL
  nNon0.max.mat = nNon0.min.mat = nNon0.median.mat = NULL
  
  kk = 0
  for (lam1 in lamb1.vec){
    
    kk = kk + 1
    cat(kk, lam1, date(), "\n")
    
    # Create 10 equally size folds
    folds   = cut(seq(1,nrow(Data)), breaks=nfold, labels=FALSE)
    MSE.vec = rep(NA, nfold)
    nNon0.max.vec = nNon0.min.vec = nNon0.median.vec = rep(NA, nfold)
    
    # Perform cross validation
    for(i in 1:nfold){
      
      cat(i, date(), "\n")
      
      testIndexes = which(folds==i)
      testData    = Data[testIndexes, ]
      trainData   = Data[-testIndexes, ]
      
      if (imethod==1){
        result = space.neighbor(trainData, lam1=lam1, lam2=l2)
      } else if (imethod==2){ 
        result = space.joint(trainData, lam1=nn*lam1, lam2=l2, lam3=lam3, iter=iter)
      } else if (imethod==3){
        result = space.joint(trainData, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=1, iter=iter)
      } else if (imethod==4){
        result = space.joint(trainData, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=2, iter=iter)
      }
      
      pCor   = result$ParCor
      nNon0  = rowSums(pCor!=0) - 1
      
      nNon0.max.vec[i] = max(nNon0)
      nNon0.min.vec[i] = min(nNon0)
      nNon0.median.vec[i] = median(nNon0)
      
      sigma.mat = matrix(NA, nrow=ncol(trainData), ncol=ncol(trainData))
      
      for (j in 1:ncol(trainData)){
        sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
      }
      
      Pred = testData %*% (pCor*sigma.mat) - testData
      MSE.vec[i] = sum((c(Pred) - c(testData))^2)/length(c(Pred))
    }
    
    MSE.mat = cbind(MSE.mat, MSE.vec)
    nNon0.max.mat    = cbind(nNon0.max.mat, nNon0.max.vec)
    nNon0.min.mat    = cbind(nNon0.min.mat, nNon0.min.vec)
    nNon0.median.mat = cbind(nNon0.median.mat, nNon0.median.vec)
  }
  
  colSums(MSE.mat)
  
  w2use = which.min(colSums(MSE.mat))
  lam1 = lamb1.vec[w2use]
  
  if (imethod==1){
    result = space.neighbor(Data, lam1=lam1, lam2=l2)
  } else if (imethod==2){ 
    result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,iter=iter)
  } else if (imethod==3){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,weight=1, iter=iter)
  } else if (imethod==4){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,weight=2, iter=iter)
  }
 
  list(result=result, imethod=imethod, lambda1=lam1, lambdas=lamb1.vec, 
       MSE=MSE.mat, nNon0.max=nNon0.max.mat, nNon0.min=nNon0.min.mat, 
       nNon0.median=nNon0.median.mat) 
}

### use BIC to select tuning parameter
space.joint.BIC <- function(Data, imethod, lamb1.vec=exp(seq(-5, 0.5, 0.1)), lam2=0, lamb3.vec=exp(seq(-2,2,0.5)), seed=777){

  ## make sure Data have been standardized 
  #Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  BIC.vec  = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  result.mat = NULL
  
  kk=0
  lam1.ind=NULL
  lam3.ind=NULL
  
  for (lam1 in lamb1.vec){
      for (lam3 in lamb3.vec){
        kk=kk+1
          lam1.ind=c(lam1.ind,lam1)
          lam3.ind=c(lam3.ind,lam3)
          cat(kk, lam1, date(), "\n")
          cat(kk, lam3, date(), "\n")
    
    # BIC
      if (imethod==1){
        result = space.neighbor(Data, lam1=lam1, lam2=l2)
      } else if (imethod==2){ 
        result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, iter=iter)
      } else if (imethod==3){
        result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=1, iter=iter)
      } else if (imethod==4){
        result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=2, iter=iter)
      }

      result.mat[[kk]]=result
      
      pCor   = result$ParCor
      nNon0  = rowSums(pCor!=0) - 1
      
      sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
      
      for (j in 1:ncol(Data)){
        sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
      }
      
      Pred = Data %*% (pCor*sigma.mat) - Data
      Rsq=cor(c(Pred), c(Data))^2
      a=(Pred-Data)
      b=apply(a, 2, function(x) sum(x^2))
      BIC=sum(nn*log(b)+log(nn)*nNon0)
      
      Rseq.vec=c(Rseq.vec, Rsq)
      Non0.vec=c(Non0.vec, nNon0)
      BIC.vec=c(BIC.vec, BIC)
      }
  }

  
  w2use = which.min(BIC.vec)
  lam1 = lam1.ind[w2use]
  lam3 = lam3.ind[w2use]
  
  if (imethod==1){
    result = space.neighbor(Data, lam1=lam1, lam2=l2)
  } else if (imethod==2){ 
    result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, iter=iter)
  } else if (imethod==3){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=1, iter=iter)
  } else if (imethod==4){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=2, iter=iter)
  }
  
  pCor   = result$ParCor
  nNon0  = rowSums(pCor!=0) - 1
  sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
  
  for (j in 1:ncol(Data)){
    sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
  }
  
  Pred = Data %*% (pCor*sigma.mat) - Data
  Rsq=cor(c(Pred), c(Data))^2
  a=(Pred-Data)
  b=apply(a, 2, function(x) sum(x^2))
  BIC=sum(nn*log(b)+log(nn)*nNon0)
 
  list(result=result, result.mat=result.mat, nNon0=nNon0,R=Rsq, BIC=BIC, imethod=imethod, lambda1=lam1, lambda3=lam3, lambdas=cbind(lam1.ind, lam3.ind),
       BIC.vec=BIC.vec, Non0.vec=Non0.vec, R.vec=Rseq.vec) 
}

### use extBIC to select tuning parameter
space.joint.extBIC <- function(Data, imethod, lamb1.vec, lam2, lamb3.vec,seed){

  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  p    = ncol(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  BIC.vec  = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  result.mat = NULL
  
  kk=0
  lam1.ind=NULL
  lam3.ind=NULL
  
  for (lam1 in lamb1.vec){
      for (lam3 in lamb3.vec){
        kk=kk+1
          lam1.ind=c(lam1.ind,lam1)
          lam3.ind=c(lam3.ind,lam3)
          cat(kk, lam1, date(), "\n")
          cat(kk, lam3, date(), "\n")
          
          # BIC
          if (imethod==1){
              result = space.neighbor(Data, lam1=lam1, lam2=l2)
          } else if (imethod==2){
              result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, iter=iter)
          } else if (imethod==3){
              result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=1, iter=iter)
          } else if (imethod==4){
              result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=2, iter=iter)
          }

                result.mat[[kk]]=result
      
      pCor   = result$ParCor
      nNon0  = rowSums(pCor!=0) - 1
      
      sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
      
      for (j in 1:ncol(Data)){
        sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
      }
      
      Pred = Data %*% (pCor*sigma.mat) - Data
      Rsq=cor(c(Pred), c(Data))^2
      a      =  (Pred-Data)
      b      =  apply(a, 2, function(x) sum(x^2))
      gamma  =  1-1/(2*log(p)/log(nn))
      logtau =  lchoose((p-1), nNon0)
      BIC=sum(nn*log(b)+log(nn)*nNon0)+sum(2*gamma*logtau)
      #BIC=sum(nn*log(b)+log(nn)*nNon0)
      
      Rseq.vec=c(Rseq.vec, Rsq)
      Non0.vec=c(Non0.vec, nNon0)
      BIC.vec=c(BIC.vec, BIC)      
    }
  }
  
  w2use = which.min(BIC.vec)
  lam1 = lam1.ind[w2use]
  lam3 = lam3.ind[w2use]
  
  if (imethod==1){
      result = space.neighbor(Data, lam1=lam1, lam2=l2)
  } else if (imethod==2){
      result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, iter=iter)
  } else if (imethod==3){
      result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=1, iter=iter)
  } else if (imethod==4){
      result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3, weight=2, iter=iter)
  }
  
  pCor   = result$ParCor
  nNon0  = rowSums(pCor!=0) - 1
  
  sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
  
  for (j in 1:ncol(Data)){
    sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
  }
  
  Pred = Data %*% (pCor*sigma.mat) - Data
  Rsq=cor(c(Pred), c(Data))^2
  a      =  (Pred-Data)
  b      =  apply(a, 2, function(x) sum(x^2))
  gamma  =  1-1/(2*log(p)/log(nn))
  logtau =  lchoose((p-1), nNon0)
  BIC=sum(nn*log(b)+log(nn)*nNon0)+sum(2*gamma*logtau)
  
  list(result=result, result.mat=result.mat, nNon0=nNon0,R=Rsq, BIC=BIC, imethod=imethod, lambda1=lam1, lambda3=lam3, lambdas=cbind(lam1.ind, lam3.ind),
  BIC.vec=BIC.vec, Non0.vec=Non0.vec, R.vec=Rseq.vec)
}


############################## (b). calling JSRM_log adaptive shooting functions in c
############# for shooting: data  is standardized to norm 1
############# for update sigma: data is standardized to sd 1;  both are done within the c code

#######JSRM_log
###### call C function directly in R

jsrm<-function(Y,sig.use,n,p,lam1,lam2,lam3, n_iter=1000)
{
lambda1=lam1
lambda2=lam2
lambda3=lam3
sigma_sr=sig.use^0.5

Beta_output<-rep(0, p*p)
Beta_1<-rep(0, p*p)
Beta_2<-rep(0, p*p)
Beta_3<-rep(0, p*p)
Beta_4<-rep(0, p*p)
Y_data<-as.vector(t(Y))
n_iter=5

#dyn.load("JSRM.so") ### compiled from "JSRM.c"

iter_count=0

junk<-.C("JSRM_log",
           as.integer(n),
           as.integer(p),
           as.single(lambda1),
           as.single(lambda2),
           as.single(lambda3),
           as.single(Y_data),
           as.single(sigma_sr),
           as.integer(n_iter),
           iter.count=as.integer(iter_count),
           beta.estimate=as.single(Beta_output),
           beta.1=as.single(Beta_1),
           beta.2=as.single(Beta_2),
           beta.3=as.single(Beta_3),
           beta.4=as.single(Beta_4)
         )

#print(junk$iter.count)
beta.v<-junk$beta.estimate
beta.m<-matrix(beta.v, p,p, byrow=T)
beta.1=matrix(junk$beta.1,p,p,byrow=T)
beta.2=matrix(junk$beta.2,p,p,byrow=T)
beta.3=matrix(junk$beta.3,p,p,byrow=T)
beta.4=matrix(junk$beta.4,p,p,byrow=T)
return(list(beta.m=beta.m,beta.1=beta.1,beta.2=beta.2,beta.3=beta.3,beta.4=beta.4))
}

########################################
##### Estimate diagonal sigma
##use the fact that 1/sig^{ii} is the residual variance in one versus all others setting
# the same function from space package (Peng et al., 2007)


InvSig.diag.new<-function(Y, Beta){
################### parameters
### Y:    n by p data matrix; should be standardized to sd 1;
### Beta: beta matrix for the regression model
 p=ncol(Y)
 Beta.temp=Beta
 diag(Beta.temp)=0
 esti.Y=Y%*%Beta.temp
 residue=Y-esti.Y
 result=apply(residue^2, 2, mean)
 return(1/result)
}


InvSig.diag.old<-function(Y, Beta){
################### parameters
### Y:    n by p data matrix; should be standardized to sd 1;
### Beta: beta matrix for the regression model
# the same function from space package (Peng et al., 2007)

 p<-ncol(Y)
 result<-numeric(p)
 for(i in 1:p){
  beta.cur<-Beta[,i]
  beta.cur<-matrix(beta.cur[-i])
  y.cur<-Y[,i]
  x.cur<-Y[,-i]
  e.cur<-y.cur-x.cur%*%beta.cur
  sig.ii<-mean(e.cur^2)
  result[i]<-1/sig.ii
 }
 return(result)
}


########################################
### given rho^{ij}, get beta[j,i]: each column for each variable
### beta[j,i]<-rho^{ij}sqrt(sig^{jj}/sig{ii})
# the same function from space package (Peng et al., 2007)


Beta.coef<-function(coef, sig.fit){
############## parameter
### coef: rho^{ij};
### sig.fit: sig^{ii}
 p<-length(sig.fit)
 result<-matrix(0,p,p)
 result[upper.tri(result)]<-coef
 result<-result+t(result)
 result<-diag(1/sqrt(sig.fit))%*%result%*%diag(sqrt(sig.fit))
 result<-t(result)
 return(result)
}


#######################################
####generate a random graph based on given degrees of nodes
# the same function from space package (Peng et al., 2007)


RanGra<-function(degree,p){
##p--number of nodes;
##degree--degrees of each node
##sum(degree) must be an even number

result<-matrix(0,p,p)
stub<-matrix(0,sum(degree),2)

count<-0
 for(i in 1:length(degree)){
  cur<-(count+1):(count+degree[i])
  stub[cur,1]<-i
  stub[cur,2]<-cur
  count<-count+degree[i]
  }

 index<-sample(1:sum(degree),sum(degree)/2)

 stub.1<-stub[index,]
 stub.2<-stub[-index,]
  for (i in 1:nrow(stub.1)){
    cur.1<-stub.1[i,1]
    cur.2<-stub.2[i,1]
    result[cur.1,cur.2]<-1
    result[cur.2,cur.1]<-1
  }
 return(result)

}



###################################
###get partial correlation matrix based on an adjacent network
# the same function from space package (Peng et al., 2007)


GenPar<-function(net.adj,umin,umax,flip=TRUE,factor=2){
##net.adj<-the adjacent network
##unim,umax,  the range of the orginal parcor
##filp=T means random sign of the parcor
 p<-nrow(net.adj)
 result<-matrix(0,p,p)
  for(i in 2:p){
   for(j in 1:(i-1)){
    cur<-runif(1,umin,umax)
      sign<-1
      if(flip) sign<-sample(c(-1,1),1)
    cur<-cur*sign
    result[i,j]<-cur*net.adj[i,j]
    result[j,i]<-result[i,j]
   }
 }

 diag.re<-matrix(apply(abs(result),1,sum),p,p)*factor+1e-6
 result<-result/diag.re
 result<-(result+t(result))/2
 diag(result)<-1

 return(result)

}


##################################################
## get covariance matrix from a partial correlation matrix
## make it p.d.
# the same function from space package (Peng et al., 2007)


GenCov<-function(ParCor.m){
temp<-(-ParCor.m)
diag(temp)<-1
Sig<-solve(ParCor.m)


p<-nrow(Sig)
for (i in 2:p){
 for (j in 1:(i-1)){
  Sig[i,j]<-Sig[i,j]/sqrt(Sig[i,i]*Sig[j,j]) ##standardize to sigma=1
  Sig[j,i]<-Sig[i,j]
 }
}
diag(Sig)<-1

#diagonose
D<-eigen(Sig)$values
if(!all(D>0)){
print("invalid covariance matrix")
return()
}

return(Sig)
}

