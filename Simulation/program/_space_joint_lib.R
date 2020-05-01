library(glmnet)

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
        result = space.neighbor(trainData, lam1=lam1, lam2=l2, lam3=lam3)
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
    result = space.neighbor(Data, lam1=lam1, lam2=l2,lam3=lam3)
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
        result = space.neighbor(Data, lam1=lam1, lam2=l2, lam3=lam3)
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
    result = space.neighbor(Data, lam1=lam1, lam2=l2,lam3=lam3)
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


space.joint.extBIC <- function(Data, imethod, lamb1.vec=exp(seq(-5, 0.5, 0.1)), lam2=0, lamb3.vec=exp(seq(-2,2,0.5)),seed=777){

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
              result = space.neighbor(Data, lam1=lam1, lam2=l2, lam3=lam3)
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
      result = space.neighbor(Data, lam1=lam1, lam2=l2,lam3=lam3)
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

lasso = function(X, y, Model.selection="ExtendedBIC"){
  g1 = glmnet(x=X, y=y)
  nn = nrow(X)
  p  = ncol(X)
  BIC.gamma = 1 - 1/(2*log(p)/log(nn))
  
  BICs = rep(NA, ncol(g1$beta))
  
  for(j in 1:ncol(g1$beta)){
    beta.j  = g1$beta[,j]
    nNon0   = length(which(beta.j != 0))
    resid.j = y - X %*% beta.j
    sigma.j = var(resid.j)*(nn-1)/nn
    
    if(Model.selection == "BIC"){
      BICs[j] = nn*log(sigma.j)+log(nn)*nNon0
    }else if(Model.selection == "ExtendedBIC"){
      tmp = 2*BIC.gamma*lchoose(p, nNon0)
      BICs[j] = nn*log(sigma.j)+log(nn)*nNon0 + tmp
    }
  }
  
  wmin = which.min(BICs)
  beta = g1$beta[,wmin]
  beta
}

# ----------------------------------------------------------------------------
# function for neighborhood selection with lasso
# ----------------------------------------------------------------------------

ne.PEN.lasso <- function(dat, V, order=FALSE, verbose=FALSE, 
                         Model.selection="ExtendedBIC")
{
  
  stopifnot(is.matrix(dat),nrow(dat)>1,ncol(dat)>1)
  p = ncol(dat)
  n = nrow(dat)
  coefNI = matrix(0, nrow=p, ncol=length(V))
  
  for (i in 1:length(V)) {
    v = V[i]  
    if (verbose) cat("variable=",v," ")
    X       = dat[,-v]
    y       = dat[,v]
    meanx   = apply(X, 2, mean)
    normx   = sqrt(rowSums((t(X)-meanx)^2)/n)
    meany   = mean(y)
    normy   = sqrt(sum((y-meany)^2)/n)
    X       = scale(X, meanx, normx)
    y       = (y - mean(y)) / normy
    
    # --------------------------------------------------------------
    # estimate skeleton by penalized regression
    # --------------------------------------------------------------
    
    if (order) {
      corr = abs(crossprod(X,y))
      o = order(corr,decreasing=T)
      o.back = order(o)
      betaHat = lasso(X[,o], y, Model.selection=Model.selection)
      coefNI[-v,i] = betaHat[o.back]
      if (verbose) cat(date(), "\n")
    }else {
      betaHat = lasso(X, y, Model.selection=Model.selection)
      coefNI[-v,i] = betaHat
      if (verbose) cat( date(), "\n")
    }
  }
  return(coefNI)
}

# Orable with FP+FN
space.joint.min.num <- function(Data, TrueData,imethod, lamb1.vec=exp(seq(-5, 0.5, 0.1)), lam2=0, lamb3.vec=exp(seq(-2,2,0.5)),seed=777){
  
  # True model
  Ak=TrueData
  w  = which(apply(Ak!=0,1,sum)>1) # nodes with co-parents
  moralA = Ak
  
  if(length(w) > 0){
    for (i in 1:length(w)) {
      moralA[t(combn(which(Ak[w[i],]!=0),2))] =1
    }
  }
  
  moralA[moralA + t(moralA)!=0]=1
  
  table(rowSums(Ak!=0)) # number of parents per node
  table(colSums(Ak!=0)) # number of parents/children per node
  
  trueGraph  = moralA
  table(trueGraph)
  
  table(rowSums( trueGraph !=0))
  table(colSums( trueGraph !=0))
  
  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  p    = ncol(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  TP.vec  = NULL
  TN.vec  = NULL
  FP.vec  = NULL
  FN.vec  = NULL
  nD.vec  = NULL
  nD0.vec = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  
  
  kk = 0
  lam1.ind=NULL
  lam3.ind=NULL
  
  for (lam1 in lamb1.vec){
    for (lam3 in lamb3.vec){
      lam1.ind=c(lam1.ind,lam1)
      lam3.ind=c(lam3.ind,lam3)
      cat(kk, lam1, date(), "\n")
      cat(kk, lam3, date(), "\n")
    
      if (imethod==1){
        result = space.neighbor(Data, lam1=lam1, lam2=l2,lam3=lam3)
      } else if (imethod==2){ 
        result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, iter=iter)
      } else if (imethod==3){
        result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=1, iter=iter)
      } else if (imethod==4){
        result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=2, iter=iter)
      }

      
      estimateGraph = (result$ParCor != 0)
      diag(estimateGraph) = 0
      diff = c(estimateGraph) - c(trueGraph)
      
      TP.vec = c(TP.vec,sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2)
      TN.vec = c(TN.vec,sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2)
      FP.vec = c(FP.vec,sum(diff == 1)/2)
      FN.vec = c(FN.vec,sum(diff == -1)/2)
      nD.vec = c(nD.vec,sum(c(estimateGraph))/2)
      nD0.vec = c(nD0.vec,sum(c(trueGraph))/2)
      
      
      pCor   = result$ParCor
      nNon0  = rowSums(pCor!=0) - 1
      
      sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
      
      for (j in 1:ncol(Data)){
        sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
      }
      
      Pred = Data %*% (pCor*sigma.mat) - Data
      Rsq=cor(c(Pred), c(Data))^2
  
      Rseq.vec=c(Rseq.vec, Rsq)
      Non0.vec=c(Non0.vec, nNon0)
    
    }
  }
    
    TPR=TP.vec/(TP.vec+FN.vec)
    FPR=FP.vec/(FP.vec+TN.vec)
    Error=FP.vec+FN.vec
    
    w2use = which.min(Error)
    lam1 = lam1.ind[w2use]
    lam3 = lam3.ind[w2use]
    
  if (imethod==1){
    result = space.neighbor(Data, lam1=lam1, lam2=l2,lam3=lam3)
  } else if (imethod==2){ 
    result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, iter=iter)
  } else if (imethod==3){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=1, iter=iter)
  } else if (imethod==4){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=2, iter=iter)
  }
  
  estimateGraph = (result$ParCor != 0)
  diag(estimateGraph) = 0
  diff = c(estimateGraph) - c(trueGraph)
  
  TP  = sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2
  TN  = sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2
  FP  = sum(diff == 1)/2
  FN  = sum(diff == -1)/2
  nD  = sum(c(estimateGraph))/2
  nD0 = sum(c(trueGraph))/2
  
  pCor   = result$ParCor
  nNon0  = rowSums(pCor!=0) - 1
  
  sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
  
  for (j in 1:ncol(Data)){
    sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
  }
  
  Pred = Data %*% (pCor*sigma.mat) - Data
  Rsq=cor(c(Pred), c(Data))^2

  list(result=result, nNon0=nNon0, R=Rsq, TP=TP, TN=TN, FP=FP, FN=FN, nD=nD, nD0=nD0, 
       imethod=imethod, lambda1=lam1, lambda3=lam3, lambdas=cbind(lam1.ind, lam3.ind), Non0.vec=Non0.vec, R.vec=Rseq.vec,
       TP.vec=TP.vec, TN.vec=TN.vec, FP.vec=FP.vec, FN.vec=FN.vec, nD.vec=nD.vec,
       nD0.vec=nD0.vec,   TPR=TPR, FPR=FPR, Error=Error)
}

# Oracle with FPR+FNR
space.joint.min.Rate <- function(Data, TrueData,imethod, lamb1.vec=exp(seq(-5, 0.5, 0.1)),lam2=0, lamb3.vec=exp(seq(-2,2,0.5)), seed=777){
  
  # True model
  Ak=TrueData
  w  = which(apply(Ak!=0,1,sum)>1) # nodes with co-parents
  moralA = Ak
  
  if(length(w) > 0){
    for (i in 1:length(w)) {
      moralA[t(combn(which(Ak[w[i],]!=0),2))] =1
    }
  }
  
  moralA[moralA + t(moralA)!=0]=1
  
  table(rowSums(Ak!=0)) # number of parents per node
  table(colSums(Ak!=0)) # number of parents/children per node
  
  trueGraph  = moralA
  table(trueGraph)
  
  table(rowSums( trueGraph !=0))
  table(colSums( trueGraph !=0))
  
  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  p    = ncol(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  TP.vec  = NULL
  TN.vec  = NULL
  FP.vec  = NULL
  FN.vec  = NULL
  nD.vec  = NULL
  nD0.vec = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  
  
  kk = 0
  lam1.ind=NULL
  lam3.ind=NULL
  
  for (lam1 in lamb1.vec){
    for (lam3 in lamb3.vec){
      lam1.ind=c(lam1.ind,lam1)
      lam3.ind=c(lam3.ind,lam3)
      cat(kk, lam1, date(), "\n")
      cat(kk, lam3, date(), "\n")
    
    if (imethod==1){
      result = space.neighbor(Data, lam1=lam1, lam2=l2, lam3=lam3)
    } else if (imethod==2){ 
      result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,iter=iter)
    } else if (imethod==3){
      result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,weight=1, iter=iter)
    } else if (imethod==4){
      result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,weight=2, iter=iter)
    }
    
    
    estimateGraph = (result$ParCor != 0)
    diag(estimateGraph) = 0
    diff = c(estimateGraph) - c(trueGraph)
    
    TP.vec = c(TP.vec,sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2)
    TN.vec = c(TN.vec,sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2)
    FP.vec = c(FP.vec,sum(diff == 1)/2)
    FN.vec = c(FN.vec,sum(diff == -1)/2)
    nD.vec = c(nD.vec,sum(c(estimateGraph))/2)
    nD0.vec = c(nD0.vec,sum(c(trueGraph))/2)
    
    
    pCor   = result$ParCor
    nNon0  = rowSums(pCor!=0) - 1
    
    sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
    
    for (j in 1:ncol(Data)){
      sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
    }
    
    Pred = Data %*% (pCor*sigma.mat) - Data
    Rsq=cor(c(Pred), c(Data))^2
    
    Rseq.vec=c(Rseq.vec, Rsq)
    Non0.vec=c(Non0.vec, nNon0)
    
    }
  }
  
  TPR=TP.vec/(TP.vec+FN.vec)
  FNR=FN.vec/(TP.vec+FN.vec)
  TNR=TN.vec/(TN.vec+FP.vec)
  FPR=FP.vec/(FP.vec+TN.vec)


  Error=FP.vec+FN.vec
  Error.rate = FPR+FNR
  
  w2use = which.min(Error.rate)
  lam1 = lam1.ind[w2use]
  lam3 = lam3.ind[w2use]
  
  if (imethod==1){
    result = space.neighbor(Data, lam1=lam1, lam2=l2,lam3=lam3)
  } else if (imethod==2){ 
    result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, iter=iter)
  } else if (imethod==3){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=1, iter=iter)
  } else if (imethod==4){
    result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=2, iter=iter)
  }
  
  estimateGraph = (result$ParCor != 0)
  diag(estimateGraph) = 0
  diff = c(estimateGraph) - c(trueGraph)
  
  TP  = sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2
  TN  = sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2
  FP  = sum(diff == 1)/2
  FN  = sum(diff == -1)/2
  nD  = sum(c(estimateGraph))/2
  nD0 = sum(c(trueGraph))/2
  
  pCor   = result$ParCor
  nNon0  = rowSums(pCor!=0) - 1
  
  sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
  
  for (j in 1:ncol(Data)){
    sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
  }
  
  Pred = Data %*% (pCor*sigma.mat) - Data
  Rsq=cor(c(Pred), c(Data))^2
  
  list(result=result, nNon0=nNon0, R=Rsq, TP=TP, TN=TN, FP=FP, FN=FN, nD=nD, nD0=nD0, 
       imethod=imethod, lambda1=lam1, lambda3=lam3, lambdas=cbind(lam1.ind, lam3.ind), Non0.vec=Non0.vec, R.vec=Rseq.vec,
       TP.vec=TP.vec, TN.vec=TN.vec, FP.vec=FP.vec, FN.vec=FN.vec, nD.vec=nD.vec,
       nD0.vec=nD0.vec, TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR,Error=Error, Error.rate=Error.rate)
}

# Oracle with F1
space.joint.F1 <- function(Data, TrueData,imethod, lamb1.vec=exp(seq(-5, 0.5, 0.1)),lam2=0, lamb3.vec=exp(seq(-2,2,0.5)), seed=777){
    
    # True model
    Ak=TrueData
    w  = which(apply(Ak!=0,1,sum)>1) # nodes with co-parents
    moralA = Ak
    
    if(length(w) > 0){
        for (i in 1:length(w)) {
            moralA[t(combn(which(Ak[w[i],]!=0),2))] =1
        }
    }
    
    moralA[moralA + t(moralA)!=0]=1
    
    table(rowSums(Ak!=0)) # number of parents per node
    table(colSums(Ak!=0)) # number of parents/children per node
    
    trueGraph  = moralA
    table(trueGraph)
    
    table(rowSums( trueGraph !=0))
    table(colSums( trueGraph !=0))
    
    ## make sure Data have been standardized
    Data = apply(Data, 2, scale)
    nn   = nrow(Data)
    p    = ncol(Data)
    
    l2   = 0
    iter = 3
    
    set.seed(seed)
    TP.vec  = NULL
    TN.vec  = NULL
    FP.vec  = NULL
    FN.vec  = NULL
    nD.vec  = NULL
    nD0.vec = NULL
    Non0.vec = NULL
    Rseq.vec = NULL
    
    
    kk = 0
    lam1.ind=NULL
    lam3.ind=NULL
    
    for (lam1 in lamb1.vec){
        for (lam3 in lamb3.vec){
            lam1.ind=c(lam1.ind,lam1)
            lam3.ind=c(lam3.ind,lam3)
            cat(kk, lam1, date(), "\n")
            cat(kk, lam3, date(), "\n")
            
            if (imethod==1){
                result = space.neighbor(Data, lam1=lam1, lam2=l2, lam3=lam3)
            } else if (imethod==2){
                result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,iter=iter)
            } else if (imethod==3){
                result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,weight=1, iter=iter)
            } else if (imethod==4){
                result = space.joint(Data, lam1=nn*lam1, lam2=l2, lam3=lam3,weight=2, iter=iter)
            }
            
            
            estimateGraph = (result$ParCor != 0)
            diag(estimateGraph) = 0
            diff = c(estimateGraph) - c(trueGraph)
            
            TP.vec = c(TP.vec,sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2)
            TN.vec = c(TN.vec,sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2)
            FP.vec = c(FP.vec,sum(diff == 1)/2)
            FN.vec = c(FN.vec,sum(diff == -1)/2)
            nD.vec = c(nD.vec,sum(c(estimateGraph))/2)
            nD0.vec = c(nD0.vec,sum(c(trueGraph))/2)
            
            
            pCor   = result$ParCor
            nNon0  = rowSums(pCor!=0) - 1
            
            sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
            
            for (j in 1:ncol(Data)){
                sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
            }
            
            Pred = Data %*% (pCor*sigma.mat) - Data
            Rsq=cor(c(Pred), c(Data))^2
            
            Rseq.vec=c(Rseq.vec, Rsq)
            Non0.vec=c(Non0.vec, nNon0)
            
        }
    }
    
    TPR=TP.vec/(TP.vec+FN.vec)
    FNR=FN.vec/(TP.vec+FN.vec)
    TNR=TN.vec/(TN.vec+FP.vec)
    FPR=FP.vec/(FP.vec+TN.vec)
    
    
    Error=FP.vec+FN.vec
    Error.rate = FPR+FNR
    precision = TP.vec/(TP.vec+FP.vec)
    recall = TPR
    F1 = 2*precision*recall/(precision+recall)
    
    w2use = which.max(F1)
    lam1 = lam1.ind[w2use]
    lam3 = lam3.ind[w2use]
    
    if (imethod==1){
        result = space.neighbor(Data, lam1=lam1, lam2=l2,lam3=lam3)
    } else if (imethod==2){
        result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, iter=iter)
    } else if (imethod==3){
        result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=1, iter=iter)
    } else if (imethod==4){
        result = space.joint(Data, lam1=nn*lam1, lam2=l2,lam3=lam3, weight=2, iter=iter)
    }
    
    estimateGraph = (result$ParCor != 0)
    diag(estimateGraph) = 0
    diff = c(estimateGraph) - c(trueGraph)
    
    TP  = sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2
    TN  = sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2
    FP  = sum(diff == 1)/2
    FN  = sum(diff == -1)/2
    nD  = sum(c(estimateGraph))/2
    nD0 = sum(c(trueGraph))/2
    
    pCor   = result$ParCor
    nNon0  = rowSums(pCor!=0) - 1
    
    sigma.mat = matrix(NA, nrow=ncol(Data), ncol=ncol(Data))
    
    for (j in 1:ncol(Data)){
        sigma.mat[,j] = sqrt(result$sig.fit/result$sig.fit[j])
    }
    
    Pred = Data %*% (pCor*sigma.mat) - Data
    Rsq=cor(c(Pred), c(Data))^2
    
    list(result=result, nNon0=nNon0, R=Rsq, TP=TP, TN=TN, FP=FP, FN=FN, nD=nD, nD0=nD0,
    imethod=imethod, lambda1=lam1, lambda3=lam3, lambdas=cbind(lam1.ind, lam3.ind), Non0.vec=Non0.vec, R.vec=Rseq.vec,
    TP.vec=TP.vec, TN.vec=TN.vec, FP.vec=FP.vec, FN.vec=FN.vec, nD.vec=nD.vec,
    nD0.vec=nD0.vec, TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR,Error=Error, F1=F1, w2use=w2use)
}

# Oracle NS-log with FP+FN or F1
PEN.min.num <- function(coefs, Data, TrueData, option='F1', seed=777){
  
  # True model
  Ak=TrueData
  w  = which(apply(Ak!=0,1,sum)>1) # nodes with co-parents
  moralA = Ak
  
  if(length(w) > 0){
    for (i in 1:length(w)) {
      moralA[t(combn(which(Ak[w[i],]!=0),2))] =1
    }
  }
  
  moralA[moralA + t(moralA)!=0]=1
  
  table(rowSums(Ak!=0)) # number of parents per node
  table(colSums(Ak!=0)) # number of parents/children per node
  
  trueGraph  = moralA
  table(trueGraph)
  
  table(rowSums( trueGraph !=0))
  table(colSums( trueGraph !=0))
  
  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  p    = ncol(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  TP.vec  = NULL
  TN.vec  = NULL
  FP.vec  = NULL
  FN.vec  = NULL
  nD.vec  = NULL
  nD0.vec = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  
  
  kk = 0
  for (lam1 in 1:dim(coefs)[3]){
    kk = kk + 1
    cat(kk, lam1, date(), "\n")
    
    result=coefs[,,lam1]
   
    
    estimateGraph = (result != 0)
    diag(estimateGraph) = 0
    diff = c(estimateGraph) - c(trueGraph)
    
    TP.vec = c(TP.vec,sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2)
    TN.vec = c(TN.vec,sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2)
    FP.vec = c(FP.vec,sum(diff == 1)/2)
    FN.vec = c(FN.vec,sum(diff == -1)/2)
    nD.vec = c(nD.vec,sum(c(estimateGraph))/2)
    nD0.vec = c(nD0.vec,sum(c(trueGraph))/2)
    
    pCor   = result
    nNon0  = rowSums(pCor!=0) - 1
    
    Pred = Data %*% (pCor) - Data
    Rsq=cor(c(Pred), c(Data))^2
    
    Rseq.vec=c(Rseq.vec, Rsq)
    Non0.vec=c(Non0.vec, nNon0)
    
  }
  
  TPR=TP.vec/(TP.vec+FN.vec)
  PPV=TP.vec/(TP.vec+FP.vec)
  NPV=TN.vec/(TN.vec+FN.vec)
  FPR=FP.vec/(FP.vec+TN.vec)
  
  FNR=1-TPR
  TNR=1-FPR

  
  Error=FP.vec+FN.vec
  Error.rate = FPR+FNR
  precision = PPV
  recall = TPR
  F1 = 2*precision*recall/(precision+recall)
  
  if (option=='F1'){
  w2use = which.max(F1)
  } else if (option=='FP+FN'){
  w2use = which.min(Error)  
  }
  
  
  list(Non0.vec=Non0.vec, R.vec=Rseq.vec,
       TP.vec=TP.vec, TN.vec=TN.vec, FP.vec=FP.vec, FN.vec=FN.vec, nD.vec=nD.vec,
       nD0.vec=nD0.vec, TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR,Error=Error, Error.rate=Error.rate, w2use=w2use)
}


ne.PEN.givenLambda <-
function (dat, nlambda, ntau, V, order = FALSE, verbose = FALSE, 
          Model.selection = "ExtendedBIC") 
{
  stopifnot(is.matrix(dat), nrow(dat) > 1, ncol(dat) > 1)
  p = ncol(dat)
  n = nrow(dat)
  meanx = apply(dat, 2, mean)
  normx = sqrt(rowSums((t(dat) - meanx)^2)/n)
  nDat = scale(dat, meanx, normx)
  corr = abs(crossprod(nDat, nDat))
  diag(corr) = 0
  lamMax = max(corr)
  thresh = 2 * exp(seq(log(lamMax), log(1/n), len = nlambda))
  lambda = thresh/10
  tau = 10^(seq(-6, -1, length.out = ntau))
  lambda = rep(lambda, each = ntau)
  tau = rep(tau, times = nlambda)
  key0 = sprintf("%.5e_%.5e", lambda, tau)
  if (Model.selection == "None") {
    coefNI = array(dim = c(p+1, length(V), length(lambda)))
  }
  else {
    coefNI = matrix(0, nrow = p+1, ncol = length(V))
  }
  for (i in 1:length(V)) {
    v = V[i]
    if (verbose) 
      cat("variable=", v, " ")
    X = nDat[, -v]
    y = nDat[, v]
    corXy = abs(drop(cor(X, y)))
    if (order) {
      o = order(corXy, decreasing = T)
      o.back = order(o)
      wp = PEN(X[, o], y, family = "gaussian", penalty = "LOG", 
               lambda = lambda, tau = tau, Model.selection = Model.selection)
      key1 = sprintf("%.5e_%.5e", wp$lambda, wp$tau)
      if (Model.selection == "None") {
        #wpBetas = wp$beta[-1, ]
        wpBetas = wp$beta
        coefNI[-(v+1), i, match(key1, key0)] = wpBetas[c(1,o.back+1),
                                                       ]
      }
      else {
        wpBetas = wp$beta
        coefNI[-(v+1), i] = wpBetas[o.back]
      }
      if (verbose) 
        cat(date(), "\n")
    }
    else {
      wp = PEN(X, y, family = "gaussian", penalty = "LOG", 
               lambda = lambda, tau = tau, Model.selection = Model.selection)
      key1 = sprintf("%.5e_%.5e", wp$lambda, wp$tau)
      if (Model.selection == "None") {
        coefNI[-(v+1), i, match(key1, key0)] = wp$beta
      }
      else {
        coefNI[-(v+1), i] = wp$beta
      }
      if (verbose) 
        cat(date(), "\n")
    }
  }
  return(list(coefNI=coefNI, lambda.vec=lambda, tau.vec=tau))
}

# Oracle NS-log by each gene with F1 or FN+FP

PEN.min.num.byrow <- function(coefs, ind, Data, TrueData, option='F1', seed=777){
  
  # True model
  Ak=TrueData
  w  = which(apply(Ak!=0,1,sum)>1) # nodes with co-parents
  moralA = Ak
  
  if(length(w) > 0){
    for (i in 1:length(w)) {
      moralA[t(combn(which(Ak[w[i],]!=0),2))] =1
    }
  }
  
  moralA[moralA + t(moralA)!=0]=1
  
  table(rowSums(Ak!=0)) # number of parents per node
  table(colSums(Ak!=0)) # number of parents/children per node
  
  trueGraph  = moralA[,ind]
 
  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  p    = ncol(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  TP.vec  = NULL
  TN.vec  = NULL
  FP.vec  = NULL
  FN.vec  = NULL
  nD.vec  = NULL
  nD0.vec = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  
  
  kk = 0
  for (lam1 in 1:dim(coefs)[3]){
    kk = kk + 1
    cat(kk, lam1, date(), "\n")
    
    result=coefs[,ind,lam1]
    result[is.na(result)]=0
    
    estimateGraph = (result != 0)
    diff = c(estimateGraph) - c(trueGraph)
    
    TP.vec = c(TP.vec,sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2)
    TN.vec = c(TN.vec,sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2)
    FP.vec = c(FP.vec,sum(diff == 1)/2)
    FN.vec = c(FN.vec,sum(diff == -1)/2)
    nD.vec = c(nD.vec,sum(c(estimateGraph))/2)
    nD0.vec = c(nD0.vec,sum(c(trueGraph))/2)
  }
  
  TPR=TP.vec/(TP.vec+FN.vec)
  PPV=TP.vec/(TP.vec+FP.vec)
  NPV=TN.vec/(TN.vec+FN.vec)
  FPR=FP.vec/(FP.vec+TN.vec)
  
  FNR=1-TPR
  TNR=1-FPR
  
  
  Error=FP.vec+FN.vec
  Error.rate = FPR+FNR
  precision = PPV
  recall = TPR
  F1 = 2*precision*recall/(precision+recall)
  
  if (option=='F1'){
    w2use = which.max(F1)
  } else if (option=='FP+FN'){
    w2use = which.min(Error)  
  }
  
  
  list(Non0.vec=Non0.vec, R.vec=Rseq.vec,
       TP.vec=TP.vec, TN.vec=TN.vec, FP.vec=FP.vec, FN.vec=FN.vec, nD.vec=nD.vec,
       nD0.vec=nD0.vec, TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR,Error=Error, Error.rate=Error.rate, w2use=w2use)
}

# Oracle NS-lasso by each gene with F1 or FN+FP
LASSO.min.num.byrow <- function(coefs, ind, Data, TrueData, option='F1', seed=777){
  
  # True model
  Ak=TrueData
  w  = which(apply(Ak!=0,1,sum)>1) # nodes with co-parents
  moralA = Ak
  
  if(length(w) > 0){
    for (i in 1:length(w)) {
      moralA[t(combn(which(Ak[w[i],]!=0),2))] =1
    }
  }
  
  moralA[moralA + t(moralA)!=0]=1
  
  table(rowSums(Ak!=0)) # number of parents per node
  table(colSums(Ak!=0)) # number of parents/children per node
  
  trueGraph  = moralA[,ind]
  
  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  p    = ncol(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  TP.vec  = NULL
  TN.vec  = NULL
  FP.vec  = NULL
  FN.vec  = NULL
  nD.vec  = NULL
  nD0.vec = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  
  
  kk = 0
  result=rep(NA, length=p)
  for (lam1 in 1:dim(coefs)[2]){
    kk = kk + 1
    cat(kk, lam1, date(), "\n")
    
    result[-ind]=coefs[,lam1]
    result[ind]=0
    
    estimateGraph = (result != 0)
    diff = c(estimateGraph) - c(trueGraph)
    
    TP.vec = c(TP.vec,sum(estimateGraph ==1 & trueGraph==1 & diff == 0))
    TN.vec = c(TN.vec,sum(estimateGraph ==0 & trueGraph==0 & diff == 0))
    FP.vec = c(FP.vec,sum(diff == 1))
    FN.vec = c(FN.vec,sum(diff == -1))
    nD.vec = c(nD.vec,sum(c(estimateGraph)))
    nD0.vec = c(nD0.vec,sum(c(trueGraph)))
  }
  
  TPR=TP.vec/(TP.vec+FN.vec)
  PPV=TP.vec/(TP.vec+FP.vec)
  NPV=TN.vec/(TN.vec+FN.vec)
  FPR=FP.vec/(FP.vec+TN.vec)
  
  FNR=1-TPR
  TNR=1-FPR
  
  
  Error=FP.vec+FN.vec
  Error.rate = FPR+FNR
  precision = PPV
  recall = TPR
  F1 = 2*precision*recall/(precision+recall)
  
  if (option=='F1'){
    w2use = which.max(F1)
  } else if (option=='FP+FN'){
    w2use = which.min(Error)  
  }
  
  
  list(Non0.vec=Non0.vec, R.vec=Rseq.vec,
       TP.vec=TP.vec, TN.vec=TN.vec, FP.vec=FP.vec, FN.vec=FN.vec, nD.vec=nD.vec,
       nD0.vec=nD0.vec, TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR,Error=Error, Error.rate=Error.rate, w2use=w2use)
}


ne.Lasso.givenLambda <-
  function (dat, nlambda, V, order = FALSE, verbose = FALSE, 
            Model.selection = "None") 
  {
    stopifnot(is.matrix(dat), nrow(dat) > 1, ncol(dat) > 1)
    p = ncol(dat)
    n = nrow(dat)
    meanx = apply(dat, 2, mean)
    normx = sqrt(rowSums((t(dat) - meanx)^2)/n)
    nDat = scale(dat, meanx, normx)
    corr = abs(crossprod(nDat, nDat))
    diag(corr) = 0
    lamMax = max(corr)
    thresh = 2 * exp(seq(log(lamMax), log(1/n), len = nlambda))
    lambda = thresh/10
    key0 = lambda
    if (Model.selection == "None") {
      coefNI = array(dim = c(p, length(V), length(lambda)))
    }
    else {
      coefNI = matrix(0, nrow = p, ncol = length(V))
    }
    for (i in 1:length(V)) {
      v = V[i]
      if (verbose) 
        cat("variable=", v, " ")
      X = nDat[, -v]
      y = nDat[, v]
      corXy = abs(drop(cor(X, y)))
      if (order) {
        o = order(corXy, decreasing = T)
        o.back = order(o)
        betaHat = glmnet(X[,o], y, lambda=lambda)
        lambda.vec = betaHat$lambda
        if (Model.selection == "None") {
          for (j in 1:length(lambda.vec)){
            coefNI[-v,i,j] = betaHat$beta[o.back,j]
          }
        }
      }  else {
          betaHat = glmnet(X, y, lambda=lambda)
          wpBetas = betaHat$beta
          for (j in 1:length(lambda.vec)){
            coefNI[-v,i,j] = wpBetas[,j]
          }
        }
        if (verbose) 
          cat(date(), "\n")
      }
    return(list(coefNI=coefNI, lambda.vec=lambda))
  }

# Oracle NS-lasso with F1 or FN+FP

LASSO.min.num <- function(coefs, Data, TrueData, option='F1', seed=777){
  
  # True model
  Ak=TrueData
  w  = which(apply(Ak!=0,1,sum)>1) # nodes with co-parents
  moralA = Ak
  
  if(length(w) > 0){
    for (i in 1:length(w)) {
      moralA[t(combn(which(Ak[w[i],]!=0),2))] =1
    }
  }
  
  moralA[moralA + t(moralA)!=0]=1
  
  table(rowSums(Ak!=0)) # number of parents per node
  table(colSums(Ak!=0)) # number of parents/children per node
  
  trueGraph  = moralA
  table(trueGraph)
  
  table(rowSums( trueGraph !=0))
  table(colSums( trueGraph !=0))
  
  ## make sure Data have been standardized 
  Data = apply(Data, 2, scale)
  nn   = nrow(Data)
  p    = ncol(Data)
  
  l2   = 0
  iter = 3
  
  set.seed(seed)
  TP.vec  = NULL
  TN.vec  = NULL
  FP.vec  = NULL
  FN.vec  = NULL
  nD.vec  = NULL
  nD0.vec = NULL
  Non0.vec = NULL
  Rseq.vec = NULL
  
  
  kk = 0
  for (lam1 in 1:dim(coefs)[3]){
    kk = kk + 1
    cat(kk, lam1, date(), "\n")
    
    result=coefs[,,lam1]
    
    estimateGraph = matrix(0,p,p)
    estimateGraph[result!=0 | t(result)!=0] =1
    diff = c(estimateGraph) - c(trueGraph)
    
    TP.vec = c(TP.vec,sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2)
    TN.vec = c(TN.vec,sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2)
    FP.vec = c(FP.vec,sum(diff == 1)/2)
    FN.vec = c(FN.vec,sum(diff == -1)/2)
    nD.vec = c(nD.vec,sum(c(estimateGraph))/2)
    nD0.vec = c(nD0.vec,sum(c(trueGraph))/2)
  }
  
  TPR=TP.vec/(TP.vec+FN.vec)
  PPV=TP.vec/(TP.vec+FP.vec)
  NPV=TN.vec/(TN.vec+FN.vec)
  FPR=FP.vec/(FP.vec+TN.vec)
  
  FNR=1-TPR
  TNR=1-FPR
  
  
  Error=FP.vec+FN.vec
  Error.rate = FPR+FNR
  precision = PPV
  recall = TPR
  F1 = 2*precision*recall/(precision+recall)
  
  if (option=='F1'){
    w2use = which.max(F1)
  } else if (option=='FP+FN'){
    w2use = which.min(Error)  
  }
  
  
  list(Non0.vec=Non0.vec, R.vec=Rseq.vec,
       TP.vec=TP.vec, TN.vec=TN.vec, FP.vec=FP.vec, FN.vec=FN.vec, nD.vec=nD.vec,
       nD0.vec=nD0.vec, TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR,Error=Error, Error.rate=Error.rate,F1=F1,w2use=w2use)
}

#---AR covariance structure---
ARcov = function(p, rho){
  Cov = matrix(0, p, p)
  for (i in 1 : p){
    for (j in 1 : p){
      Cov[i, j] = rho^(abs(i - j))
    }
  }
  return(Cov)
}

#---Block diagonal covariance structure 2---
BD2 = function(p, k, tau2 = 1, rho){
  # k is the number of blocks; rho is a vector, giving the lower and upper limits of the coefficients
  C = matrix(0, p, p)
  d = p / k
  for (m in 1 : k){
    rhotemp = runif(1, rho[1], rho[2])
    for (i in ((m - 1) * d + 1) : (m * d)) {
      for (j in ((m - 1) * d + 1) : (m * d)){
        if (i == j) C[i, j] = tau2
        else C[i, j] = rhotemp
      }	
    }
  }
  return(C)
}

#---Function for permuting the indices of the matrices---
sigP = function(Sigma){
  p = dim(Sigma)[1]
  permutation = sample(c(1 : p))
  res = matrix(0, p, p)
  for (i in 1 : p){
    for (j in 1 : p){
      res[i, j] = Sigma[permutation[i], permutation[j]]
    }
  }
  return(res)
} 


Eres<-function(coef.mat,nDat,cDat)
{
    n=dim(nDat)[1]
    p=dim(nDat)[2]
    CoefMatrix = matrix(0,p,p-1)
    Eresidual  = matrix(0,n,p)
    for (i in 1:p){
        x=nDat[,-i]
        y=cDat[,i]
        beta = coef.mat[,i]
        Eresidual[,i] = y - nDat%*%beta
        CoefMatrix[i,] = beta[-i] / apply(cDat[,-i], 2, sd) #  question: both X and Y has been scaled. I don't think we need this step....
        #CoefMatrix[i,] = beta
    }
    return(a=list(Eresidual=Eresidual, CoefMatrix=CoefMatrix))
}

