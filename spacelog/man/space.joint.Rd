\name{space.joint}
\alias{space.joint}
\title{A function to estimate partial correlations using the Joint Sparse Regression Model with LASSO, Elastic Net or Log penalty}
\description{A function to estimate partial correlations using the Joint Sparse Regression Model with LASSO, Elastic Net or Log penalty}
   
\usage{
space.joint(Y.m, lam1, lam2=0, lam3=0, sig=NULL, weight=NULL,iter=2)
}

\arguments{
    \item{Y.m}{numeric matrix. Columns are for variables and rows are for samples. 
                     Missing values are not allowed. It's recommended to first standardize each column to have mean 0 and \eqn{l_2} norm 1. 
                     }
    \item{lam1}{numeric value. This is the \eqn{l_1} norm penalty parameter. If the columns of Y.m have norm one, 
                then the suggested range of lam1 is: \eqn{O(n^{3/2}\Phi^{-1}(1-\alpha/(2p^2)))} 
                for small \eqn{\alpha} such as 0.1.
                } 
    \item{lam2}{numeric value. If not specified, lasso regression is used in the Joint Sparse Regression Model (JSRM). 
                Otherwise, elastic net regression is used in JSRM and \var{lam2} serves as the \eqn{l_2} norm penalty parameter.
                } 

    \item{lam3}{numeric value. If not specified lam3=0, it's the same as space.joint function in space package for Joint Sparse Regression Model (JSRM) using lasso or elastic net penalty. 
                Otherwise, log penalty is used in JSRM and \var{lam3} serves as the \eqn{l_0} norm penalty parameter.
                }          
    
    \item{sig}{numeric vector. Its length should be the same as the number of columns of \code{Y.m}. It is the vector of \eqn{\sigma^{ii}} (the diagonal of the inverse covariance matrix).
                If not specified, \eqn{\sigma^{ii}} will be estimated during the model fitting with initial 
                values \var{rep(1,p)}. The number of the iteration of the model fitting (\var{iter}) will then be at least 2. 
                Note, the scale of \var{sig} does not matter.
               } 
    \item{weight}{numeric value or vector. It specifies the weights or the type of weights used for each regression in JSRM. 
                  The default value is NULL, which means all regressions will be weighted equally in the joint model. If \var{weight}\eqn{=1},
                  residue variances will be used for weights. If \var{weight}\eqn{=2}, the estimated degree of each variable will be used for weights. 
                  Otherwise, it should be a positive numeric vector, whose length is equal to the number of columns of \var{Y.m}.
                  }           
    \item{iter}{integer. It is the total number of interactions in JSRM for estimating \eqn{\sigma^{ii}} and partial correlations. 
                 When \var{sig}\eqn{=NULL} and/or \var{weight}\eqn{=NULL} or 2,
                \var{iter} should be at least 2.}              
          }        
\details{
    \code{space.joint} adapted log penalty in partial correlations under the high-dimension-low-sample-size setting (Wu, Sun, and Hsu, 2020). It's an extension to the space paper (Peng and et al., 2007), which uses LASSO or Elastic Net penalty. 
    }

\value{
    A list with two components
    \item{ParCor}{the estimated partial correlation matrix.}
    \item{sig.fit}{numeric vector of the estimated diagonal \eqn{\sigma^{ii}}.} 
       }

\references{
J. Peng, P. Wang, N. Zhou, J. Zhu (2007), Partial Correlation Estimation by Joint Sparse
Regression Model. 
V.Q. Wu, W. Sun, L. Hsu (2020), SpaceLog: an R package for inferring gene-gene networks using SPACE model with log penalty
}

\examples{
library(spacelog)

data(BA_simu_n_400_p_100_e_1) # there are 100 different simulation run in X and A
k=1
Data = X[[k]]
TrueData=A[[k]]
n=dim(Data)[1]
p=dim(Data)[2]
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
imethod=2

# lam2=0 is LASSO penalty 

r1 = space.joint(nDat, lam1=lamb1.vec[1], lam2=0, lam3=0, weight=2,iter=2)
}


