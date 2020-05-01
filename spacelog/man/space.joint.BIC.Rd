\name{space.joint.BIC}
\alias{space.joint.BIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to estimate partial correlations using the Joint Sparse Regression Model with LASSO, Elastic Net or Log penalty and use BIC to select the tuning parameter}
\description{A function to estimate partial correlations using the Joint Sparse Regression Model with LASSO, Elastic Net or Log penalty and use BIC to select the tuning parameter}

\usage{
space.joint.BIC(Data, imethod, lamb1.vec, lam2, lamb3.vec, seed)
}
\arguments{
  \item{Data}{numeric matrix. Columns are for variables and rows are for samples. 
                     Missing values are not allowed. It's recommended to first standardize each column to have mean 0 and \eqn{l_2} norm 1.}
 \item{imethod}{numeric value 1-4. It is the same weights used in Peng et al., 2007. imethod=1 will use space.neighbor function which is using neighborhood selection approach (single gene regression). imethod=2 is space_no, which means all regressions are weighted equally. imethod=3 is space_res, which uses residual variances, and imethod=4 is space_df, which uses the estimated degree of each variable.} 
    \item{lamb1.vec}{a vector, numeric value. It's the tuning parameters for penalized regression. It uses BIC to select the tuning parameter. } 

    \item{lam2}{a vector, numeric value. If lam2=0, lasso regression is used in the Joint Sparse Regression Model (JSRM). Otherwise, elastic net regression is used in JSRM and \var{lam2} serves as the \eqn{l_2} norm penalty parameter.} 

    \item{lamb3.vec}{a vector, numeric value. If not specified lam3=0, it's the same as space.joint function in space package for Joint Sparse Regression Model (JSRM) using lasso or elastic net penalty. Otherwise, log penalty is used in JSRM_log and \var{lam3} serves as the \eqn{l_0} norm penalty parameter.
                }             
    \item{seed}{numeric vector. This is for random seed}           
}  
\details{
    \code{space.joint.BIC} adapted log penalty in partial correlations under the high-dimension-low-sample-size setting (Wu, Sun, and Hsu, 2020). It's an extension to the space paper (Peng and et al., 2007), which uses LASSO or Elastic Net penalty. It use BIC to select tuning parameters. 
    }

\value{
    A list with two components, adapt the same name using output from Peng et al., 2007
    \item{ParCor}{the estimated partial correlation matrix.}
    \item{sig.fit}{numeric vector of the estimated diagonal \eqn{\sigma^{ii}}.} 
       }

\references{
J. Peng, P. Wang, N. Zhou, J. Zhu (2007), Partial Correlation Estimation by Joint Sparse
Regression Model. 
V.Q. Wu, W. Sun, L. Hsu (2020), SpaceLog: an R package for inferring gene-gene networks using SPACE model with log penalty
}



\examples{
#############################################################################################
# (A) The simulated example in Wu et al.
#############################################################################################
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

r1 = space.joint.BIC(nDat, imethod, lamb1.vec, lam2=0, lamb3.vec, seed)

}

