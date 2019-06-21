\name{space.neighbor}
\alias{space.neighbor}
\title{A function to estimate partial correlations using the neighborhood selection approach}
\description{A function to estimate partial correlations using the neighborhood selection approach}
   
\usage{
space.neighbor(Y.m, lam1, lam2=0)
}

\arguments{
    \item{Y.m}{numeric matrix. Each column is for one variable and each row is for one sample. 
                     Missing values are not allowed. It's recommended to first standardize each column to have mean 0 and norm 1. 
                     }
    \item{lam1}{numeric value. This is the \eqn{l_1} norm penalty parameeter. If the columns of Y.m have norm one, 
                then the suggested range of lam1 is: \eqn{O(n^{1/2}\Phi^{-1}(1-\alpha/(2p^2)))} for small \eqn{\alpha} such as 0.1.
                } 
    \item{lam2}{numeric value. If not specified, lasso regression is used in the neighborhood selection.
                Otherwise, elastic net regression is used and \var{lam2} serves as the \eqn{l_2} norm penalty parameter.
                } 
 }        
\details{
    \code{space.neighbor} estimate partial correlations using the neighborhood selection approach (Meinshausen and Buhlmann, 2006).}

\value{
    A list with two components
    \item{ParCor}{the estimated partial correlation matrix.}
    \item{sig.fit}{numeric vector of the estimated \eqn{\sigma^{ii}}} 
       }

\references{
J. Peng, P. Wang, N. Zhou, J. Zhu (2007). Partial Correlation Estimation by Joint Sparse
Regression Model.


Meinshausen, N., and Buhlmann, P. (2006), High Dimensional Graphs and Variable Selection with the Lasso,  Annals of Statistics, 34, 1436-1462.
}

\author{J. Peng, P. Wang, N. Zhou, J. Zhu}

\keyword{methods}

\examples{
#############################################################################################
############################ (A) The simulated Hub.net example in Peng et. al. (2007).
#############################################################################################
data(spaceSimu)

n=nrow(spaceSimu$Y.data)
p=ncol(spaceSimu$Y.data)
true.adj=abs(spaceSimu$ParCor.true)>1e-6

#################### view the network corresponding to the parcial correlation matrix in the simulation example
#################### view the network corresponding to the parcial correlation matrix in the simulation example
########### the following code can run only if the "igraph" is installed in the system.
#library(igraph)
#plot.adj=true.adj
#diag(plot.adj)=0
#temp=graph.adjacency(adjmatrix=plot.adj, mode="undirected")
#temp.degree=apply(plot.adj, 2, sum)
#V(temp)$color=(temp.degree>9)+3
#plot(temp, vertex.size=3, vertex.frame.color="white",layout=layout.fruchterman.reingold, vertex.label=NA, edge.color=grey(0.5))


#################### estimate the parcial correlation matrix with various methods
alpha=1
l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
iter=3


########### the values of lam1 were selected to make the results of different methods comparable. 
#### 1. MB method
result1=space.neighbor(spaceSimu$Y.data, lam1=l1*0.7, lam2=0)
fit.adj=abs(result1$ParCor)>1e-6
sum(fit.adj==1)/2                  ##total number of edges detected      
sum(fit.adj[true.adj==1]==1)/2  ##total number of true edges detected        
 
#### 2. Joint method with no weight
result2=space.joint(spaceSimu$Y.data, lam1=l1*n*1.56, lam2=0, iter=iter)
fit.adj=abs(result2$ParCor)>1e-6
sum(fit.adj==1)/2                  ##total number of edges detected      
sum(fit.adj[true.adj==1]==1)/2  ##total number of true edges detected        

#### 3. Joint method with residue variance based weights
result3=space.joint(spaceSimu$Y.data, lam1=l1*n*1.86, lam2=0, weight=1, iter=iter)
fit.adj=abs(result3$ParCor)>1e-6
sum(fit.adj==1)/2                  ##total number of edges detected      
sum(fit.adj[true.adj==1]==1)/2  ##total number of true edges detected        

#### 4. Joint method with degree based weights
result4=space.joint(spaceSimu$Y.data, lam1=l1*n*1.61, lam2=0, weight=2, iter=iter)
fit.adj=abs(result4$ParCor)>1e-6
sum(fit.adj==1)/2                  ##total number of edges detected      
sum(fit.adj[true.adj==1]==1)/2  ##total number of true edges detected        

}