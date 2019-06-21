\name{spaceSimu}
\alias{spaceSimu}
\title{A simulated data example for inferring the partial correlation matrix}
\description{A list containing a simulated data example for package \code{space}}

\value{
spaceSimu is a list of two components: 
\item{Y.data}{a numeric matrix consisting of 500 rows and 250 columns.}
\item{ParCor.true}{a numeric matrix consisting of 500 rows and 500 columns.} 
}

\details{
\code{Y.data} is a simulated array data with columns corresponding to samples and rows corresponding genes/clones. 
\code{ParCor.true} is the true partical correlation matrix used to simulate \code{Y.data}. }

\references{
J. Peng, P. Wang, N. Zhou, J. Zhu (2007). Partial Correlation Estimation by Joint Sparse
Regression Model.
}

\usage{data(spaceSimu)}

\keyword{datasets}

