\name{myfdamethod}
\alias{print.myfda}
\alias{plot.myfda}
	
\title{Methods for Functional Data Fit}

\description{
 These are the plot and print methods which are associated with the objects of class "myfda".
}
\usage{
\method{plot}{myfda}(x, which=NULL, addpoints=TRUE, npoints=100, ...)

\method{print}{myfda}(x, ...)

}
\arguments{
  \item{x}{Object of class "myfda" (see \code{\link{coefEst}}).}

  \item{which}{Which curve shall we plot}

  \item{addpoints}{If \code{TRUE} and \code{which} is not null,
  observations are added to the graph}
  
\item{npoints}{The number of points used to plot the curves}

\item{...}{Other argument that is passed to the \code{plot} or
  \code{print} methods}
}

\examples{
data(GDPv56) 
t <- seq(0,1,len=nrow(GDPv56))

## Linear estimation
res <- coefEst(y=GDPv56, t=t, lam=0.003, k=15)
print(res)
plot(res)
plot(res,which=1:3)

}

