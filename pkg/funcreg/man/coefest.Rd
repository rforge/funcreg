\name{coefEst}
\alias{coefEst}
\alias{nlCoefEst}
\alias{makeLinFda}
\alias{logFda}
\alias{fitted.myfda}
\alias{as.fd.myfda}

\title{Functional Data Fit}

\description{
 These are the main functions to fit a functional data object to a
 dataset. The model can be linear or nonlinear. The smoothing parameter
 and the number of basis is assumed to be known. 
}
\usage{
coefEst(y, t, lam, k, L=2,
        create_basis=create.bspline.basis, ...)

nlCoefEst(y, t, lam, k, L=2,
          create_basis=create.bspline.basis, maxit=100, tol=1e-8,
          typels=c("brent","grid"), addDat=FALSE, rangeval=NULL, ...)

makeLinFda(obj, npoints=200)

logFda(obj)

\method{as.fd}{myfda}(x, fdnames=NULL, npoints=200, ...)

\method{fitted}{myfda}(object, t=NULL, ...)
}
\arguments{
\item{y}{A \eqn{T \times n} matrix of data points, when \eqn{T} is the number
  of periods and \eqn{n} the number of individuals.}

\item{t}{A numeric vector of time units. In \code{coefEst} and
  \code{nlCoefEst}, the number of elements is equal
  to the number of rows of \code{y}. For the \code{fitted} method, it
  can be a fine grid, and the range must be inside the range of the basis.}

\item{lam}{The regularization parameter that penalizes for the absence of
  smoothness.}

\item{k}{The number of basis.}

\item{L}{Either a nonnegative integer defining an order of a
          derivative or a linear differential operator (see
          \code{\link{eval.penalty}}).}
	
\item{create_basis}{The function used to create the basis object (see
  \code{\link{create.bspline.basis}}).}

\item{maxit}{The maximum number of iteration for the Newton method}

\item{tol}{The tolerance parameter for the stopping rule of the Newton
  method}

\item{...}{Other argument that is passed to \code{create_basis}  (see
  \code{\link{create.bspline.basis}}), or to other methods in the case
  of \code{as.fd.myfda}.}

\item{obj}{Object of class "myfda"}

\item{npoints}{Number of points to fit the linear FDA to the non-linear
  one}

\item{x}{Object of class "myfda"}

\item{object}{Object of class "myfda"}

\item{fdnames}{Names of elements of the "fd" object (see
  \code{\link{fd}})}

\item{typels}{The \code{nlCoefEst} is based on Newton with line
  search. The argument refers to the method used for the line search}

\item{addDat}{If TRUE, fake observations are added before and after
  assuming stationarity. for the estimation, the time span is expanded
  accordingly.}

\item{rangeval}{This can be set if the basis is defined on a range other
  than \code{range(t)}.} 
}

\details{
In \code{coefEst}, the functional data is a linear combination of the
basis, while in \code{nlCoefEst} it is the exponential of the linear
combination of the basis. The latter is used when we want to force the
functional data curve to be positive,

\code{makeLinFda} will transform the nonlinear fit to a linear FDA by
fitting the nonlinear curve to a linear combination of the same
basis. The curve is almost identical and it is transformed this way
by \code{\link{funcreg}}

\code{logFDA} transforms the non-linear FDA into a linear FDA. It only
redefine the \code{link} equation and take the log of the observed data.
}

\value{
It returns an object of 'class' "myfda" 

The object of class "myfda" is a list containing at least:

\item{coefficients}{\eqn{k\times n} matrix of coefficients.}

\item{y}{The dataset matrix \code{y}.}

\item{t}{The input vector \code{t}.}

\item{basis}{The basis object}

\item{link}{The function that links \code{y} to the linear combination
  of the basis.}

\item{lambda}{The value of \code{lam}.}

\item{convergence}{Information on the convergence of the Newton method.}

 }

 \references{
   Ramsay, James O., & Silverman, Bernard W. (2005), 
   Functional Data Analysis, \emph{Springer}, New York.
 }
 
\examples{
data(GDPv56) 
t <- seq(0,1,len=nrow(GDPv56))

## Linear estimation
res <- coefEst(y=GDPv56, t=t, lam=0.003, k=15)
res

cbind(data=GDPv56[,1], fitted=fitted(res)[,1])

## Nonlinear Estimation
res <- nlCoefEst(y=GDPv56, t=t, lam=0.003, k=15)
res

resLin <- makeLinFda(res)
resLin

## Take the logarithm of the non-linear functional data
logres <- logFda(res)
plot(res,1:2)

## Transform the myfda object into an fd object
resfd <- as.fd(resLin)

## Nonlinear Estimation with zeros
data(simData)
t <- seq(0,1,length.out=5)
res <- nlCoefEst(y=simData, t=t, lam=0.003, k=7)
res
res <- nlCoefEst(y=simData, t=t, lam=0.003, k=7, addDat=TRUE)
res


}

