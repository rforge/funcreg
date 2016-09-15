\name{fdaCV}
\alias{fdaCV}
\alias{nlFdaCV}
	
\title{Cross-Validation}

\description{ The functions compute the leave-one-out cross-validation of functional
 data. It applies also to already estimated models.}
\usage{ 
nlFdaCV(y, t, lamvec, kvec, L=2,
        create_basis=create.bspline.basis, maxit=100, tol=1e-8,
         obj=NULL, typels=c("brent","grid"), addDat=FALSE, ...)

fdaCV(y, t, lamvec, kvec, L=2,
      create_basis=create.bspline.basis, obj=NULL, ...)
}
\arguments{
\item{y}{A \eqn{T \times n} matrix of data points, when \eqn{T} is the
  number of periods and \eqn{n} the number of individuals.}

\item{t}{A numeric vector of time units. The number of elements is equal to the
number of rows of \code{y}. It must have a range inside
\code{rangeval}.}

\item{lamvec}{A vector of regularization parameters
  that penalizes for the absence of smoothness.}

\item{kvec}{A vector of
  integers that indicates the number of basis.}

\item{L}{Either a
nonnegative integer defining an order of a derivative or a linear
differential operator (see \code{\link{eval.penalty}}).}

\item{create_basis}{The function used to create the basis object (see
  \code{\link{create.bspline.basis}}).}

\item{maxit}{The maximum number
  of iteration for the Newton method}

\item{tol}{The tolerance parameter for the stopping rule of the Newton
  method}

\item{obj}{Object of class "myfda" if we want the cross-validation of an already fitted model}

\item{typels}{The \code{nlCoefEst} is based on Newton with line
  search. The argument refers to the method used for the line search}

\item{addDat}{If TRUE, fake observations are added before and after
  assuming stationarity. for the estimation, the time span is expanded accordingly.}

\item{...}{Other argument that is passed to \code{create_basis} (see \code{\link{create.bspline.basis}}).}
}

\details{
It returns a matrix of cross-validations. The rows are for the lambda,
and the columns for the number of basis.
 }

\value{\code{fdaCV} returns a matrix of cross-validations. The rows are for the lambda,
and the columns for the number of basis.

\code{nlFdaCV} returns a list with the following items:

\item{cv}{The cross-validation matrix (see details)}

\item{convergence}{A matrix returning \code{TRUE} if there was no
  convergence problem while computing the cross-validation.}

\item{kvec}{The vector of integers that indicates the number of basis.}

\item{lamvec}{The vector of regularization parameters.}

\item{info}{An array of convergence codes for all estimations in the
  process of computing the cross-validations}

  }
 \references{
   Ramsay, James O., & Silverman, Bernard W. (2005), 
   Functional Data Analysis, \emph{Springer}, New York.
 }
 
\examples{
data(GDPv56) 
t <- seq(0,1,len=nrow(GDPv56))

## Linear estimation
####################

# From an object
res <- coefEst(y=GDPv56, t=t, lam=0.003, k=15)
fdaCV(obj=res)

# CV with the estimation
cv <- fdaCV(y=GDPv56,t,10^seq(-5,1,len=4), k=5:7)
cv


## Nonlinear Estimation
#######################

res <- nlCoefEst(y=GDPv56, t=t, lam=0.003, k=10)
nlFdaCV(obj=res)

# CV with the estimation
cv <- nlFdaCV(y=GDPv56,t,10^seq(-5,1,len=4), k=5:7)

## Nonlinear Estimation with zeros
##################################

data(simData)
t <- seq(0,1,length.out=5)
cv <- nlFdaCV(y=simData, t=t, lam=c(1e-5, 1e-2), k=c(6:8), addDat=TRUE)
print(cv, conv=TRUE)

}

