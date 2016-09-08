\name{getLandK}
\alias{nlGetLam}
\alias{getLam}
\alias{nlGetLandK}
\alias{nlGetLandKopt}
\alias{getLandK}
\alias{getLandKopt}
	
\title{Cross-Validation}

\description{Several functions that computes the leave-one-out
  cross-validation of functional data. They include cases of linear and
  nonlinear fit.}

\usage{
nlGetLam(y, t, lamInt, k, L=2, rangeval = c(0,1),
         create_basis=create.bspline.basis, maxit=100, tol=1e-8,
         maxitalgo=100, tolalgo=1e-7, type=c("Brent","GS"),
         typels=c("brent","grid"), ...)

nlGetLandKopt(y, t, lamInt, kvec, L=2, rangeval=c(0,1),
              create_basis=create.bspline.basis, maxit=100, tol=1e-7,
              maxitalgo=100, tolalgo=1e-6, type=c("Brent", "GS"),
              typels=c("brent","grid"), ...)

nlGetLandK(y, t, lamvec, kvec, L=2, rangeval=c(0,1),
           create_basis=create.bspline.basis, maxit=100, tol=1e-8,
           typels=c("brent","grid"), ...)

getLam(y, t, lamInt, k, L=2, rangeval = c(0,1),
       create_basis=create.bspline.basis, 
       maxitalgo=100, tolalgo=1e-7, type=c("Brent","GS"), ...)

getLandK(y, t, lamvec, kvec, L=2, rangeval=c(0,1),
         create_basis=create.bspline.basis, ...)

getLandKopt(y, t, lamInt, kvec, L=2, rangeval=c(0,1),
            create_basis=create.bspline.basis, 
            maxitalgo=100, tolalgo=1e-6, type=c("Brent", "GS"), ...)
}
\arguments{
\item{y}{A \eqn{T \times n} matrix of data points, when \eqn{T} is the
  number of periods and \eqn{n} the number of individuals.}

\item{t}{A numeric vector of time units. The number of elements is equal to the
number of rows of \code{y}. It must have a range inside
\code{rangeval}.}

\item{lamvec}{A vector of regularization parameters
  that penalizes for the absence of smoothness.}

\item{lamInt}{A vector that contains the lower and upper bounds for
  Lambda that is passed to the Brent or Golden Section method.}

\item{kvec}{A vector of
  integers that indicates the number of basis.}

\item{k}{A scalar that indicates the number of basis.}

\item{L}{Either a
nonnegative integer defining an order of a derivative or a linear
differential operator (see \code{\link{eval.penalty}}).}

\item{rangeval}{The time interval over which the basis are defined (see
  \code{\link{create.bspline.basis}}).}

\item{create_basis}{The function used to create the basis object (see
  \code{\link{create.bspline.basis}}).}

\item{maxit}{The maximum number
  of iteration for the Newton method}

\item{tol}{The tolerance parameter for the stopping rule of the Newton
  method}

\item{maxitalgo}{Maximum number of iteration for the brent and
  Golden-Section methods}

\item{tolalgo}{Tolerance level for the stopping rule of the Brent or
  Golden-Section method.}

\item{type}{Which algorithm should we use to get the lambda.}

\item{typels}{The \code{nlCoefEst} is based on Newton with line
  search. The argument refers to the method used for the line search}

\item{...}{Other argument that is passed to \code{create_basis}}  (see
\code{\link{create.bspline.basis}}).}

\details{Function that starts with 'nl' are the nonlinear functional
  data (see \code{\link{coefEst}}).

  \code{getLam} and \code{nlGetLam}, finds the Lambda that
  minimizes the cross-validation for a given \code{k} using either the
  Brent method or the Golden-Section method.

  \code{nlGetLandK} and \code{getLandK} find the optimal Lambda and
  \code{k} using a grid search. The grid is provided by the user through
  \code{lamvec} and \code{kvec}.

  \code{getLandKopt} and \code{nlGetLandKopt} find the optimal Lambda and
  \code{k} using a grid search for \code{k} only. Lambda is obtained
  using the Brent or Golden-Section method. 
 }

 \value{
   \code{getLam} and \code{nlGetLam} returns a list with the
   following items:
   
   \item{lam}{The optimal Lambda}
   
   \item{info}{A convergence code for the optimization method (0 for
   normal convergence).}
   
   \item{iter}{The number of iterations for the Brent or Golden-Section method.}
   
   \item{cv}{The minimum cross-validation}
   
   \code{getLandK}, code{nlGetLandK}, \code{getLandKopt}, and
   code{nlGetLandKopt}, returns an object of class "myfda" (see
   \code{\link{coefEst}}). In addition to the usual element, the list
   inludes information about the search.  

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

# Using a grid
nlam <- 5
lamvec <- 10^seq(-5,1,length.out=nlam)
kvec <- 5:7
res1 <- getLandK(GDPv56, t, lamvec, kvec)
res1

# Using Brent
res2 <- getLandKopt(GDPv56, t, lamvec, kvec, tolalgo=1e-5)
res2

## Nonlinear estimation
###########################

# Using a grid
res3 <- nlGetLandK(GDPv56, t, lamvec, kvec)
res3

# Using Brent
res4 <- nlGetLandKopt(GDPv56, t, lamvec, kvec, tolalgo=1e-5)
res3

}

