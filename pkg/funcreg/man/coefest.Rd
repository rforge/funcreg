\name{coefEst}
\alias{coefEst}
\alias{nlCoefEst}
	
\title{Functional Data Fit}

\description{
 These are the main functions to fit a functional data object to a
 dataset. The model can be linear or nonlinear. The smoothing parameter
 and the number of basis is assumed to be known. 
}
\usage{
coefEst(y, t, lam, k, L=2, rangeval = c(0,1),
        create_basis=create.bspline.basis, ...)
nlCoefEst(y, t, lam, k, L=2, rangeval = c(0,1),
          create_basis=create.bspline.basis, maxit=100, tol=1e-8, ...)
}
\arguments{
\item{y}{A \eqn{T \times n} matrix of data points, when \eqn{T} is the number
  of periods and \eqn{n} the number of individuals.}

\item{t}{A numeric vector of time units. The number of elements is equal
  to the number of rows of \code{y}. It must have a range inside
  \code{rangeval}.}

\item{lam}{The regularization parameter that penalizes for the absence of
  smoothness.}

\item{k}{The number of basis.}

\item{L}{Either a nonnegative integer defining an order of a
          derivative or a linear differential operator (see
          \code{\link{eval.penalty}}).}
	
\item{rangeval}{The time interval over which the basis are defined (see
  \code{\link{create.bspline.basis}}).}

\item{create_basis}{The function used to create the basis object (see
  \code{\link{create.bspline.basis}}).}

\item{maxit}{The maximum number of iteration for the Newton method}

\item{tol}{The tolerance parameter for the stopping rule of the Newton
  method}

\item{...}{Other argument that is passed to \code{create_basis}  (see \code{\link{create.bspline.basis}})}
}

\details{
In \code{coefEst}, the functional data is a linear combination of the
basis, while in \code{nlCoefEst} it is the exponential of the linear
combination of the basis. The latter is used when we want to force the
functional data curve to be positive}

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

## Nonlinear Estimation
res <- nlCoefEst(y=GDPv56, t=t, lam=0.003, k=15)
res


}

