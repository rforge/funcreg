\name{crossval}
\alias{funcregCV}
\alias{getFuncregLam}

\title{Leave-one-out cross-validation}

\description{
  It computes the cross-validation for a given lambda and k or find the
  optimal lambda for a given k (see below).
  }

\usage{
funcregCV(form, create_basis=create.bspline.basis, LD=2, lambda,
                      k, CstInt=FALSE, obj=NULL, ...)

getFuncregLam(form, create_basis=create.bspline.basis, LD=2, lam0,
                          k, CstInt=FALSE, ..., loglam=FALSE,
                          method="BFGS", optimArg=list())
}

\arguments{
  \item{obj}{An object of class "funcreg". If \code{NULL}, all other
    arguments must be set. }

 \item{form}{A formula for the functional regression. Arguments are
    objects of class "myfda" obtained either by \code{\link{coefEst}} or
    \code{\link{nlCoefEst}}, and must be in the local environment.}

 \item{create_basis}{The function used to create the basis object (see
   \code{\link{create.bspline.basis}}) for the functional parameters.}

 \item{LD}{Either a nonnegative integer defining an order of a
   derivative or a linear differential operator (see
   \code{\link{eval.penalty}}).}

 \item{lambda}{A vector of regularization parameters that penalizes for the absence of
   smoothness. The dimention is 2 times the number of regressors plus 1
   if there is an intercept and if \code{CstInt=FALSE}. The first is for
   the intercept, the second and third for the first regressor and so
   on. If it is a scalar, it is assumed that the smoothing parameter is
   the same for all functional parameters.}

 \item{lam0}{Initial value of lambda for the minimization of the cross-validation}
 
 \item{k}{A vector of integers representing the number of basis for each
   functional parameter. The first is for the intercept (if there is an
   intercept and that \code{Cstint=FALSE}), the second for the first
   functional parameter and so on. If it is a scalar, it is assumed that
   the number of basis is the same for all functional parameters.}
  
 \item{CstInt}{If \code{TRUE}, the intercept is assumed constant. If
   \code{FALSE}, it is a one dimensional functional parameter.}

 \item{loglam}{If \code{TRUE}, lam0 is \code{log10(lambda)} and the
   function minimizes \code{CV(10^lam)}. It prevents the lambda's from
   being negative.}

 \item{method}{The algorithm used by \code{link{optim}}}

 \item{optimArg}{A list of additional arguments to pass to \code{link{optim}}}
 
 \item{...}{Additional arguments to be passed to \code{create_basis}
   (see \code{\link{create.bspline.basis}}).}

}

\value{
'funregCV' returns the value of the leave-one-out cross-validation. The
attributes 'convergence', returns convergence code for each
estimation.

'getFuncregLam' returns the fitted object of class 'funcreg' using the
optimal lambda. The argument 'optimRes' gives the \code{\link{optim}}
output from the minimization of the cross-validation. 
}
\references{
   Ramsay, James O., & Silverman, Bernard W. (2005), 
   Functional Data Analysis, \emph{Springer}, New York.
 }
 
\examples{
data(GDPv56) 

## We just create response and a covariate artificialy from the GDP
## series

y <- GDPv56[,1:30]
x <- GDPv56[,31:60]
t <- seq(0,1,len=nrow(y))

## First we create the "myfda" objects
yfd <- coefEst(y, t, .0004, 15, 2)
xfd <- coefEst(x, t, .0004, 15, 2)


## with the object
res <- funcreg(yfd~xfd, k=c(5,5), lambda=c(.001,.001,.001))
funcregCV(obj=res)

## without object
funcregCV(yfd~xfd, k=c(5,5), lambda=c(.001,.001,.001))

## With the default values
getFuncregLam(yfd~xfd, k=c(5,5), lam0=c(.001,.001,.001))

## With upper and lower bound
getFuncregLam(yfd~xfd, k=c(5,5), lam0=c(.001,.001,.001), method="L-BFGS-B",
             optimArg=list(lower=c(0,0,0),upper=c(200,200,200))) 

## With loglam
res <- getFuncregLam(yfd~xfd, k=c(5,5), lam0=c(-3,-3,-3), loglam=TRUE)
res

## Get info about the optim()
res$optimRes
}

