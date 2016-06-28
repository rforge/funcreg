\name{crossval}
\alias{funcregCV}
	
\title{Leave-one-out cross-validation}

\description{
  It computes the cross-validation for a given lambda and k (see below).
  }

\usage{
funcregCV(form, create_basis=create.bspline.basis, LD=2, lambda,
                      k, CstInt=FALSE, obj=NULL, ...)
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

 \item{k}{A vector of integers representing the number of basis for each
   functional parameter. The first is for the intercept (if there is an
   intercept and that \code{Cstint=FALSE}), the second for the first
   functional parameter and so on. If it is a scalar, it is assumed that
   the number of basis is the same for all functional parameters.}
  
 \item{CstInt}{If \code{TRUE}, the intercept is assumed constant. If
   \code{FALSE}, it is a one dimensional functional parameter.}

 \item{...}{Additional arguments to be passed to \code{create_basis}
   (see \code{\link{create.bspline.basis}}).}

}

\value{
It returns the value of the leave-one-out cross-validation. The
attributes 'convergence', returns convergence code for each
estimation. 
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


}

