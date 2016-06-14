\name{funcreg}
\alias{funcreg}
\alias{summary.funcreg}
\alias{print.funcreg}
\alias{fitted.funcreg}

\title{Linear Functional Regression}

\description{
 Estimation of the functional parameters in a multiple functional
 regression. The smoothing parameters and the number of basis is assumed
 to be known.
}

\usage{
funcreg(form, create_basis=create.bspline.basis, LD=2, lambda,
                    k, regularized=TRUE, CstInt=FALSE, ...)

\method{print}{funcreg}(x, ...)

\method{summary}{funcreg}(object, ...)

\method{fitted}{funcreg}(object, ...)
}

\arguments{
  \item{x}{An object of class "funcreg"}
   
  \item{object}{An object of class "funcreg"}
  
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

 \item{regularized}{If \code{TRUE}, the system is regularized if it is
 badly conditioned. The value is slowly increases until the linear
 system of equation can be solved.}
	
 \item{CstInt}{If \code{TRUE}, the intercept is assumed constant. If
   \code{FALSE}, it is a one dimensional functional parameter.}

 \item{...}{Other argument that is passed to \code{create_basis}
   (see \code{\link{create.bspline.basis}})}
}

\value{
It returns an object of 'class' "funcreg" 

The object of class "funcreg" is a list containing at least:

\item{X}{The list of functional regressors.}

\item{Y}{The functional response.}

\item{data}{A list with element "t" for time units of the observed data
  and "Y" for the matrix of observed response varable.}
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


## we just set lambda and k to arbitrary values
res <- funcreg(yfd~xfd, k=c(5,5), lambda=c(.001,.001,.001))
res

summary(res)

plot(fitted(res))
}

