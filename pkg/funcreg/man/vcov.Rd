\name{vcov}
\alias{vcov.funcreg}
	
\title{Variance of functional coefficients}

\description{
 It computes the functional variance of the estimated functional
 parameters obtained by \code{\link{funcreg}}
}
\usage{
\method{vcov}{funcreg}(object, which, type = c("beta", "beta_t", "beta_s", "beta_st"),
                         s = NULL, t = NULL, ...)
}

\arguments{
  \item{object}{An object of class "funcreg"}
  
 \item{which}{The variance of which functional parameters shall we
   compute? 0 for the intercept, 1 for the first functional slope and so on.}

 \item{type}{What type of variance to compute (see details).}

 \item{s}{If type is "beta_s", the variance is computed at each element
   of \code{s}. If \code{NULL}, it is computed at the same time units as
   for the observed observations. }

 \item{t}{If type is "beta_t", the variance is computed at each element
   of \code{s}. If \code{NULL}, it is computed at the same time units as
   for the observed observations. }

 \item{...}{Argument to pass to other objects}
}

\details{
The function computed the variance of the integral of the functional
parameters (single for the intercept and double for the slopes) when
\code{type='beta'}.

If \code{type='beta_t'}, it is the pointwise variance of the functional
intercept if \code{which=0} and the pointwise variance of the integral
of the slopes which respect to \code{s}. By default (\code{t=NULL}) the
time units of the observed data is used. Alternatively, we can provide
a vector of points.

If \code{type='beta_s'}, (not available for \code{which=0}) it is the
pointwise variance of the integral of the slopes which respect to
\code{y}.  By default (\code{s=NULL}), equally spaced grid points are
generated using the range of s-basis and the number of time units equals
to the one from the observed observations.

If \code{type='beta_st'} (not available for \code{which=0}), a two
dimensional pointwise variance is computed.

For the last three cases, \code{vcov} returns either a one dimension
functional data object or a two-dimensional one. 
}

\value{
It returns a list containing at least:

\item{V}{Either a numeric value when \code{type='beta'} or an 'fda'
  object when \code{type='beta_t'} or \code{type='beta_s'}, or a matrix
  when \code{type='beta_st'}.}
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

vcov(res, which=0)

plot(vcov(res, 1, type="beta_t")$V)

}

