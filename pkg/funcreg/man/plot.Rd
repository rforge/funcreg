\name{plotCoef}
\alias{plotCoef}
\alias{plotConfInt}
	
\title{Plot functions for functional coefficients}

\description{
 It plot the estimated functional coefficients with confidence
 bands. The coefficients are estimated by \code{\link{funcreg}}.
}
\usage{
plotCoef(obj, which, n = c(50, 50), fixeds=NULL, fixedt=NULL, plotWF=TRUE, ...)

plotConfInt(obj, which, type = c("beta_t", "beta_s", "beta_st"),
                       n = c(50, 50), level = 0.95, plotWF = TRUE, beta=NULL)
}

\arguments{
  \item{obj}{An object of class "funcreg"}
  
 \item{which}{Which functional coefficient we want to plot. 0 for the
 intercept, 1 for the first functional slope and so on.}

\item{n}{The number of grid points for the plot. The first element is for
  the \code{t} dimension and the second for the \code{s} dimension.}

 \item{fixeds}{If \code{NULL}, a three dimensional plot is produced.
 Alternatively, we can set a value for \code{s} and plot beta(t,s) for a
 given \code{s}. It has no effect when \code{which=0}.}

 \item{fixedt}{If \code{NULL}, a three dimensional plot is produced.
 Alternatively, we can set a value for \code{t} and plot beta(t,s) for a
 given \code{t}. It has no effect when \code{which=0}.}

 \item{plotWF}{If \code{FALSE}, the \code{wireframe} is not shown. In
 that case, only the output of the function is returned. Only applicable
 for three dimensional plot.}

 \item{type}{The type of plot (see details).}
 
 \item{level}{The level for the confidence interval}

 \item{beta}{A one or two argument function (see details).}

 \item{...}{Other argument. to pass to either \code{plot} or \code{wireframe}.}
}

\details{
If \code{type="beta_t"}, the default, the integral with respect to
\code{s} is shown with pointwise confidence bands. It is simply the functional
intercept with its pointwise confidence bands if \code{which=0}.

If \code{type="beta_s"} (not applicable for \code{which=0}), the
integral with respect to \code{t} is shown with its pointwise confidence bands.

If \code{type="beta_st"} (not applicable for \code{which=0}), the whole
three dimensional curve is shown with its pointwise confidence
surfaces.

If \code{beta} is provided, the estimated curve is replaced by the one
produced by the function \code{beta}, but the same confidence bands are
produced. Useful for pointwise hypothesis testing. If
\code{type="beta_s"} or \code{type="beta_t"}, \code{beta} must be a
function of \code{t} or \code{s}. If \code{type="beta_st"}, it must be a
function of both. 
}

\value{
If \code{plotWF=FALSE} and the expected plot is a surface, it returns an
object of class  "trellis" produced by \code{\link{wireframe}}.
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

plotCoef(res, 1, type="beta_st")

plotConfInt(res, 1, type="beta_t")

## Assume we what to do pointwise hypothesis testing
plotConfInt(res, 1, type="beta_t", beta=function(t) .6-.7*t)

## Nice 3-D confidence interval
plotConfInt(res, 1, type="beta_st")

}

