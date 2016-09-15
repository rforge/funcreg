#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/



#### Function to get the coefficients
#### Y is a matrix, not a list (txn)

coefEst <- function(y, t, lam, k, L=2,
                    create_basis=create.bspline.basis, ...)
    {
        type <- "Linear FDA Fit"
        y <- as.matrix(y)
        n <- ncol(y)
        T <- nrow(y)    
        if (length(t) != T)
            stop("nrow(y) must be equal to length(t)")
        rangeval <- range(t)
        basis <- create_basis(rangeval,k, ...)
        basisval <- eval.basis(t,basis)
        pen <- eval.penalty(basis,L)    
        res <- .Fortran("coefest", as.double(y), as.double(basisval), as.integer(n),
                        as.integer(k), as.integer(T), as.double(lam), as.double(pen),
                        coef = double(k*n), PACKAGE="funcreg")
        ans <- list(coefficients=matrix(res$coef, k, n), convergence=NULL,
                    t=t, y=y, basis=basis, link=function(x) x,
                    lambda=lam, L=L, type=type)
        class(ans) <- "myfda"
        ans
    }

#### Main function to get the CV matrix (nlam x nK)
#### lam and k are vectors
####################################################
fdaCV <- function(y, t, lamvec, kvec, L=2,
                  create_basis=create.bspline.basis, obj=NULL, ...)
    {
        if (!is.null(obj))
            {
                if (class(obj) != "myfda")
                    stop("object must be of class myfda")
                y <- obj$y; t <- obj$t; kvec <- obj$basis$nbasis; n <- ncol(y)
                L <- obj$L; T <- nrow(y); nlam <- 1; lamvec <- obj$lambda
            } else {
                y <- as.matrix(y)
                n <- ncol(y)
                T <- nrow(y)
                nlam <- length(lamvec)
                if (length(t) != T)
                    stop("nrow(y) must be equal to length(t)")
            }
        rangeval <- range(t)
        f <- function(k, basis=NULL)
            {
                if (is.null(basis))
                    basis <- create_basis(rangeval,k, ...)
                basisval <- eval.basis(t,basis)
                pen <- eval.penalty(basis,L)    
                res <- .Fortran("crval", as.double(y), as.double(basisval),
                                as.integer(n), as.integer(k), as.integer(T),
                                as.integer(nlam), as.double(lamvec),
                                as.double(pen), cv = double(nlam), PACKAGE="funcreg")
                res$cv
            }
        if (is.null(obj))
            ans <- sapply(kvec, f)
        else
            ans <- f(kvec, obj$basis)
        ans
    }

## This function finds the optimal Lambda for a given k
## using either the Brent method or the Golden Section method
## maxitalgo and tolalgo are the parameters for the minimization of the
## cross-validation. The linear Case
 #################################################
getLam <- function(y, t, lamInt, k, L=2,
                   create_basis=create.bspline.basis, 
                   maxitalgo=100, tolalgo=1e-7, type=c("Brent","GS"), ...)
    {
        type = match.arg(type)
        type = ifelse(type=="Brent", "lambrent", "lamgs")
        rangeval <- range(t)
        basis <- create_basis(rangeval,k, ...)
        basisval <- eval.basis(t,basis)
        pen <- eval.penalty(basis,L)    
        T <- nrow(y)
        n <- ncol(y)
        res <- .Fortran(type, as.double(y), as.double(basisval),
                        as.integer(n), as.integer(k), as.integer(T),
                        as.double(lamInt[1]), as.double(lamInt[2]), as.double(pen),
                        as.integer(maxitalgo), as.double(tolalgo), 
                        info = integer(1), iter = integer(1),
                        lam = double(1), cv=double(1), PACKAGE="funcreg")
        c(lam=res$lam, info=res$info, iter=res$iter, cv=res$cv)                        
    }


### This function search over kvec and
### pick the one with minimum CV. For each K, the optimal Lambda is obtained
### numerically either by the Brent Method or the Golden Section Method.
### maxitalgo and tolalgo are the parameters for the minimization of the
### cross-validation. logLamInt is a 2x1 vector with the upper and lower bound
### for the search.
### It returns an object of class myfda
### The linear case
#########################################################################
getLandKopt <- function(y, t, lamInt, kvec, L=2,
                        create_basis=create.bspline.basis, 
                        maxitalgo=100, tolalgo=1e-6, type=c("Brent", "GS"), ...)
    {
        type <- match.arg(type)
        res <- sapply(kvec, function(k) getLam(y=y, t=t, lamInt,
                                               k=k, L=L,
                                               create_basis=create_basis,
                                               maxitalgo=maxitalgo,
                                               tolalgo=tolalgo, type=type, ...))
        w <- which(res[4,]==min(res[4,],na.rm=TRUE))
        lambda <- res[1,w]
        k=kvec[w]
        obj <- coefEst(y=y, t=t, lam=lambda, k=k, L=L,
                       create_basis=create_basis, ...)
        obj$cv <- res[4,w]
        obj$algo <- list(lamInter=lamInt, convergence=res[2,], numIter=res[3,])
        obj
    }

### Function to get the minimum CV
#####################################
getLandK <- function(y, t, lamvec, kvec, L=2,
                     create_basis=create.bspline.basis, ...)
    {
        cv <- fdaCV(y, t, lamvec, kvec, L, create_basis, ...)
        w <- which(cv==min(cv,na.rm=TRUE), arr.ind=TRUE)
        lambda <- lamvec[w[1]]
        k=kvec[w[2]]
        res <- coefEst(y=y, t=t, lam=lambda, k=k, L=L,
                         create_basis=create_basis, ...)
        res$cv <- min(cv,na.rm=TRUE)
        res$grid <- list(lamgrid=lamvec, kgrid=kvec, cv=cv)
        res
    }

