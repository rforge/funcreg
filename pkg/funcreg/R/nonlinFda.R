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
#### Y is a matrix
nlCoefEst <- function(y, t, lam, k, L=2, rangeval = c(0,1),
                      create_basis=create.bspline.basis, maxit=100, tol=1e-8, ...)
    {
        type = "Non-Linear FDA: Exp() fit"
        y <- as.matrix(y)
        n <- ncol(y)
        T <- nrow(y)    
        if (length(t) != T)
            stop("nrow(y) must be equal to length(t)")
        if (any(t>rangeval[2]) | any(t<rangeval[1]))
            stop("Values of t are not in the rangeval provided")
        basis <- create_basis(rangeval,k, ...)
        basisval <- eval.basis(t,basis)
        pen <- eval.penalty(basis,L)
        res <- .Fortran("nlcoefest", as.double(y), as.double(basisval), as.integer(n),
                        as.integer(k), as.integer(T), as.double(lam), as.double(pen),
                        as.double(tol), as.integer(maxit), info = integer(n),
                        coef = double(k*n), PACKAGE="funcreg")
        ans <- list(coefficients=matrix(res$coef, k, n), convergence=res$info,
                    t=t, y=y, basis=basis, link=function(x) exp(x),
                    lambda=lam, L=L, type=type, maxi=maxit, tol=tol)
        class(ans) <- c("myfda")
        ans
    }

### This function puts all R and basisvalue for all K
### into an array
#########################################################
.getPenArray <- function(kvec, t, L, create_basis, rangeval, ...)
    {
        pen <- array(0,dim=c(max(kvec),max(kvec), length(kvec)))
        basisval <- array(0,dim=c(length(t), max(kvec), length(kvec)))
        for (i in 1:length(kvec))
            {
                b <- create_basis(rangeval,kvec[i], ...)
                tmp <- eval.penalty(b, L)
                tmp2 <- eval.basis(t, b)
                pen[1:kvec[i], 1:kvec[i], i] <- tmp
                basisval[,1:kvec[i], i] <- tmp2
            }
        list(pen=pen, basisval=basisval)
    }

#### Main function to get the CV matrix (nlam x nK)
#### lam and k are vectors
####################################################
nlFdaCV <- function(y, t, lamvec, kvec, L=2, rangeval = c(0,1),
                    create_basis=create.bspline.basis, maxit=100, tol=1e-8,
                    obj=NULL, ...)
    {
        if (!is.null(obj))
            {
                if (class(obj) != "myfda")
                    stop("object must be of class myfda")
                y <- obj$y; t <- obj$t; kvec <- obj$basis$nbasis; n <- ncol(y)
                L <- obj$L; maxit <- obj$maxit; tol <- obj$tol
                T <- nrow(y); nlam <- 1; lamvec <- obj$lambda
                nk <- 1; maxk=kvec
                obj <- list(pen=eval.penalty(obj$basis, L),
                            basisval <- eval.basis(t, obj$basis))
            } else {        
                obj <- .getPenArray(kvec, t, L, create_basis, rangeval, ...)
                maxk <- max(kvec)
                nk <- length(kvec)
                nlam <- length(lamvec)
                T <- nrow(y)
                n <- ncol(y)
            }
        res <- .Fortran("nlcrval", as.double(y), as.double(obj$basisval),
                        as.integer(n), as.integer(kvec), as.integer(T),
                        as.integer(nlam), as.integer(nk),as.integer(maxk),
                        as.double(lamvec),  as.double(obj$pen), as.double(tol),
                        as.integer(maxit), info = integer(T*n*nlam*nk),
                        cv = double(nlam*nk), PACKAGE="funcreg")
        convergence <- array(res$info, c(T, n, nlam, nk))
        dimnames(convergence) <- list(paste("t",1:T,sep=""),
                                      paste("Y",1:n, sep=""),
                                      paste("Lambda",1:nlam,sep=""),
                                      paste("K",1:nk,sep=""))
        cv <- matrix(res$cv,nlam,nk)
        dimnames(cv) <- list(paste("Lambda",1:nlam,sep=""),
                             paste("K",1:nk,sep=""))
        chk <- .chkFct(convergence)
        ans <- list(cv=cv, convergence=chk, kvec=kvec, lamvec=lamvec,
                    info=convergence)
    }

### This function search over loglam and kvec and
### pick the one with minimum CV
### it returns an object of class myfda
#########################################################################
nlGetLandK <- function(y, t, lamvec, kvec, L=2, rangeval=c(0,1),
                       create_basis=create.bspline.basis, maxit=100, tol=1e-8, ...)
    {
        cv <- nlFdaCV(y=y, t=t, lamvec=lamvec, kvec=kvec, L=L,
                       rangeval=rangeval, create_basis=create_basis,
                       maxit=maxit, tol=tol, ...)
        chk <- cv$convergence
        if (any(chk))
            mess = "Some estimations failed to converge"
        else
            mess = "All estimations converged"        
        cv$cv[chk] <- NA
        if (all(is.na(cv$cv)))
            return("All estimations failed to converge")            
        w <- which(cv$cv==min(cv$cv,na.rm=TRUE), arr.ind=TRUE)
        lambda <- lamvec[w[1]]
        k=kvec[w[2]]
        res <- nlCoefEst(y=y, t=t, lam=lambda, k=k, L=L, rangeval = rangeval,
                         create_basis=create_basis, maxit=maxit, tol=tol, ...)
        res$cv <- min(cv$cv,na.rm=TRUE)
        res$grid <- list(lamgrid=lamvec, kgrid=kvec, cv=cv$cv, mess=mess)
        res
    }

## This function checks for any bad convergence
## while computing the cross-validation
 #################################################
.chkFct <- function(obj)
    {
        nk <- dim(obj)[4]
        nlam <- dim(obj)[3]
        chk <- sapply(1:nk, function(i) sapply(1:nlam, function(j) any(obj[,,j,i]>0)))
    }

## This function finds the optimal Lambda for a given k
## using either the Brent method or the Golden Section method
## maxit and tol are the tuning parameters for the estimation of the coefficients
## maxitalgo and tolalgo are the parameters for the minimization of the
## cross-validation. 
 #################################################
nlGetLam <- function(y, t, lamInt, k, L=2, rangeval = c(0,1),
                   create_basis=create.bspline.basis, maxit=100, tol=1e-8,
                   maxitalgo=100, tolalgo=1e-7, type=c("Brent","GS"), ...)
    {
        type <- match.arg(type)
        type <- ifelse(type=="Brent", "nllambrent", "nllamgs")
        obj <- .getPenArray(k, t, L, create_basis, rangeval, ...)
        maxk <- k
        nk <- 1
        T <- nrow(y)
        n <- ncol(y)
        res <- .Fortran(type, as.double(y), as.double(obj$basisval),
                        as.integer(n), as.integer(k), as.integer(T),
                        as.double(lamInt[1]), as.double(lamInt[2]), as.double(obj$pen),
                        as.double(tol), as.integer(maxit), as.integer(maxitalgo),
                        as.double(tolalgo), info = integer(1), iter = integer(1),
                        lam = double(1), cv=double(1), PACKAGE="funcreg")
        c(lam=res$lam, info=res$info, iter=res$iter, cv=res$cv)                        
    }

### This function search over kvec and
### pick the one with minimum CV. For each K, the optimal Lambda is obtained
### numerically either by the Brent Method or the Golden Section Method.
### maxit and tol are the tuning parameters for the estimation of the coefficients
### maxitalgo and tolalgo are the parameters for the minimization of the
### cross-validation. logLamInt is a 2x1 vector with the upper and lower bound for the search.
### It returns an object of class myfda
#########################################################################
nlGetLandKopt <- function(y, t, lamInt, kvec, L=2, rangeval=c(0,1),
                          create_basis=create.bspline.basis, maxit=100, tol=1e-7,
                          maxitalgo=100, tolalgo=1e-6, type=c("Brent", "GS"), ...)
    {
        type <- match.arg(type)
        res <- sapply(kvec, function(k) nlGetLam(y=y, t=t, lamInt=lamInt, k=k, L=L,
                                                 rangeval=rangeval,
                                                 create_basis=create_basis,
                                                 maxit=maxit, tol=tol,
                                                 maxitalgo=maxitalgo,
                                                 tolalgo=tolalgo, type=type, ...))
        w <- which(res[4,]==min(res[4,],na.rm=TRUE))
        lambda <- res[1,w]
        k=kvec[w]
        obj <- nlCoefEst(y=y, t=t, lam=lambda, k=k, L=L, rangeval = rangeval,
                         create_basis=create_basis, maxit=maxit, tol=tol, ...)
        obj$cv <- res[4,w]
        obj$algo <- list(lamInter=lamInt, convergence=res[2,], numIter=res[3,])
        obj
    }


makeLinFda <- function(obj, npoints=200)
    {
        if (class(obj) != "myfda")
            stop("Only applicable to objects of class myfda")
        type <- strtrim(obj$type, 3)
        if (type=="Lin")
            return(obj)
        rangeval <- obj$basis$rangeval
        t <- seq(rangeval[1], rangeval[2], length.out=npoints)
        xhat <- eval.basis(t, obj$basis)
        xhat <- xhat%*%obj$coefficients
        xhat <- obj$link(xhat)
        xhat <- smooth.basis(t, xhat, fdPar(obj$basis))$fd
        obj$coefficients <- xhat$coefs
        obj$link <- function(x) x
        obj$type <- "Linear FDA fitted from a Non-linear FDA"
        obj
    }
    
logFda <- function(obj)
    {
        if (class(obj) != "myfda")
            stop("Only applicable to objects of class myfda")
        if (obj$type=="Lin")
            stop("Only applicable to Non-linear FDA")        
        obj$link <- function(x) x
        obj$type <- "Linear FDA obtained by taking the log od Non-linear FDA"
        obj$y <- log(obj$y)
        obj
    }
