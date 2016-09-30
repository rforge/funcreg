funcreg <- function(form, create_basis=create.bspline.basis, LD=2, lambda,
                     k, regularized=TRUE, CstInt=FALSE, CV=FALSE, ...)
    {
        tr <- terms(form)
        call <- match.call()
        if (attr(tr, "response") != 1)
            stop("You cannot run a functional regression without response variable")
        namey <- rownames(attr(tr, "factors"))[1]
        namex <- colnames(attr(tr, "factors"))
        all <- eval(attr(tr, "variables"),  )
        Yfd <- all[[1]]
        if (strtrim(Yfd$type, 3)=="Non")
            Yfd <- makeLinFda(Yfd)
        data <- list(t=Yfd$t, Y=Yfd$y)
        if (!is.null(colnames(Yfd$y)))
            {
                tmp <- colnames(Yfd$y)                
            } else {
                tmp <- paste(namey, "_", 1:ncol(Yfd$y), sep="")
            }
        Yfd <- as.fd(Yfd, fdnames=list("time", tmp, namey))
        Xfd <- all[-1]
        names(Xfd) <- namex
        Intercept <- attr(tr, "intercept")
        nx <- length(Xfd)
        for (i in 1:nx)
            {
                if (!is.null(colnames(Xfd[[i]]$y)))
                    tmp <- colnames(Xfd[[i]]$y)
                else 
                    tmp <- paste(namex[i], "_", 1:ncol(Xfd[[i]]$y), sep="")
                Xfd[[i]] <- as.fd(Xfd[[i]], fnames=list("time", tmp, namex[i]))
            }
        rangeval <- Yfd$basis$rangeval        
        chk <- sapply(1:nx, function(i) all(Xfd[[i]]$basis$rangeval != rangeval))
        if (any(chk))
            stop("The rangeval of some regressors is not like the one from the response")
        nk0 <- !CstInt*Intercept
        nk <- nx + nk0
        nlam <- 2*nx + nk0
        if (length(k) == 1)
            {
                k <- rep(k,nk)
            } else {
                if (nk != length(k))
                    stop("Wrong dimension of k")
            }                
        if (length(lambda) == 1)
            {
                warning("Only one Lambda provided. it is assumed they are all constant")
                lambda <- rep(lambda, nlam)
            } else {                
                if ((2*nx+nk0) != length(lambda))                    
                    stop("Wrong number of smoothing parameters")
            }
        betaList <- list()
        if (Intercept)
            {
                if (CstInt)
                    {
                        b <- create.constant.basis(rangeval)
                        betaList[[1]] <- fdPar(b, LD)                
                    } else {
                        b <- create_basis(rangeval, k[1], ...)
                        betaList[[1]] <- fdPar(b, LD, lambda[1])
                    }
            }
        start <- nk0+1  
        for (i in 1:nx)
            {
                b <- create_basis(rangeval, k[i+nk0], ...)
                bb <- bifd(matrix(0, k[i+nk0], k[i+nk0]), b, b)
                betaList[[i+Intercept]] <- bifdPar(bb, LD, LD, lambda[start+1],
                                                   lambda[start])
		start <- start + 2		
            }
	chk <- TRUE
        obj <- .prepMat(Xfd, Yfd, betaList)
        obj$data <- data
        resFin <- .funcregFit(obj)
        if (CV)
            {
                res2 <- .funcregCV(obj)
                resFin$CV <- res2$CV
                resFin$cvInfo <- res2$cvInfo
            }
        ans <- list(res = resFin, lambda = lambda,
                    X = Xfd, Y = Yfd, LD=LD, data=data,
                    Intercept=Intercept, call=call,obj=obj)        
        class(ans) <- "funcreg"
        ans
    }


### Computes all what is needed for the Fortran code.
### All of it is independent on the lambda's
####################################################
.prepMat <- function (xfdobjList, yfdobj, betaList)    
    {
        nx <- length(xfdobjList)
        n <- ncol(yfdobj$coefs)
        intercept <- nx == (length(betaList)-1)
        if (intercept)
            {
                k0 = betaList[[1]]$fd$basis$nbasis
                lam0 <- betaList[[1]]$lambda
                R0 <- eval.penalty(betaList[[1]]$fd$basis, betaList[[1]]$Lfd)
                ipaa <- inprod(betaList[[1]]$fd$basis, betaList[[1]]$fd$basis)
                ipya <- inprod(yfdobj$basis, betaList[[1]]$fd$basis)
                Fmat <- crossprod(yfdobj$coefs, ipya)
            }else {
                k0 <- 0; lam0 <- 0
                R0 <- 0; ipaa <- 0
                ipya <- 0
                Fmat <- 0
            }
        kb <- vector()
        kx <- vector()
        ky <- yfdobj$basis$nbasis
        lam <- vector()
        if (nx == 0)
            return(list(k0=k0, lam0=lam0, R0=R0, ipaa=ipaa, ipya=ipya, ky=ky,
                        kx=kx))
        for (i in 1:nx)
            {
                kx <- c(kx, xfdobjList[[i]]$basis$nbasis)
                kb <- rbind(kb, c(betaList[[i+intercept]]$bifd$tbasis$nbasis,
                                  betaList[[i+intercept]]$bifd$sbasis$nbasis))
                lam <- rbind(lam, c(betaList[[i+intercept]]$lambdat,
                                    betaList[[i+intercept]]$lambdas))
            }
        if (intercept)
            {
                ipabt <- array(0, c(k0, max(kb[,2]), nx))
                
            } else {
                ipabt <- 0

            }
        iplbs <- ipbs <- array(0, c(max(kb[,2]), max(kb[,2]), nx))
        iplbt <- ipbt <- array(0, c(max(kb[,1]), max(kb[,1]), nx))
        H <- array(0, c(n, max(kb[,2]), nx))
        ## ir order for say 3 X's
        ## b1b1', b2b1', b3b1', b2b2', b3b2', b3b3'
        ipbbt <- array(0, c(max(kb[,1]), max(kb[,1]), nx*(nx+1)/2))
        bx <- array(0, c(max(kx), n, nx))
        Rs <- Rt <- array(0, c(max(kb[,1])*max(kb[,2]),
                               max(kb[,1])*max(kb[,2]), nx))
        G <- array(0, c(n, max(kb[,1]),nx))
        by <- yfdobj$coefs
        l <- 1
        for (i in 1:nx)
            {
                if (intercept)
                    {
                        ipabt[,1:kb[i,1],i] <- inprod(betaList[[1]]$fd$basis,
                                                      betaList[[i+1]]$bifd$tbasis)
                    }
                for (j in i:nx)
                    {
                        ipbbt[1:kb[j,1], 1:kb[i,1], l] <-
                            inprod(betaList[[intercept+j]]$bifd$tbasis,
                                   betaList[[intercept+i]]$bifd$tbasis)
                        l <- l+1
                    }
                iplbs[1:kb[i,2],1:kb[i,2],i] <-
                    eval.penalty(betaList[[intercept+i]]$bifd$sbasis,
                                 betaList[[intercept+i]]$Lfds)
                iplbt[1:kb[i,1],1:kb[i,1],i] <-
                    eval.penalty(betaList[[intercept+i]]$bifd$tbasis,
                                 betaList[[intercept+i]]$Lfdt)
                ipbs[1:kb[i,2],1:kb[i,2],i] <-
                    inprod(betaList[[intercept+i]]$bifd$sbasis,
                           betaList[[intercept+i]]$bifd$sbasis)
                ipbt[1:kb[i,1],1:kb[i,1],i] <-
                    inprod(betaList[[intercept+i]]$bifd$tbasis,
                           betaList[[intercept+i]]$bifd$tbasis)
                Bi <- inprod(xfdobjList[[i]]$basis,
                             betaList[[intercept+i]]$bifd$sbasis)
                H[,1:kb[i,2],i] <- crossprod(xfdobjList[[i]]$coefs, Bi)
                bx[1:kx[i],,i] <- xfdobjList[[i]]$coefs
                nR <- kb[i,1]*kb[i,2]
                Rs[1:nR,1:nR,i] <- kronecker(iplbs[1:kb[i,2],1:kb[i,2],i],
                                              ipbt[1:kb[i,1],1:kb[i,1],i])
                Rt[1:nR,1:nR,i] <- kronecker(ipbs[1:kb[i,2],1:kb[i,2],i],
                                             iplbt[1:kb[i,1],1:kb[i,1],i])
                Bi <- inprod(yfdobj$basis,
                             betaList[[intercept+i]]$bifd$tbasis)
                G[,1:max(kb[,1]),i] <- crossprod(yfdobj$coefs, Bi)
            }
        ipyy <- inprod(yfdobj$basis, yfdobj$basis)
        ipyy <- sapply(1:n, function(i) crossprod(yfdobj$coefs[,i],
                                                  ipyy%*%yfdobj$coefs[,i]))
        list(ipya=ipya, Rt=Rt, Rs=Rs, R0=R0, H=H, G=G, ipaa=ipaa, ipabt=ipabt,
             ipbbt=ipbbt, lam=lam, lam0=lam0, k0=k0, kb=kb, kx=kx, ky=ky,Fmat=Fmat,
             bx=bx, by=by, n=n, nx=nx, ncoef=k0+sum(kb[,1]*kb[,2]),ipyy=ipyy,
             xfdobjList=xfdobjList, yfdobj=yfdobj, betaList=betaList)
    }

### The main estimation function for multivariate functional
### regression with functional response and
###  two dimensional functional parameters
##############################################################


.funcregFit <- function (obj)
    {
        res <- .Fortran("funcregwv",as.double(obj$G),as.double(obj$H),as.double(obj$R0),
                        as.double(obj$Rt),as.double(obj$Rs),as.double(obj$ipaa),
                        as.double(obj$ipabt),as.double(obj$ipbbt),
                        as.double(obj$Fmat),as.integer(obj$n),
                        as.integer(obj$nx),as.integer(obj$k0),as.integer(obj$kb),
                        as.integer(max(obj$kb[,1])),as.integer(max(obj$kb[,2])),
                        as.integer(obj$ncoef), as.double(obj$lam0),
                        as.double(obj$lam), coef=double(obj$ncoef), info=integer(1),
                        cmat=double(obj$ncoef^2), dmat=double(obj$ncoef),
                        cmat0=double(obj$ncoef^2))
        coefvec <- res$coef
        info <- res$info
        Intercept <- obj$k0>0
        # to modify to get more flexibility
        tfine <- seq(obj$yfdobj$basis$rangeval[1],
                     obj$yfdobj$basis$rangeval[2], length.out=300)
        if (Intercept)
            {
	    	alphacoef = coefvec[1:obj$k0]
	    	alphafdnames = obj$yfdobj$fdnames
	    	alphafdnames[[3]] = "Intercept"
	    	alphafd = fd(alphacoef, obj$betaList[[1]]$fd$basis, alphafdnames)
		start <- obj$k0+1
	        yhatmat = eval.fd(tfine, alphafd) %*% matrix(1, 1, obj$n) 
            } else {
		start <- 1
            }
	betacoef <- list()
	betafd <- list()
	xbetacoef <- list()
	xbetafd <- list()
	for (i in 1:obj$nx)
            {
		nb <- prod(obj$kb[i,])
		ind <- start:(start+nb-1)
		start <- start+nb	
		betacoef[[i]] <- matrix(coefvec[ind], obj$kb[i,1], obj$kb[i,2])
		betafdnames = obj$xfdobjList[[i]]$fdnames
	        betafdnames[[3]] = paste("Reg. Coefficient_", i, sep="")
	        betafd[[i]] = bifd(t(betacoef[[i]]),
                          obj$betaList[[i+Intercept]]$bifd$sbasis,
                          obj$betaList[[i+Intercept]]$bifd$tbasis, betafdnames)
                H = obj$H[,1:obj$kb[i,2],i]
	        xbetacoef[[i]] = betacoef[[i]] %*% t(H)
	        xbetafd[[i]] = fd(xbetacoef[[i]], obj$betaList[[Intercept+i]]$bifd$tbasis)
		if ((i==1) & !Intercept)	
                    yhatmat <- eval.fd(tfine, xbetafd[[i]])
		else  
                    yhatmat <- yhatmat + eval.fd(tfine, xbetafd[[i]])  
            }
	yhatfd = smooth.basis(tfine, yhatmat, obj$yfdobj$basis)$fd
        yhatfd$fdnames <- obj$yfdobj$fdnames
        yhatfd$fdnames[[3]] <- "fitted"
	linmodList = list(betaestbifdList = betafd, yhatfdobj = yhatfd,
            xbetafd = xbetafd, info=info)
	if (Intercept)
            linmodList$beta0estfd <- alphafd
        resid <- obj$yfdobj-yhatfd
        resid$fdnames <- obj$yfdobj$fdnames
        resid$fdnames[[3]] <- "residuals"
        linmodList$residuals <- resid

        ## Compute the Variance
        yhat <- eval.fd(obj$data$t, yhatfd)
        e <- obj$data$Y-yhat
        ## Need to think about a more general covariance matrix
        ## The one computed is the sandwich matrix under homoscedasticity
        Cmat <-  matrix(res$cmat, ncol=obj$ncoef)
        meat <-  matrix(res$cmat0, ncol=obj$ncoef)
        ### Only the lower triangular part is referenced in Fortran
        Cmat[upper.tri(Cmat)] <- t(Cmat)[upper.tri(Cmat)]
        meat[upper.tri(meat)] <- t(meat)[upper.tri(meat)]
        sigma <-  var(c(e))
        Vcoef <- solve(Cmat, meat)
        Vcoef <- solve(Cmat, t(Vcoef))
        Sig <- list()
        if (Intercept)
            {
                Sig[[1]] <- Vcoef[1:obj$k0,1:obj$k0]
                start <- obj$k0+1
            } else {
                start <- 1
            }	
        for (i in 1:obj$nx)
            {
                ind <- start:(start+prod(obj$kb[i,])-1)	
                Sig[[i+Intercept]] <- Vcoef[ind, ind]
                start <- start+prod(obj$kb[i,])
            }
        linmodList$covPar <- Sig
        linmodList$data <- obj$data
        return(linmodList)
    }

.funcregCV <- function (obj)
    {
        res <- .Fortran("funcregcv",as.double(obj$G),as.double(obj$H),as.double(obj$R0),
                        as.double(obj$Rt),as.double(obj$Rs),as.double(obj$ipaa),
                        as.double(obj$ipabt),as.double(obj$ipbbt),
                        as.double(obj$Fmat),as.integer(obj$n),
                        as.integer(obj$nx),as.integer(obj$k0),as.integer(obj$kb),
                        as.integer(max(obj$kb[,1])),as.integer(max(obj$kb[,2])),
                        as.integer(obj$ncoef), as.double(obj$lam0),
                        as.double(obj$lam), as.double(obj$ipyy),
                        cv=double(1), info=integer(obj$n))
        list(CV = res$cv, cvInfo = res$info)
    }

.dfuncregCV <- function (obj)
    {
        res <- .Fortran("dfuncregcv",as.double(obj$G),as.double(obj$H),as.double(obj$R0),
                        as.double(obj$Rt),as.double(obj$Rs),as.double(obj$ipaa),
                        as.double(obj$ipabt),as.double(obj$ipbbt),
                        as.double(obj$Fmat),as.integer(obj$n),
                        as.integer(obj$nx),as.integer(obj$k0),as.integer(obj$kb),
                        as.integer(max(obj$kb[,1])),as.integer(max(obj$kb[,2])),
                        as.integer(obj$ncoef), as.double(obj$lam0),
                        as.double(obj$lam), as.double(obj$ipyy),
                        dcv=double(1+2*obj$nx))
        res$dcv
    }


funcregCV <- function(form, create_basis=create.bspline.basis, LD=2, lambda,
                      k, CstInt=FALSE, obj=NULL, ...)
    {
        if (!is.null(obj))
            {
                if (class(obj) != "funcreg")
                    stop("obj must be of class funcreg")
                if (!is.null(obj$res$CV))
                    res <- obj$res[c("CV","cvInfo")]
                else
                    res <- .funcregCV(obj$obj)
            } else {
                obj <- funcreg(form=form, create_basis=create_basis, LD=LD,
                               lambda=lambda, k=k, CstInt=CstInt, CV=TRUE, ...)
                res <- obj$res[c("CV","cvInfo")]
            }
        ans <- res$CV
        attr(ans, "convergence") <- res$cvInfo
        ans
    }





vcov.funcreg <- function(object, which, type = c("beta", "beta_t", "beta_s", "beta_st"),
                         s = NULL, t = NULL, ...)
    {
        obj <- object$res
	type <- match.arg(type)
	Intercept <- !is.null(obj$beta0estfd)
	if (is.null(obj$data))
            stop("obj was created without data. Variances are therefore not available")
        
	if (which == 0)
            {
		if (!(type %in% c("beta","beta_t")))
                    stop("Only type beta_t or beta are available for beta0")
            }
	if (which==0 & !Intercept)
            stop("There is no intercept in the model")
	if (which > (Intercept + length(obj$betaestbifdList)))
            stop("which is greater than the number of estimated coefficients")
        
	if (which == 0)	
            {
		if (type == "beta_t")
                    {
			if (is.null(t))
                            t <- obj$data$t
			phit <- eval.basis(t, obj$beta0estfd$basis) 
			V <- sapply(1:length(t), function(i) t(phit[i,])%*%obj$covPar[[1]]%*%
                                        phit[i,])
			Vfd <- smooth.basis(t, c(V), obj$beta0estfd$basis)$fd
			fdnames <- list("t","","Variance of beta0(t)")
			Vfd$fdnames <- fdnames 
			return(list(V = Vfd))
                    }
		if (type == "beta")
                    {
			Temp <- inprod(obj$beta0estfd$basis,obj$beta0estfd$basis)
			V <- crossprod(c(t(Temp)), c(obj$covPar[[1]]))
			return(list(V=c(V)))
                    }			
            }
	else 
            {
		tbasis <- obj$betaestbifdList[[which]]$tbasis
		sbasis <- obj$betaestbifdList[[which]]$sbasis
		if (type == "beta_t")
                    {
			if (is.null(t))
                            t <- obj$data$t
			phit <- eval.basis(t, tbasis)
			phis <- c(inprod(sbasis))
			Temp <- kronecker(t(phis), phit)
			V <- sapply(1:length(t), function(i) t(Temp[i,]) %*%
                                        obj$covPar[[Intercept+which]]%*%Temp[i,])
			Vfd <- smooth.basis(t, c(V), tbasis)$fd 
			nvar <- paste("Variance of beta", which, "(t)",sep="")
			fdnames <- list("t","",nvar)
			Vfd$fdnames <- fdnames 
			return(list(V = Vfd))
                    }
		if (type == "beta_s")
                    {
			if (is.null(s))
                            s <- seq(sbasis$rangeval[1], sbasis$rangeval[2],
                                     length.out=length(obj$data$t))
			phis <- eval.basis(s, sbasis)
			phit <- c(inprod(tbasis))
			Temp <- kronecker(phis, t(phit))
			V <- sapply(1:length(s), function(i) t(Temp[i,]) %*%
                                        obj$covPar[[Intercept+which]]%*%Temp[i,]) 
			Vfd <- smooth.basis(s, c(V), sbasis)$fd
			nvar <- paste("Variance of beta", which, "(s)",sep="")
			fdnames <- list("s","",nvar)
			Vfd$fdnames <- fdnames 
			return(list(V = Vfd))
                    }
		if (type == "beta")
                    {
			phis <- c(inprod(sbasis))
			phit <- c(inprod(tbasis))
			Temp <- c(kronecker(t(phis), t(phit)))
			V <- t(Temp)%*%obj$covPar[[Intercept+which]]%*%Temp
			return(list(V=c(V)))
                    }
		if (type == "beta_st")
                    {
			if (is.null(s))
                            s <- seq(sbasis$rangeval[1], sbasis$rangeval[2],
                                     length.out=length(obj$data$t))
			if (is.null(t))			
                            t <- obj$data$t
			phis <- eval.basis(s, sbasis)
			phit <- eval.basis(t, tbasis)
			f <- function(i)
                            {
				Temp <- kronecker(phis,t(phit[i,]))
				V <- sapply(1:length(s), function(j) t(Temp[j,]) %*%
                                                obj$covPar[[Intercept+which]]%*%Temp[j,])
				c(V)	
                            }
			V <- sapply(1:length(t), f)
			return(list(V=V, s=s, t=t))
                    }

            }

    }

.meanFuncreg <- function(obj, which, type = c("beta", "beta_t", "beta_s"),
                         s = NULL, t = NULL)
    {
        res <- obj
        obj <- obj$res
	type <- match.arg(type)
	Intercept <- !is.null(obj$beta0estfd)
	if (which==0 & !Intercept)
            stop("There is no intercept in the model")
	if (which > (Intercept + length(obj$betaestbifdList)))
            stop("which is greater than the number of estimated coefficients")
	if (which == 0)
            {
		rng <- obj$yhatfdobj$basis$rangeval
		mbeta <- t(c(inprod(obj$beta0estfd$basis)))%*%obj$beta0estfd$coefs
		V <- vcov(res, which, type, s, t)$V
		return(list(mean=c(mbeta), sd=sqrt(V)))
            } else {
		tbasis <- obj$betaestbifdList[[which]]$tbasis
		sbasis <- obj$betaestbifdList[[which]]$sbasis
		rngt <- tbasis$rangeval
		rngs <- sbasis$rangeval
		if (type == "beta")
                    {
			Temp <- inprod(sbasis)%*%t(inprod(tbasis))
			mbeta <- crossprod(c(t(Temp)), c(obj$betaestbifdList[[which]]$coefs))
			V <- vcov(res, which, type, s, t)$V
			return(list(mean=c(mbeta), sd=sqrt(V)))
                    }
		if (type == "beta_s")
                    {
			Temp <- c(inprod(tbasis))
			if (is.null(s))
                            s <- sbasis$rangeval[1]:sbasis$rangeval[2]
			Temp <- obj$betaestbifdList[[which]]$coefs%*%Temp
			mbeta <- fd(c(Temp), sbasis)
			fdn <- paste("Beta",which, "(s)", sep="")
			mbeta$fdnames <- list("s","",fdn)
			V <- vcov(res, which, type, s, t)$V
			return(list(mean=mbeta, Vfd=V))
                    }
		if (type == "beta_t")
                    {
			Temp <- c(inprod(sbasis))
			if (is.null(t))
                            t <- tbasis$rangeval[1]:tbasis$rangeval[2]
			Temp <- crossprod(Temp,obj$betaestbifdList[[which]]$coefs)
			mbeta <- fd(c(Temp), tbasis)
			fdn <- paste("Beta",which, "(t)", sep="")
			mbeta$fdnames <- list("t","",fdn)
			V <- vcov(res, which, type, s, t)$V
			V <- fd(V$coefs, V$basis, V$fdnames)
			return(list(mean=mbeta, Vfd=V))
                    }
            }
    }

summary.funcreg <- function(object, ...)
    {        
        meanVec <- vector()
        sdVec <- vector()
        tVec <- vector()
        pval <- vector()
        cname <- vector()
        Intercept <- !is.null(object$res$beta0estfd)
        if (Intercept)
            {
                cname <- c(cname, c("Intercept"))
                rb0 <- .meanFuncreg(object, 0)
                meanVec <- c(meanVec, rb0$mean)
                sdVec <- c(sdVec, rb0$sd)
                tVec <- c(tVec, meanVec[1]/sdVec[1])
                pval <- c(pval, 2*pnorm(-abs(tVec[1]))) 
            }
        cname <- c(cname, names(object$X))
        for (i in 1:length(object$res$betaestbifdList))
            {
                rb1 <- .meanFuncreg(object, i)
                meanVec <- c(meanVec, rb1$mean)
                sdVec <- c(sdVec, rb1$sd)
                tVec <- c(tVec, meanVec[Intercept+i]/sdVec[Intercept+i])
                pval <- c(pval, 2*pnorm(-abs(tVec[Intercept+i])))
            }
        ans <- cbind(meanVec, sdVec, tVec, pval)
        dimnames(ans) <- list(cname, c("Means","Stdev","t-ratio","P-val"))
        ans
    }

plotConfInt <- function(obj, which, type = c("beta_t", "beta_s", "beta_st"),
                       n = c(50, 50), level = 0.95, plotWF = TRUE, beta=NULL)
    {
        if (class(obj) != "funcreg")
            stop("obj must be an object of class funcreg")
	z <- abs(qnorm((1-level)/2))
	type <- match.arg(type)
	Intercept <- !is.null(obj$res$beta0estfd)
	if (is.null(obj$data))
            stop("obj was created without data. Variances are therefore not available")
	if (which==0 & !Intercept)
            stop("There is no intercept in the model")
	if (which > (Intercept + length(obj$res$betaestbifdList)))
            stop("which is greater than the number of estimated coefficients")
	if (which == 0)
            {
		ranget <- range(obj$res$data$t)
		t <- seq(ranget[1],ranget[2], len=n[1])
		Vfd <- vcov(obj, 0, type = c("beta_t"), t = t)$V
		sdVec <- c(sqrt(eval.fd(t, Vfd)))
		b0 <- c(eval.fd(t, obj$res$beta0estfd))
		b0_up <- b0+z*sdVec
		b0_down <- b0-z*sdVec
		if (!is.null(beta))
                    {
			Tb0 <- beta(t)
			ylim <- range(c(Tb0, b0_up, b0_down))
			plot(t, Tb0, lwd=2, type="l",
                             main="Pointwise confident interval	for the Intercept \nwith true coefficient",
                             ylim=ylim, xlab="t",ylab=expression(hat(beta)[0](t)))	
                    } else {
                        ylim <- range(c(b0, b0_up, b0_down))
                        plot(t, b0, lwd=2, type="l",
                             main="Pointwise confident interval	for the Intercept", 
                             ylim=ylim, xlab="t",ylab=expression(hat(beta)[0](t)))
                    }
		lines(t, b0_up, lty=3) 
		lines(t, b0_down, lty=3) 
		abline(h=0)
            }
	else
            {
		ranges <- obj$res$betaestbifdList[[which]]$sbasis$rangeval
		s <- seq(ranges[1],ranges[2], len=n[2])
		ranget <- range(obj$data$t)
		t <- seq(ranget[1],ranget[2], len=n[1])
		if (type == "beta_s")
                    {
			res <- .meanFuncreg(obj, which, type, s = s, t = t)
			bj <- c(eval.fd(s, res$mean))
			sdVec <- c(sqrt(eval.fd(s, res$Vfd)))
			bj_up <- bj+z*sdVec
			bj_down <- bj-z*sdVec
			if (!is.null(beta))
                            {
				Tb <- beta(s)	
				ylim <- range(c(Tb, bj_up, bj_down))
				plot(s, Tb, lwd=2, type="l", ylim=ylim, xlab="s",
                                     ylab=expression(hat(beta)(s)))
				title(paste("Pointwise confident interval for the beta",
                                            which, "(s)\n with the true curve",sep=""))
                            } else {
				ylim <- range(c(bj, bj_up, bj_down))
				plot(s, bj, lwd=2, type="l", ylim=ylim, xlab="s",
                                     ylab=expression(hat(beta)(s)))
				title(paste("Pointwise confident interval for the beta",
                                            which, "(s)",sep=""))
                            }
			lines(s, bj_up, lty=3) 
			lines(s, bj_down, lty=3) 
			abline(h=0)
                    }
		if (type == "beta_t")
                    {
			res <- .meanFuncreg(obj, which, type, s = s, t = t)
			bj <- c(eval.fd(t, res$mean))
			sdVec <- c(sqrt(eval.fd(t, res$Vfd)))
			bj_up <- bj+z*sdVec
			bj_down <- bj-z*sdVec
			if (!is.null(beta))
                            {
                                Tb <- beta(t)	
                                ylim <- range(c(Tb, bj_up, bj_down))
                                plot(t, Tb, lwd=2, type="l", ylim=ylim, xlab="s",
                                     ylab=expression(hat(beta)(s)))
                                title(paste("Pointwise confident interval for the beta",
                                            which, "(t)\n with the true curve",sep=""))
                            } else {
				ylim <- range(c(bj, bj_up, bj_down))
				plot(t, bj, lwd=2, type="l", ylim=ylim, xlab="t",
                                     ylab=expression(hat(beta)(t)))
				title(paste("Pointwise confident interval for the beta",
                                            which, "(t)",sep=""))
                            }
			lines(t, bj_up, lty=3) 
			lines(t, bj_down, lty=3) 
			abline(h=0)
                    }
		if (type == "beta_st")
                    {
                        V <- vcov(obj, which, type, s =s, t = t)$V
                        B <- eval.bifd(s, t, obj$res$betaestbifdList[[which]])
                        z1 <- c(B)
                        zup <- z1+z*c(sqrt(V))
                        zd <- z1-z*c(sqrt(V))
                        if (!is.null(beta))
                            {
                                st <- expand.grid(s,t)
                                z1 <- sapply(1:nrow(st), function(i) beta(st[i,1],
                                                                          st[i,2]))
                                wfTitle <- paste("Confidence surface for beta", 
                                                 which,"(s,t)\nwith the true curve",
                                                 sep="")
                            } else {
                                wfTitle <- paste("Confidence surface for beta",
                                                 which,"(s,t)",sep="")
                            }
                        f <- expand.grid(s=s,t=t,gr=c(1,2,3))
                        f$z <- c(z1, zup, zd)
                        P <- wireframe(z~s*t, data=f, groups=f$gr,row.values = f$s,
                                       column.values = f$t, scales = list(arrows = FALSE), 
                                       drape = TRUE, colorkey = TRUE, main=wfTitle)
                        if (plotWF)
                            print(P)		
                        else
                            return(P)
                    }
            }
    }

plotCoef <- function(obj, which, n = c(50, 50), fixeds=NULL, fixedt=NULL,
                     plotWF=TRUE, ...)
    {
        if (class(obj) != "funcreg")
            stop("obj must be an object of class funcreg")
	if (which == 0)
            {
		if (is.null(obj$res$beta0estfd))
                    stop("There is no Intercept")
		plot(obj$res$beta0estfd)
            }
	if(which > 0)
            {
		if (is.null(obj$res$betaestbifdList))
                    {
			if (which>1)
                            stop("There is only 1 beta1")
			betabifd <- obj$res$beta1estbifd
                    } else {
			if (which>length(obj$res$betaestbifdList))
                            stop("which is greater than the number of beta1")
			betabifd <- obj$res$betaestbifdList[[which]]
                    }
		rt <- betabifd$tbasis$rangeval
		rs <- betabifd$sbasis$rangeval
		if (is.null(fixeds) & is.null(fixedt))
                    {
			t <- seq(rt[1],rt[2],len=n[1])
			s <- seq(rt[1],rt[2],len=n[2])
			B <- eval.bifd(s, t, betabifd)
			P <- wireframe(B, row.values = s, column.values = t,
                                       scales = list(arrows = FALSE), 
                                       xlab="s", ylab="t",
                                       zlab= as.expression(substitute(beta[x](s,t),
                                           list(x=which))), 
                                       drape = TRUE, colorkey = TRUE,...)
			if (plotWF)
                            print(P)
			else
                            return(list(Graph = P, betaMat = c(B)))
                    }
		else if (!is.null(fixeds))
                    {
			t <- seq(rt[1],rt[2],len=n[1])
			s <- fixeds
			B <- eval.bifd(s, t, betabifd)
			ylim1 <- ifelse(range(B)[1]<0, range(B)[1]*1.2, range(B)[1]*.8)
			ylim2 <- ifelse(range(B)[2]<0, range(B)[2]*.8, range(B)[2]*1.2)
			plot(t, c(B),xlab="t",ylab=expression(beta(bar(s),t)),
                             main=paste("bivariate coefficient with s fixed to ",
                                 fixeds,sep=""),type="l", ylim=c(ylim1,ylim2),...)
                    }
		else
                    {
			s <- seq(rt[1],rt[2],len=n[1])
			t <- fixedt
			B <- eval.bifd(s, t, betabifd)
			ylim1 <- ifelse(range(B)[1]<0, range(B)[1]*1.2, range(B)[1]*.8)
			ylim2 <- ifelse(range(B)[2]<0, range(B)[2]*.8, range(B)[2]*1.2)
			plot(s, c(B),xlab="t",ylab=expression(beta(s,bar(t))),
                             main=paste("bivariate coefficient with t fixed to ",
                                 fixedt,sep=""),type="l", ylim=c(ylim1,ylim2),...)
                    }

            }
    }

print.funcreg <- function(x, ...)
    {
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        cat("Mean of the functional parameters:\n")
        m <- summary(x)[,1,drop=FALSE]
        print.default(t(m))
        if (!is.null(x$res$CV))
            cat("Cross-Validation: ", format(x$res$CV, nsmall=5), "\n", sep="")
        invisible(x)
    }

fitted.funcreg <- function(object, ...)
    object$res$yhatfdobj

residuals.funcreg <- function(object, ...)
    object$res$residuals


as.fd.myfda <- function(x, fdnames=NULL, npoints=200, ...)
    {
        if (strtrim(x$type,3) == "Non")
            x <- makeLinFda(x, npoints=npoints)
        xfd <- fd(x$coef, x$basis, fdnames=fdnames)
        xfd
    }


getFuncregLam <- function(form, create_basis=create.bspline.basis, LD=2, lam0,
                          k, CstInt=FALSE, ..., loglam=FALSE,
                          method="BFGS", optimArg=list())
    {
        if(loglam)
            lamtmp <- 10^lam0
        else
            lamtmp <- lam0
        obj <- funcreg(form, create_basis=create_basis, LD=LD, lambda=lamtmp,
                           k=k, CstInt=CstInt, ...)$obj
        obj$loglam <- loglam
        f <- function(lam, obj)
            {
                if (obj$loglam)
                    lam <- 10^lam
                if (obj$k0>0)
                    {
                        obj$lam0 <- lam[1]
                        obj$lam = matrix(lam[-1], ncol=2)
                    } else {
                        obj$lam0 <- 0
                        obj$lam <- matrix(lam,ncol=2)
                    }
                cv <- .funcregCV(obj)$CV
                cv
            }
        df <- function(lam, obj)
            {
                if (obj$loglam)
                    lam <- 10^lam
                if (obj$k0>0)
                    {
                        obj$lam0 <- lam[1]
                        obj$lam = matrix(lam[-1], ncol=2)
                    } else {
                        obj$lam0 <- 0
                        obj$lam <- matrix(lam,ncol=2)
                    }
                dcv <- .dfuncregCV(obj)
                if (obj$k0 == 0)
                    dcv <- dcv[-1]
                if (obj$loglam)
                    dcv <- dcv*lam*log(10)
                dcv
            }
        optimArg <- c(optimArg, list(method=method, par=lam0, fn=f, gr=df,
                                     obj=obj))
        res <- do.call(optim, optimArg)
        if (loglam)
            lambda <- 10^res$par
        else
            lambda <- res$par
        res2 <- funcreg(form, create_basis=create_basis, LD=LD, lambda=lambda,
                        k=k, CstInt=CstInt, CV=TRUE, ...)
        res2$optimRes <- res
        res2
    }
