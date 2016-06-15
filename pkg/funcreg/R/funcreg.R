funcreg <- function(form, create_basis=create.bspline.basis, LD=2, lambda,
                    k, regularized=TRUE, CstInt=FALSE, ...)
    {
        tr <- terms(form)
        call <- match.call()
        if (attr(tr, "response") != 1)
            stop("You cannot run a functional regression without response variable")
        namey <- rownames(attr(tr, "factors"))[1]
        namex <- colnames(attr(tr, "factors"))
        all <- eval(attr(tr, "variables"))
        Yfd <- all[[1]]
        data <- list(t=Yfd$t, Y=Yfd$y)
        Yfd <- fd(Yfd$coef, Yfd$basis)
        Xfd <- all[-1]
        names(Xfd) <- namex
        Intercept <- attr(tr, "intercept")
        nx <- length(Xfd)
        for (i in 1:nx)
            Xfd[[i]] <- fd(Xfd[[i]]$coef, Xfd[[i]]$basis)
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
                betaList[[i+Intercept]] <- bifdPar(bb, LD, LD, lambda[start],
                                                   lambda[start+1])
		start <- start + 2		
            }
	chk <- TRUE
	while (chk)
            {
                resFin <- try(.funcregFit(Xfd, Yfd, betaList,
                                          Intercept=(Intercept==1),data=data,
                                          regularized=regularized),silent=TRUE) 
                if (class(resFin) == "try-error")
                    {
                        regularized = regularized*10
                    } else {
                        chk <- FALSE
                    }
            }
        ans <- list(res = resFin, lambda = lambda,
                    X = Xfd, Y = Yfd, LD=LD, regularized=regularized, data=data,
                    Intercept=Intercept, call=call)
        class(ans) <- "funcreg"
        ans
    }


### The main estimation function for multivariate functional
### regression with functional response and
###  two dimensional functional parameters
##############################################################

.funcregFit <- function (xfdobjList, yfdobj, betaList, data=NULL,
                         vcov = c("iid", "noniid"), wtvec = NULL,
                         Intercept = TRUE, regularized = TRUE, sandwich=TRUE,
                         getGCV=FALSE, multicore=0) 
    {
### data is a list with the matrix of observations named as follows:
### t (T vector) for the index of time of the observed data, and
### Y (TxN) for the matrix of observed response. It is used to compute the 
### covariance matrix of the error term. If NULL, no covariane matrix can be computed

	Nx <- length(xfdobjList)
	vcov <- match.arg(vcov)

	if ( !(Intercept %in% c(0,1)))
            stop("Intercept is either TRUE of FALSE")

### Collecting information about Y
### Checking if it is of the right format
#########################################
	ybasis = yfdobj$basis
	ynbasis = ybasis$nbasis
	ranget = ybasis$rangeval
	nfine = max(c(201, 10 * ynbasis + 1))
	tfine = seq(ranget[1], ranget[2], len = nfine)
	coefy = yfdobj$coef
	coefdy = dim(coefy)
	ncurves = coefdy[2]
	if (!is.fd(yfdobj)) 
            stop("YFD is not a functional data object.")

### Collecting information about the X's
### Checking if they are of the right format
############################################
	check <- sapply(1:Nx, function(i) !is.fd(xfdobjList[[i]]))
	if (any(check)) 
            stop("On or more XFD are not functional data objects.")
	xbasis <- list()
	ranges <- list()
	coefx <- list()
	coefdx <- list()
	for (i in 1:Nx)
            {
		xbasis[[i]] <- xfdobjList[[i]]$basis
	        ranges[[i]] = xbasis[[i]]$rangeval
		coefx[[i]] <- xfdobjList[[i]]$coef
		coefdx[[i]] <- dim(coefx[[i]])
            }
	check <- sapply(1:Nx, function(i) coefdx[[i]][2] != ncurves)
	if (any(check)) 
            stop("Numbers of observations in one or more X's and Y do not match.")

### Collecting information on the beta's
### Checking that they have the right format
#############################################

	if (!inherits(betaList, "list")) 
            stop("betaList is not a list object.")
	if (length(betaList) != (Nx+Intercept)) 
            stop("The number of beta's do not match the number of regressors.")
	nalpha <- 0
	if (Intercept)
            {
		alphafdPar = betaList[[1]]
	        alphafd = alphafdPar$fd
	        alphaLfd = alphafdPar$Lfd
	        alphalambda = alphafdPar$lambda
	        alphabasis = alphafd$basis
	        alpharange = alphabasis$rangeval
	        alphanbasis = alphabasis$nbasis
	        if (!inherits(alphafdPar, "fdPar")) 
                    stop("Alpha = BETACELL[[1]] is not a fdPar object.")
	        if (alpharange[1] != ranget[1] || alpharange[2] != ranget[2]) 
                    stop("Range of ALPHAFD coefficient and YFD not compatible.")
		nalpha <- alphanbasis
            }
	betabifdPar <- list()
	betasLfd <- list()
	betatLfd <- list()
	betaslambda = list()
	betatlambda = list()
	betabifd = list()
	betasbasis = list()
	betasrange = list()
	betatbasis = list()
	betatrange = list()
	betasnbasis = list()
	betatnbasis = list()
	nbeta <- 0

	for (i in 1:Nx)
            {
	        betabifdPar[[i]] = betaList[[Intercept+i]]
	        if (!inherits(betabifdPar[[i]], "bifdPar")) 
                    stop("One of the beta(s,t) is not a bifdPar object.")
		betasLfd[[i]] = betabifdPar[[i]]$Lfds
	        betatLfd[[i]] = betabifdPar[[i]]$Lfdt
	        betaslambda[[i]] = betabifdPar[[i]]$lambdas
	        betatlambda[[i]] = betabifdPar[[i]]$lambdat
	        betabifd[[i]] = betabifdPar[[i]]$bifd
	        betasbasis[[i]] = betabifd[[i]]$sbasis
	        betasrange[[i]] = betasbasis[[i]]$rangeval
	        betatbasis[[i]] = betabifd[[i]]$tbasis
	        betatrange[[i]] = betatbasis[[i]]$rangeval
	        betasnbasis[[i]] = betasbasis[[i]]$nbasis
	        betatnbasis[[i]] = betatbasis[[i]]$nbasis
	        if (betasrange[[i]][1] != ranges[[i]][1] || betasrange[[i]][2] != ranges[[i]][2]) 
                    stop("Range of one of the BETASFD coefficient not compatible with its XFD.")
	        if (betatrange[[i]][1] != ranget[1] || betatrange[[i]][2] != ranget[2]) 
                    stop("Range of one of the BETATFD coefficient and YFD not compatible.")
		nbeta <- nbeta+betasnbasis[[i]]*betatnbasis[[i]]
            }
### Starting to compute lists of matrices
#########################################

### The list of penalty matrices
##############################
	betassmat <- list()
	betattmat <- list()
	RList <- list()
	if(Intercept)
            {
		if (alphalambda>0)
                    RList[[1]] <- eval.penalty(alphabasis, alphaLfd)
            }

	for (i in 1:Nx)
            {
	        betassmat[[i]] = inprod(betasbasis[[i]], betasbasis[[i]])
	        betattmat[[i]] = inprod(betatbasis[[i]], betatbasis[[i]])
		RList[[i+Intercept]] <- list(s=0,t=0)
		if(betaslambda[[i]]>0)
                    RList[[i+Intercept]]$s <- eval.penalty(betasbasis[[i]], betasLfd[[i]])
		if(betatlambda[[i]]>0)
                    RList[[i+Intercept]]$t <- eval.penalty(betatbasis[[i]], betatLfd[[i]])
            }

### Constructing F, G and H inprod
#####################################
	if (Intercept)
            {
	        Finprod = inprod(ybasis, alphabasis)
	        Fmat = crossprod(coefy, Finprod)
            }

	Ginprod <- list()
	Hinprod <- list()
	Gmat <- list()
	Hmat <- list()
	H1CP <- list()
	HHCP = list()

	for (i in 1:Nx)
            {
	        Ginprod[[i]] = inprod(ybasis, betatbasis[[i]])
	        Hinprod[[i]] = inprod(xbasis[[i]], betasbasis[[i]])
	        Gmat[[i]] = crossprod(coefy,Ginprod[[i]])
	        Hmat[[i]] = crossprod(coefx[[i]], Hinprod[[i]])
		H1CP[[i]] = as.matrix(colSums(Hmat[[i]]))	
	        HHCP[[i]] = crossprod(Hmat[[i]])
            }

### Building the Cmat matrix and the Dmat
##########################################

	.Matfct <- function(Sigma=NULL, isPen=TRUE, both=FALSE)
            {
                Dmat <- vector()
                if (Intercept)
                    Dmat <- c(Dmat, colSums(Fmat))
                for (i in 1: Nx)
                    {
                        HGCP = crossprod(Gmat[[i]], Hmat[[i]])
                        Dmat <- c(Dmat, c(HGCP))
                    }
                Cmat <- matrix(0, ncoef, ncoef)
                if (both)
                    Cmat_noPen <- matrix(0, ncoef, ncoef)
                if (Intercept)
                    {
                        if (is.null(Sigma)) {
                            Cst <- ncurves
                        } else {
                            if (is.null(dim(Sigma)))
				Cst <- ncurves*Sigma
                            else
				Cst <- sum(Sigma)
                        }
                        Cmat[1:alphanbasis, 1:alphanbasis] <- Cst *
                            inprod(alphabasis, alphabasis)
                        if (both)
                            Cmat_noPen[1:alphanbasis, 1:alphanbasis] <-
                                Cmat[1:alphanbasis, 1:alphanbasis]
                        if (alphalambda > 0 & isPen)
                            Cmat[1:alphanbasis, 1:alphanbasis] <-
                                Cmat[1:alphanbasis, 1:alphanbasis] +
				alphalambda*RList[[1]]
                        start <- alphanbasis+1
                        for (i in 1: Nx)
                            {
                                nb <- betatnbasis[[i]]*betasnbasis[[i]]
                                betalttmat = inprod(betatbasis[[i]], alphabasis)
                                if (!is.null(Sigma))
                                    {
                                        if (is.null(dim(Sigma)))
                                            {
                                                H1CP2 <- Sigma*H1CP[[i]]
                                            } else {
                                                H1CP2 <- crossprod(Hmat[[i]], colSums(Sigma))
                                            }
                                    } else {
                                        H1CP2 <- H1CP[[i]]	
                                    }

                                Cmat[1:alphanbasis, start:(start+nb-1)] <-
                                    t(kronecker(H1CP2, betalttmat))
                                Cmat[start:(start+nb-1), 1:alphanbasis] <-
                                    t(Cmat[1:alphanbasis, start:(start+nb-1)])
                                if (both)
                                    {
                                        Cmat_noPen[1:alphanbasis, start:(start+nb-1)] <-
                                            Cmat[1:alphanbasis, start:(start+nb-1)]
                                        Cmat_noPen[start:(start+nb-1), 1:alphanbasis] <-
                                            Cmat[start:(start+nb-1), 1:alphanbasis]
                                    }
                                start <- start+nb
                            }
                        colstart <- alphanbasis+1
                        rowstart <- alphanbasis+1
                    } else {
                        colstart <- 1
                        rowstart <- 1
                    }

                for (i in 1:Nx)
                    {
                        if (!is.null(Sigma))
                            {
				if (is.null(dim(Sigma)))
                                    HHCP2 <- Sigma*HHCP[[i]]
				else
                                    HHCP2 <- crossprod(Hmat[[i]],Sigma%*%Hmat[[i]])
                            } else {
                                HHCP2 <- HHCP[[i]]	
                            }
                        tempMat = kronecker(HHCP2, betattmat[[i]])
                        rowind <- rowstart:(rowstart+nrow(tempMat)-1)
                        colind <- colstart:(colstart+ncol(tempMat)-1)
                        rowstart <- rowstart+nrow(tempMat)
                        colstart <- colstart+ncol(tempMat)
                        if (both)
                            Cmat_noPen[rowind,colind] <- tempMat
                        if (betaslambda[[i]] > 0 & isPen) 
                            tempMat = tempMat + betaslambda[[i]] *
                                kronecker(RList[[Intercept+i]]$s, betattmat[[i]])
                        if (betatlambda[[i]] > 0 & isPen) 
                            tempMat = tempMat + betatlambda[[i]] *
                                kronecker(betassmat[[i]], RList[[Intercept+i]]$t)
                        Cmat[rowind,colind] <- tempMat
                        if (i < Nx)
                            {
                                for (j in (i+1):(Nx))
                                    {
                                        if (!is.null(Sigma))
                                            {
                                                if (is.null(dim(Sigma)))
                                                    HHCP2 <- Sigma*crossprod(Hmat[[j]],
                                                                             Hmat[[i]])
                                                else
                                                    HHCP2 <- crossprod(Hmat[[j]],
                                                                       Sigma%*%Hmat[[i]])
                                            } else {
						HHCP2 <- crossprod(Hmat[[j]],Hmat[[i]])
                                            }
                                        betaij <- inprod(betatbasis[[j]],betatbasis[[i]])
                                        tempMat <- kronecker(HHCP2, betaij)
                                        colind <- (max(colind)+1):(max(colind)+nrow(tempMat))
                                        Cmat[rowind,colind] <- t(tempMat)
                                        Cmat[colind,rowind] <- tempMat
                                        if (both)
                                            {
                                                Cmat_noPen[rowind,colind] <- t(tempMat)
                                                Cmat_noPen[colind,rowind] <- tempMat
                                            }
                                    }
                                
                            }

                    }
                if (both)
                    {	
                        return(list(Cmat=Cmat, Dmat=Dmat,Cmat_noPen=Cmat_noPen))
                    } else {
                        return(list(Cmat=Cmat, Dmat=Dmat))
                    }
            }
	ncoef <- nalpha+nbeta
	if (getGCV) {
            AllMat <- .Matfct(both=TRUE)
            Cmat_noPen <- AllMat$Cmat_noPen
	} else {
            AllMat <- .Matfct()}
	Cmat <- AllMat$Cmat
	Dmat <- AllMat$Dma

### Solving the system
### and building the coef matrices
### If eigen values are too small, we regularized the solution	
	if (regularized == TRUE)
            {
                eigCmat <- eigen(Cmat)$val 
                if (min(eigCmat)<1e-6)
                    Cmat <- Cmat + diag(ncol(Cmat))*(0.1*min(eigCmat[eigCmat>1e-5]))
            }
	if(is.numeric(regularized))
            {
 		Cmat <- Cmat + diag(ncol(Cmat))*regularized
            }	
	
	coefvec = symsolve(Cmat, Dmat)

	if (Intercept)
            {
	    	alphacoef = coefvec[1:alphanbasis]
	    	alphafdnames = yfdobj$fdnames
	    	alphafdnames[[3]] = "Intercept"
	    	alphafd = fd(alphacoef, alphabasis, alphafdnames)
		start <- alphanbasis+1
	        yhatmat = eval.fd(tfine, alphafd) %*% matrix(1, 1, ncurves) 
            } else {
		start <- 1
            }
	betacoef <- list()
	betafd <- list()
	xbetacoef <- list()
	xbetafd <- list()
	for (i in 1:Nx)
            {
		nb <- betatnbasis[[i]]*betasnbasis[[i]]
		ind <- start:(start+nb-1)
		start <- start+nb	
		betacoef[[i]] <- matrix(coefvec[ind], betatnbasis[[i]], betasnbasis[[i]])
		betafdnames = xfdobjList[[i]]$fdnames
	        betafdnames[[3]] = paste("Reg. Coefficient_", i, sep="")
	        betafd[[i]] = bifd(t(betacoef[[i]]), betasbasis[[i]],
                          betatbasis[[i]], betafdnames)
	        xbetacoef[[i]] = betacoef[[i]] %*% t(Hmat[[i]])
	        xbetafd[[i]] = fd(xbetacoef[[i]], betatbasis[[i]])
		if ((i==1) & !Intercept)	
                    yhatmat <- eval.fd(tfine, xbetafd[[i]])
		else  
                    yhatmat <- yhatmat + eval.fd(tfine, xbetafd[[i]])  
            }
	yhatfd = smooth.basis(tfine, yhatmat, ybasis)$fd
	linmodList = list(betaestbifdList = betafd, yhatfdobj = yhatfd,
            xbetafd = xbetafd, Hj = Hinprod, Cmat = Cmat)
	if (Intercept)
            linmodList$beta0estfd <- alphafd

### Computing variances
####################################
	if (!is.null(data))
            {
		yhat <- eval.fd(data$t, yhatfd)
		e <- data$Y-yhat
		if (vcov != "iid")
                    stop("non iid not yet implemented")
                sigma <-  var(c(e))
		if(sandwich)
                    {
			meat <- .Matfct(isPen=FALSE, Sigma=sigma)$Cmat
		        Vcoef <- solve(Cmat, meat)
		        Vcoef <- solve(Cmat, t(Vcoef))
                    } else {
			Vcoef <- sigma*solve(Cmat)	
                    }
	        Sig <- list()
        	if (Intercept)
                    {
			Sig[[1]] <- Vcoef[1:alphanbasis,1:alphanbasis]
			start <- alphanbasis+1
                    } else {
			start <- 1
                    }	
		for (i in 1:Nx)
                    {
			ind <- start:(start+betasnbasis[[i]]*betatnbasis[[i]]-1)	
			Sig[[i+Intercept]] <- Vcoef[ind, ind]
			start <- start+betasnbasis[[i]]*betatnbasis[[i]]
                    }
		linmodList$covPar <- Sig
		linmodList$data <- data
            }

### Computing the generalized cross-validation
##############################################
	if (getGCV)
            {
### warning("The computation of the GCV is experimental for now")
		Integrate_n <- 30
		e <- yfdobj-yhatfd
		ssr <- sum(inprod(e^2))
		P <- function(t)
                    {
			e <- c(eval.fd(t,yfdobj-yhatfd))
			ssr <- sum(e^2)	
			Wmat_t <- vector()
			if (Intercept)
                            Wmat_t <- rep(1,ncurves)%x%eval.basis(t,alphabasis)					
			for (i in 1:Nx)
                            {
				phi_t <- eval.basis(t,betatbasis[[i]])
				Wmat_t <- cbind(Wmat_t, phi_t%x%Hmat[[i]])
                            }
			Pmat <- Wmat_t%*%symsolve(Cmat, t(Wmat_t))
			df <- sum(diag(Pmat))
			GCV <- ncurves*ssr/(ncurves-df)^2
                    }
		Simpson <- function(f, a, b, n, ...) 
                    {
			n <- floor(n/2) * 2 + 1
			x <- seq(a, b, len = n)
			z <- rep(c(4, 2), (n - 3)/2)
			z <- c(1, z, 4, 1)
			h <- x[2] - x[1]
			sum(z * f(x, ...) * h/3)
                    }
		Pf <- function(t)
                    {
                        res <- lapply(t,P)
			simplify2array(res)
                    }
### Have to think about it...                 
### GCV <- Simpson(Pf, yfdobj$basis$rangeval[1],yfdobj$basis$rangeval[2], Integrate_n)
### df <- Simpson(Pf, yfdobj$basis$rangeval[1],yfdobj$basis$rangeval[2], Integrate_n)
                
		df <- sum(diag(symsolve(Cmat, Cmat_noPen)))
		n <- ncurves
		GCV <- n*ssr/(n-df)^2
		attr(GCV,"n") <- n
		attr(GCV,"df") <- df
		attr(GCV,"ssr") <- ssr
		linmodList$GCV <- GCV
            }
	return(linmodList)
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
        for (i in 1:length(object$res$betaestbifdList))
            {
                cname <- c(cname, names(object$X))
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
        invisible(x)
    }

fitted.funcreg <- function(object, ...)
    {
        object$res$yhatfdobj
    }

residuals.funcreg <- function(object, ...)
    {
        e <- object$Y-object$res$yhatfdobj
        e
    }
