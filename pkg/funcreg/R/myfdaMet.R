print.myfdaCV <- function(x, digits=4, conv=FALSE, ...)
    {
        cat("Cross-Validation of Non-linear FDA Fit\n")
        cat("**************************************\n")
        if (all(dim(x$cv)==1))
            {
                cat("CV = ", round(c(x$cv), digits), "\n")
            } else {
                print.default(format(x$cv, digits=digits), quote=FALSE,
                              print.gap=2L)
            }
        if (any(!x$convergence) & !conv)
            cat("Warning: Some estimations did not converge. Print with conv=TRUE for details\n")
        if (conv)
            {
                cat("\n Proportion of Y used in the CV\n")
                cat("********************************\n")
                t <- dim(x$info)[1]
                n <- dim(x$info)[2]
                nlam <- dim(x$info)[3]
                nk <- dim(x$info)[4]
                info <- x$info
                res <- sapply(1:nk, function(i) sapply(1:nlam, function(j)
                    sum(info[,,j,i]==0)))
                res <- res/(t*n)
                res <- matrix(c(res), nlam, nk)
                dimnames(res) <- dimnames(x$cv)
                print.default(format(res, digits=digits), quote=FALSE,
                              print.gap=2L)
            }
        invisible(x)
    }

print.myfda <- function(x, ...)
    {
        cat(x$type, "\n")
        cat(ncol(x$y), " Y, ", nrow(x$y), " time units\n", sep="")
        cat(x$basis$nbasis, " number of basis and Lambda = ", x$lambda,
            "\n",sep="")
        if (!is.null(x$convergence))
            cat("Convergence: ", all(x$convergence==0), "\n",sep="")
    }


plot.myfda <- function(x, which=NULL, addpoints=TRUE, npoints=100, ...)
    {
        #        t <- seq(x$basis$rangeval[1], x$basis$rangeval[2], len=npoints)
        ranget <- range(x$t)
        t <- seq(ranget[1], ranget[2], len=npoints)
        if (is.null(which))
            {
                
                basisval <- eval.basis(t,x$basis)
                yhat <- basisval%*%x$coefficients
                yhat <- x$link(yhat)
                ylim <- range(yhat, na.rm=TRUE)
                plot(t, yhat[,1], col=1, lty=2, xlab="t", ylab="Y(t)",
                     ylim = ylim, type="l", ...)
                for (i in 2:ncol(x$y))
                    lines(t,yhat[,i], col=i, lty=2)
            } else {                
                basisval <- eval.basis(t,x$basis)
                yhat <- basisval%*%x$coefficients[,which]
                yhat <- as.matrix(yhat)
                yhat <- x$link(yhat)                
                ylim <- range(c(yhat, x$y[,which]), na.rm=TRUE)
                plot(t, yhat[,1], col=1, lty=2, xlab="t", ylab="Y(t)",
                     ylim = ylim, type="l", ...)
                if (addpoints)
                    points(x$t, x$y[,which[1]], pch=21, col=1)                
                for (i in which[-1])
                    {
                        lines(t,yhat[,i], col=i, lty=2)                        
                        if (addpoints)
                            points(x$t, x$y[,i], pch=21, col=i)                
                    }
            }
    }

