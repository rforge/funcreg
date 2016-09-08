c     Finding the initial value of C for the positive fda estimation
c     The initial value is the starting point for the Newton Algorithm
c     *****************************************************************            
      subroutine initialc(y, basisval, n, k, t, c0)
      integer n, t, k, i
      double precision y(t,n), basisval(t,k), c0(k,n), s

      s = sum(basisval)/t
      do i=1,n
         c0(:,i) = log(sum(y(:,i))/t)/s
      end do

      end

c     Compute the outer product of phi_i because it does not depend on the coef
c     It is therefore computed only once for all y and lambda
c     Only the upper triangle part of phiphi(:,:,i) are referenced
c     *************************************************************************            
      subroutine prepopt(basisval, k, t, phiphi)
      integer k, t, i
      double precision basisval(t,k), phiphi(k,k,t)

      phiphi = 0.0d0
      do i=1,t
         call dsyr('U', k, 1.0d0, basisval(i,:), 1, phiphi(:,:,i), k)
      end do
      end
               
c     Hessian and Jacobian for a single Y:
c     It requires phiphi from subroutine prepopt and lpen = 2*lambda*R. 
c     Only the uppper triangular part of the hessian is referenced
c     *****************************************************************      
      subroutine hesjac(y, basisval, k, t, coef, lpen, phiphi, hes,
     *     jac, pssr)
      integer t, k, i, j, l
      double precision y(t), basisval(t,k), coef(k), lpen(k,k)
      double precision cphi(t), ei, phiphi(k,k,t)
      double precision hes(k,k), jac(k), yi1, yi2, pssr

      hes = lpen
      call dgemv('N', t, k, 1.0d0, basisval, t, coef, 1, 0.0d0, cphi, 1)
      call dsymv('U', k, 1.0d0, lpen, k, coef, 1, 0.0d0, jac, 1)

      pssr = sum(coef*jac)/2.0d0
      
      do i=1,t
         ei =  exp(cphi(i))
         yi1 = -2*(y(i)-2*ei)*ei
         yi2 = -2*(y(i)-ei)*ei
         pssr = pssr + (y(i)-ei)**2
         jac = jac + yi2*basisval(i,:)
         do j=1,k
            do l=1,j
               hes(l,j)=hes(l,j)+yi1*phiphi(l,j,i)
            end do
         end do

      end do
      end

      subroutine getpssr(y, basisval, k, t, coef, lpen, phiphi, pssr)
      integer t, k, i, j, l
      double precision y(t), basisval(t,k), coef(k), lpen(k,k)
      double precision cphi(t), ei, phiphi(k,k,t)
      double precision tmp(k), pssr

      call dgemv('N', t, k, 1.0d0, basisval, t, coef, 1, 0.0d0, cphi, 1)
      call dsymv('U', k, 1.0d0, lpen, k, coef, 1, 0.0d0, tmp, 1)

      pssr = sum(coef*tmp)/2.0d0
      
      do i=1,t
         ei =  exp(cphi(i))
         pssr = pssr + (y(i)-ei)**2
      end do
      end

      subroutine pssrls_gr(y, basisval, k, t, coef, dcoef, lpen, phiphi,
     *     s)
      integer t, k, i
      double precision y(t), basisval(t,k), coef(k), lpen(k,k), sv(11)
      double precision phiphi(k,k,t),  pssr0, pssr1, dcoef(k), s, ct(k)
      double precision ctm(k), pssr1m
      sv(6) = 1.0d0
      do i=1,5
         sv(6+i) = sv(6+i-1)*1.3
         sv(6-i) = sv(6-i+1)/1.2
      end do
      do i=1,11
         ct = coef+sv(i)*dcoef
         ctm = coef-sv(i)*dcoef
         call getpssr(y, basisval, k, t, ct, lpen, phiphi, pssr1)
         call getpssr(y, basisval, k, t, ctm, lpen, phiphi, pssr1m)
         if (pssr1m < pssr1) then
            sv(i) = -sv(i)
            pssr1 = pssr1m
         end if
         
         if (i==1) then
            pssr0 = pssr1
            s = sv(i)
         else if (pssr1<pssr0) then
            pssr0 = pssr1
            s = sv(i)
         end if
      end do
      end

      subroutine pssrls_br(y, basisval, k, t, coef, dcoef, lpen, phiphi,
     *     maxita, tola, ax, cx, x)
      integer t, k, i, maxita, iter, info
      double precision y(t), basisval(t,k), coef(k), lpen(k,k), ct1(k)
      double precision phiphi(k,k,t),  pssr1, pssr2, dcoef(k), ct2(k)
      double precision ax, cx, bx, tola, tol1, tol2, tmpl(2)
      double precision gs, zeps, a, b, d, etemp, fu, fv, fw, fx
      double precision p,q,r,u,v,w,x,xm,e,u2,fu2
      gs = 0.3819660
      e = 0.0d0
      zeps = epsilon(gs)*1.0e-3
c     Find the third point between ax and cx
      tmpl(1) = (1-gs)*ax+gs*cx
      tmpl(2) = gs*ax+(1-gs)*cx
      ct1 = coef+tmpl(1)*dcoef
      ct2 = coef+tmpl(2)*dcoef      
      call getpssr(y, basisval, k, t, ct1, lpen, phiphi, pssr1)
      call getpssr(y, basisval, k, t, ct2, lpen, phiphi, pssr2)
      if (pssr1 < pssr2) then
         bx = tmpl(1); cx = tmpl(2); fx = pssr1
      else
         bx = tmpl(2); ax = tmpl(1); fx = pssr2
      end if
c     ############################
      x = bx; w = bx; v = x; fv = fx; fw = fx
      i = 1; a = ax; b = cx
      do
         iter = i
         if (i>maxita) then
            info = 1
            exit
         end if
         xm = 0.5*(b+a)
         tol1 = tola*abs(x)+zeps
         tol2 = 2.0*tol1
         if (abs(x-xm) <= (tol2-0.5*(b-a))) then
            exit
         end if
         if (abs(e) > tol1) then
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q-(x-w)*r
            q = 2.0*(q-r)
            if (q>0.0) then
               p = -p
            end if
            q = abs(q)
            etemp = e
            e = d
            if (abs(p) >= abs(0.5*q*etemp) .or. p <= q*(a-x) .or.
     *           p >= q*(b-x)) then
               if (x >= xm) then
                  e = a-x
               else
                  e = b-x
               end if               
               d = gs*e
            else
               d = p/q
               u = x+d
               if (u-a < tol2 .or. b-u < tol2) then
                  d = sign(tol1,xm-x)
               end if
            end if   
         else
            if (x >= xm) then
               e = a-x
            else
               e = b-x
            end if               
            d = gs*e
         end if
         if (abs(d) >= tol1) then
            u=x+d
         else
            u = x+sign(tol1,d)
         end if
         u2 = u
         ct1 = coef+u2*dcoef         
         call getpssr(y, basisval, k, t, ct1, lpen, phiphi, fu2)
         fu = fu2
         if (fu <= fx) then
            if (u >= x) then
               a=x
            else
               b=x
            end if
            v = w; w = x; x = u; fv = fw; fw = fx; fx = fu
         else
            if (u<x) then
               a=u
            else
               b=u
            end if
            if (fu <= fw .or. w == x) then
               v = w; w = u; fv = fw; fw = fu
            else if (fu <= fv .or. v == x) then
               v = u; fv = fu
            end if
         end if
         i = i+1
      end do
      end
   
c     Main function to estimate the coefficients of the non-negative fda
c     basisval and pen are generated by the fda R package
c     info is nx1 and returns convergence codes:
c     0: normal convergence, 1: failure to solve Hx=J, 2: maxit reached
c     *******************************************************************     
      subroutine nlcoefest(y, basisval, n, k, t, lambda, pen,
     *     tol, maxit, info, coef, typels)
      integer n, k, t, maxit, p(k), lwork, info(n), info2, i, j
      integer maxita
      integer lworkmax
      parameter (lworkmax = 5000)
      double precision basisval(t,k), lambda, coef(k,n), phiphi(k,k,t)
      double precision y(t,n), pen(k,k), jac(k), hes(k,k), lpen(k,k)
      double precision c0(k,n), tol, crit, work(lworkmax), n1, n2
      double precision pssr(n), s, tola, ax, cx
      character typels
      
      lpen = 2*lambda*pen
      tola = 1.0e-4
      maxita = 50
      
      call prepopt(basisval, k, t, phiphi)
      call initialc(y, basisval, n, k, t, c0)
      do j=1,n
      i = 1      
      do
         ax = -4.0d0
         cx = 4.0d0
         call hesjac(y(:,j), basisval, k, t, c0(:,j), lpen, phiphi,
     *        hes, jac, pssr(j))
         lwork = -1
         call dsysv('U', k, 1, hes, k, p, jac, k, work, lwork, info2)
         lwork = min(lworkmax, int(work(1)))
         call dsysv('U', k, 1, hes, k, p, jac, k, work, lwork, info2)
         if (info2 > 0) then
            info(j) = 1
            exit
         end if
         if (typels == 'g') then 
            call pssrls_gr(y(:,j), basisval, k, t, c0(:,j), -jac, lpen,
     *           phiphi, s)
         else            
            call pssrls_br(y(:,j), basisval, k, t, c0(:,j), -jac, lpen,
     *           phiphi, maxita, tola, ax, cx, s)
         end if
         coef(:,j) = c0(:,j) - s*jac
         n1 = sqrt(sum(jac**2))
         n2 = sqrt(sum(c0(:,j)**2))
         crit = n1/(1+n2)
         if (crit < tol) then
            info(j) = 0
            exit
         end if
         if (i > maxit) then
            info(j) = 2
            exit
         end if
         c0(:,j) = coef(:,j)
         i = i+1
      end do
      end do
      end
      

c     Codes to get the coefficients and yhat's for non-linear fda
c     and observation i removed
c     ***********************************************************
      subroutine nlyhati(y, basisval, n, k, t, lambda, pen, i, 
     *     tol, maxit, info, yhat, typels)
      integer n, k, t, s(t-1), i, j, info(n), maxit
      double precision basisval(t,k), lambda, coef(k,n)
      double precision y(t,n), pen(k,k), yhat(n), tol
      character typels
      do j=1,(t-1)
         if (j < i) then
            s(j) = j
         else
            s(j) = j+1
         end if
      end do
      call nlcoefest(y(s,:), basisval(s,:), n, k, t-1, lambda, pen,
     *     tol, maxit, info, coef, typels)
      call dgemv('T', k, n, 1.0d0, coef, k, basisval(i,:), 1, 0.0d0,
     *     yhat, 1)
      yhat = exp(yhat)
      end 
      
c     Codes to get the leave-one-out cross-validation for non-linear fda
c     lambda and k can be vectors.
c     The result is a nlam x nk matrix of CV's.
c     ***********************************************************      
      subroutine nlcrval(y, basisval, n, k, t, nlam, nk, maxk, lambda, 
     *     pen, tol, maxit, info, cv, typels)
      integer n, t, i, j, nlam, maxit, nk, l
      integer k(nk), info(t,n,nlam,nk)
      double precision y(t,n), basisval(t,maxk,nk)
      double precision lambda(nlam), pen(maxk,maxk,nk)
      double precision cv(nlam,nk), yhat(n), tol
      character typels
      cv = 0.0d0
      do l=1,nk
         do j=1,nlam
            do i=1,t
               call nlyhati(y, basisval(:,1:k(l),l), n, k(l), t,  
     *              lambda(j), pen(1:k(l), 1:k(l), l), i, tol, maxit,
     *              info(i,:,j,l), yhat, typels)
               cv(j,l) = cv(j,l) + sum((yhat-y(i,:))**2)
            end do
         end do
      end do
      end

c     Minimization of the cross-validation with respect to
c     lambda for a given k using the Golden-Section method
c     The nonlinear case. 
c     ********************************************************* 
      subroutine nllamgs(y, basisval, n, k, t, l1, l2, 
     *     pen, tol, maxit, maxitalgo, tolalgo, info, iter, lam, cvf,
     *     typels)
      integer n, t, k, maxit, maxitalgo, i, k2(1)
      integer info2(t,n,2,1), info, iter
      double precision y(t,n), basisval(t,k), l1, l2
      double precision lam, pen(k,k), tmpl(2), cvf
      double precision cv(2), tol, gspar, tolalgo
      character typels

      gspar = 1-0.3819660
      tmpl(1) = gspar*l1+(1-gspar)*l2
      tmpl(2) = (1-gspar)*l1+gspar*l2
      info = 0
      k2(1) = k
      i = 1
      call nlcrval(y, basisval, n, k2, t, 2, 1, k, tmpl, 
     *     pen, tol, maxit, info2(:,:,:,1), cv, typels)
      do
         iter = i
         if (cv(1) < cv(2)) then
            l2  = tmpl(2)
            tmpl(2) = tmpl(1)
            cv(2) = cv(1)
            lam = tmpl(2)
            cvf = cv(2)
            tmpl(1) = gspar*l1+(1-gspar)*l2
            call nlcrval(y, basisval, n, k2, t, 1, 1, k, tmpl(1), 
     *           pen, tol, maxit, info2(:,:,1,1), cv(1), typels)
         else
            l1 = tmpl(1)
            tmpl(1) = tmpl(2)
            lam = tmpl(1)
            cv(1) = cv(2)
            cvf = cv(1)
            tmpl(2) = (1-gspar)*l1+gspar*l2
            call nlcrval(y, basisval, n, k2, t, 1, 1, k, tmpl(2), 
     *           pen, tol, maxit, info2(:,:,2,1), cv(2), typels)            
         end if
         if (abs(l2-l1) < lam*tolalgo) then
            exit
         end if
         if (i>maxitalgo) then
            info = 1
            exit
         end if
         i = i+1
      end do
      end

c     Minimization of the cross-validation with respect to
c     lambda for a given k using the Brent method
c     The nonlinear case. 
c     *********************************************************       
      subroutine nllambrent(y, basisval, n, k, t, ax, cx, 
     *     pen, tol, maxit, maxita, tola, info, iter, x, fx,
     *     typels)
      integer n, t, k, maxit, maxita, i, k2(1)
      integer info2(t,n,2,1), info, iter
      double precision y(t,n), basisval(t,k), l1, l2
      double precision pen(k,k), tmpl(2), ax, cx, bx
      double precision cv(2), tol, tola, tol1, tol2
      double precision gs, zeps, a, b, d, etemp, fu, fv, fw, fx
      double precision p,q,r,u,v,w,x,xm,e,u2(1), fu2(1)
      character typels
      gs = 0.3819660
      e = 0.0d0
      zeps = epsilon(gs)*1.0e-3
      k2(1) = k
c     Find the third point between ax and cx
      tmpl(1) = (1-gs)*ax+gs*cx
      tmpl(2) = gs*ax+(1-gs)*cx
      call nlcrval(y, basisval, n, k2, t, 2, 1, k, tmpl, 
     *     pen, tol, maxit, info2(:,:,:,1), cv, typels)
      if (cv(1) < cv(2)) then
         bx = tmpl(1); cx = tmpl(2); fx = cv(1)
      else
         bx = tmpl(2); ax = tmpl(1); fx = cv(2)
      end if
c     ############################
      x = bx; w = bx; v = x; fv = fx; fw = fx
      i = 1; a = ax; b = cx; info = 0
      do
         iter = i
         if (i>maxita) then
            info = 1
            exit
         end if
         xm = 0.5*(b+a)
         tol1 = tola*abs(x)+zeps
         tol2 = 2.0*tol1
         if (abs(x-xm) <= (tol2-0.5*(b-a))) then
            exit
         end if
         if (abs(e) > tol1) then
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q-(x-w)*r
            q = 2.0*(q-r)
            if (q>0.0) then
               p = -p
            end if
            q = abs(q)
            etemp = e
            e = d
            if (abs(p) >= abs(0.5*q*etemp) .or. p <= q*(a-x) .or.
     *           p >= q*(b-x)) then
               if (x >= xm) then
                  e = a-x
               else
                  e = b-x
               end if               
               d = gs*e
            else
               d = p/q
               u = x+d
               if (u-a < tol2 .or. b-u < tol2) then
                  d = sign(tol1,xm-x)
               end if
            end if   
         else
            if (x >= xm) then
               e = a-x
            else
               e = b-x
            end if               
            d = gs*e
         end if
         if (abs(d) >= tol1) then
            u=x+d
         else
            u = x+sign(tol1,d)
         end if
         u2(1) = u
         call nlcrval(y, basisval, n, k2, t, 1, 1, k, u2, 
     *        pen, tol, maxit, info2(:,:,1,1), fu2, typels)
         fu = fu2(1)
         if (fu <= fx) then
            if (u >= x) then
               a=x
            else
               b=x
            end if
            v = w; w = x; x = u; fv = fw; fw = fx; fx = fu
         else
            if (u<x) then
               a=u
            else
               b=u
            end if
            if (fu <= fw .or. w == x) then
               v = w; w = u; fv = fw; fw = fu
            else if (fu <= fv .or. v == x) then
               v = u; fv = fu
            end if
         end if
         i = i+1
      end do
      end

