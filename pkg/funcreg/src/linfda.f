c     ***********************************************************
c     ***********************************************************
c                             LINEAR FDA
c     ***********************************************************
c     ***********************************************************      
      
c     Codes to estimate the coefficients for linear fda
c     ***********************************************************
      subroutine coefest(y, basisval, n, k, t, lambda, pen, coef)
      integer n, k, t, p(k), infotr, infosv
      double precision basisval(t,k), lambda, coef(k,n)
      double precision y(t,n), pen(k,k), tbasis(k,t), pen2(k,k)
      
      pen2 = pen
      call dgemm('T', 'N', k, k, t, 1.0d0, basisval, t, basisval, t,
     *     lambda*t, pen2, k)
      call dgetrf(k, k, pen2, k, p, infotr)
      tbasis = transpose(basisval)
      call dgetrs('N', k, t, pen2, k, p, tbasis, k, infosv)
      
      call dgemm('N', 'N', k, n, t, 1.0d0, tbasis, k, y, t,
     *     0.0d0, coef, k)

      end

c     Codes to get the coefficients and yhat's for linear fda
c     and observation i removed
c     ***********************************************************
      subroutine yhati(y, basisval, n, k, t, lambda, pen, i, coef, yhat)
      integer n, k, t, s(t-1), i, j
      double precision basisval(t,k), lambda, coef(k,n)
      double precision y(t,n), pen(k,k), yhat(n)

      do j=1,(t-1)
         if (j < i) then
            s(j) = j
         else
            s(j) = j+1
         end if
      end do      
      call coefest(y(s,:), basisval(s,:), n, k, t-1, lambda, pen,
     *     coef)
      call dgemv('T', k, n, 1.0d0, coef, k, basisval(i,:), 1, 0.0d0,
     *     yhat, 1)
   
      end 


c     Codes to get the leave-one-out cross-validation for linear fda
c     lambda can be a vector, which will produce a vector of CV.
c     ***********************************************************      
      subroutine crval(y, basisval, n, k, t, nlam, lambda, pen, cv)
      integer n, k, t, i, j, nlam
      double precision y(t,n), basisval(t,k), lambda(nlam), pen(k,k)
      double precision cv(nlam), coef(k,n), yhat(n)
      cv = 0.0d0
      do j=1,nlam
         do i=1,t
            call yhati(y, basisval, n, k, t, lambda(j), pen, i, coef,
     *           yhat)
            cv(j) = cv(j) + sum((yhat-y(i,:))**2)
         end do
      end do
      
      end

c     Minimization of the cross-validation with respect to
c     lambda for a given k using the Golden-Section method
c     The linear case. 
c     ********************************************************* 
      subroutine lamgs(y, basisval, n, k, t, l1, l2, 
     *     pen, maxitalgo, tolalgo, info, iter, lam, cvf)
      integer n, t, k, maxitalgo, i, info, iter
      double precision y(t,n), basisval(t,k), l1, l2
      double precision lam, pen(k,k), tmpl(2), cvf
      double precision cv(2), gspar, tolalgo

      gspar = 1-0.3819660
      tmpl(1) = gspar*l1+(1-gspar)*l2
      tmpl(2) = (1-gspar)*l1+gspar*l2
      info = 0
      i = 1
      call crval(y, basisval, n, k, t, 2, tmpl, pen, cv)
      do
         iter = i
         if (cv(1) < cv(2)) then
            l2  = tmpl(2)
            tmpl(2) = tmpl(1)
            cv(2) = cv(1)
            lam = tmpl(2)
            cvf = cv(2)
            tmpl(1) = gspar*l1+(1-gspar)*l2
            call crval(y, basisval, n, k, t, 1, tmpl(1), pen, cv(1))
         else
            l1 = tmpl(1)
            tmpl(1) = tmpl(2)
            lam = tmpl(1)
            cv(1) = cv(2)
            cvf = cv(1)
            tmpl(2) = (1-gspar)*l1+gspar*l2
            call crval(y, basisval, n, k, t, 1, tmpl(2), pen, cv(2))
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
c     The linear case. 
c     *********************************************************       
      subroutine lambrent(y, basisval, n, k, t, ax, cx, 
     *     pen, maxita, tola, info, iter, x, fx)
      integer n, t, k, maxita, i, info, iter
      double precision y(t,n), basisval(t,k), l1, l2
      double precision pen(k,k), tmpl(2), ax, cx, bx
      double precision cv(2), tola, tol1, tol2
      double precision gs, zeps, a, b, d, etemp, fu, fv, fw, fx
      double precision p,q,r,u,v,w,x,xm,e,u2(1), fu2(1)
      gs = 0.3819660
      e = 0.0d0
      zeps = epsilon(gs)*1.0e-3
c     Find the third point between ax and cx
      tmpl(1) = (1-gs)*ax+gs*cx
      tmpl(2) = gs*ax+(1-gs)*cx
      call crval(y, basisval, n, k, t, 2, tmpl, pen, cv)
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
         call crval(y, basisval, n, k, t, 1, u2, pen, fu2)
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

      

