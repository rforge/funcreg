      subroutine kron(a,m,n,b,l,s,ab)
      integer m,n,l,s,i,j,i0,i1,j0,j1,t
      double precision a(m,n), b(l,s),ab0(m*n,l*s)
      double precision ab(m*l, n*s)

      ab0 = 0.0d0
      t = 1
      call dger(m*n, l*s, 1.0d0, a, 1, b, 1, ab0, m*n)
      do i=1,n
         i0 = s*(i-1)+1
         i1 = i0+s-1
         do j=1,m
            j0 = l*(j-1)+1
            j1 = j0+l-1
            ab(j0:j1, i0:i1) = reshape(ab0(t,:),(/l, s/))
            t = t+1
         end do
      end do      
      end


c     Functions to estimate a functional regression
c     **********************************************


c     Getting (W(t)'W(t)) and (W(t)'y(t))
c     ipbbt <phi_1,phi_1'>,...,<phi_nx,phi_1'>,<phi_2,phi_2'>,
c     ...<phi_nx,phi_2'>,<phi_3,.phi_3'>,.....
c     kb(:,1) are the number of t-basis and kb(:,2) the number of s-basis
      
      subroutine prepfr(g, h, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, cmat, dmat, lam0, lam)
      integer n, nx, k0, kb(nx,2), mkbt, mkbs, ncoef, j, i, start
      integer ncoefi, sc1, sc2, t, sc12, sc22, sc10
      double precision fmat(n, k0), ipbbt(mkbt,mkbt,nx*(nx+1)/2)
      double precision ipabt(k0, mkbt, nx), ipaa(k0,k0)
      double precision rs(mkbs*mkbt, mkbs*mkbt, nx), r0(k0,k0)
      double precision rt(mkbs*mkbt, mkbs*mkbt, nx), h(n,mkbs,nx)
      double precision g(n, mkbt, nx), cmat(ncoef,ncoef), dmat(ncoef)
      double precision sfmat(k0), lam0, lam(nx, 2)

      start = 1
      sc1 = 1
      sc10 = 1
      sc2 = 1
      if (k0>0) then
         sfmat = sum(fmat, dim=1)
         dmat(1:k0) = sfmat
         start = start+k0
         cmat(1:k0,1:k0) = n*ipaa + lam0*r0
         do i=1,nx
            ncoefi = product(kb(i,:))
            call vdprod(g(:,1:kb(i,1),i), n, kb(i,1), h(:,1:kb(i,2),i),
     *           kb(i,2), dmat(start:(start+ncoefi-1)))
            call kron(sum(h(:,1:kb(i,2),i), dim=1), kb(i,2), 1,
     *           transpose(ipabt(:,1:kb(i,1),i)), kb(i,1), k0,
     *           cmat(start:(start+ncoefi-1), 1:k0))
            start = start+ncoefi
         end do
         sc1 = k0+1
         sc2 = k0+1
         sc10 = k0+1
      end if

      t = 1
      do j=1,nx
         sc22 = sc2+kb(j,1)*kb(j,2)-1         
         do i=j,nx
            sc12 = sc1+kb(i,1)*kb(i,2)-1
            call cpkron(h(:,1:kb(i,2),i), n, kb(i,2), h(:,1:kb(j,2),j),
     *           kb(j,2), ipbbt(1:kb(i,1), 1:kb(j,1),t),kb(i,1),kb(j,1),
     *            cmat(sc1:sc12, sc2:sc22))
            if (i==j) then
               cmat(sc1:sc12, sc2:sc22) = cmat(sc1:sc12, sc2:sc22) +
     *         lam(i,1)*rt(1:(kb(i,1)*kb(i,2)),1:(kb(i,1)*kb(i,2)),i) +
     *         lam(i,2)*rs(1:(kb(i,1)*kb(i,2)),1:(kb(i,1)*kb(i,2)),i) 
            end if
            t = t+1
            sc1 = sc1+kb(i,1)*kb(i,2)
         end do
         sc2 = sc2+kb(j,1)*kb(j,2)
         sc10 = sc10 + kb(j,1)*kb(j,2)
         sc1 = sc10
      end do
      end


      subroutine vdprod(a,m,n,b,l,c)
      integer m, n, l
      double precision a(m,n), b(m,l), c(n*l)

      c = 0.0d0
      if (m>1) then 
         call dgemm('T','N', n, l, m, 1.0d0, a, m, b, m, 0.0d0,
     *        c, n)
      else
         call dger(n, l, 1.0d0, a(1,:), 1, b(1,:), 1, c, n)
      end if
      
      end

      subroutine cpkron(a,m,n,b,l,c,r,s,ans)
      integer m,n,l,r,s
      double precision a(m,n), b(m,l), c(r,s)
      double precision ans(n*r, l*s), h(n,l)

      if (m>1) then
         call dgemm('T','N',n,l,m,1.0d0, a, m, b, m, 0.0d0, h, n)
      else
         h = 0.0d0
         call dger(n, l, 1.0d0, a(1,:), 1, b(1,:), 1, h, n)
      end if      
      call kron(h, n, l, c, r, s, ans)
      end

c     This is the one called from R to simply estimate the model
C     And provide the information to compute the variance
      subroutine funcregwv(g, h, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, lam0, lam,
     *     coef, info, cmat, dmat, cmat0)
      integer n, nx, k0, kb(nx,2), mkbt, mkbs, ncoef, ip(ncoef)
      integer info, lworkmax, lwork, info2
      parameter (lworkmax = 5000)
      double precision fmat(n, k0), ipbbt(mkbt,mkbt,nx*(nx+1)/2)
      double precision ipabt(k0, mkbt, nx), ipaa(k0,k0), ipyy(n)
      double precision rs(mkbs*mkbt, mkbs*mkbt, nx), r0(k0,k0)
      double precision rt(mkbs*mkbt, mkbs*mkbt, nx), h(n,mkbs,nx)
      double precision g(n, mkbt, nx), cmat(ncoef,ncoef), coef(ncoef)
      double precision lam0, lam(nx, 2), work(lworkmax)
      double precision lam02, lam2(nx,2), cmat0(ncoef,ncoef),dmat(ncoef)

      call prepfr(g, h, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, cmat, coef, lam0, lam)
      info = 0
      lwork = -1
      call dsysv('L', ncoef, 1, cmat, ncoef, ip, coef, ncoef, work,
     *     lwork, info2)
      lwork = min(lworkmax, int(work(1)))
      call dsysv('L', ncoef, 1, cmat, ncoef, ip, coef, ncoef, work,
     *     lwork, info2)
      if (info2 > 0) then
         info = 1
      end if
      lam02 = 0.0d0
      lam2 = 0.0d0
      call prepfr(g, h, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, cmat0, dmat, lam02, lam2)
      
      end


      subroutine funcreg(g, h, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, lam0, lam,
     *     coef, info)
      integer n, nx, k0, kb(nx,2), mkbt, mkbs, ncoef, ip(ncoef)
      integer info, lworkmax, lwork, info2
      parameter (lworkmax = 5000)
      double precision fmat(n, k0), ipbbt(mkbt,mkbt,nx*(nx+1)/2)
      double precision ipabt(k0, mkbt, nx), ipaa(k0,k0), ipyy(n)
      double precision rs(mkbs*mkbt, mkbs*mkbt, nx), r0(k0,k0)
      double precision rt(mkbs*mkbt, mkbs*mkbt, nx), h(n,mkbs,nx)
      double precision g(n, mkbt, nx), cmat(ncoef,ncoef), coef(ncoef)
      double precision lam0, lam(nx, 2), work(lworkmax)

      call prepfr(g, h, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, cmat, coef, lam0, lam)
      info = 0
      lwork = -1
      call dsysv('L', ncoef, 1, cmat, ncoef, ip, coef, ncoef, work,
     *     lwork, info2)
      lwork = min(lworkmax, int(work(1)))
      call dsysv('L', ncoef, 1, cmat, ncoef, ip, coef, ncoef, work,
     *     lwork, info2)
      if (info2 > 0) then
         info = 1
      end if      
      end

   
      subroutine funcregi(g, h, i, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, lam0, lam,
     *     ipyy, coef, cv, info)
      integer n, nx, k0, kb(nx,2), mkbt, mkbs, ncoef, i, info
      integer j, s(n-1)
      double precision fmat(n, k0), ipbbt(mkbt,mkbt,nx*(nx+1)/2)
      double precision ipabt(k0, mkbt, nx), ipaa(k0,k0), ipyy
      double precision rs(mkbs*mkbt, mkbs*mkbt, nx), r0(k0,k0)
      double precision rt(mkbs*mkbt, mkbs*mkbt, nx), h(n,mkbs,nx)
      double precision g(n, mkbt, nx), coef(ncoef)
      double precision lam0, lam(nx, 2), cv

      do j=1,(n-1)
         if (j < i) then
            s(j) = j
         else
            s(j) = j+1
         end if
      end do      
      
      call funcreg(g(s,:,:), h(s,:,:), r0, rt, rs, ipaa, ipabt, ipbbt,
     *     fmat(s,:),
     *     n-1, nx, k0, kb, mkbt, mkbs, ncoef, lam0, lam, coef, info)
      
      call inprodresi(g(i,:,:), h(i,:,:), ipaa, ipabt, ipbbt, fmat(i,:),
     *     nx, k0, kb, mkbt, mkbs, ncoef, coef, ipyy, cv)      
      end

      
      subroutine prepfri(g, h, ipaa, ipabt, ipbbt, fmat,
     *     nx, k0, kb, mkbt, mkbs, ncoef, cmat, dmat)
      integer nx, k0, kb(nx,2), mkbt, mkbs, ncoef, j, i, start
      integer ncoefi, sc1, sc2, t, sc12, sc22, sc10
      double precision fmat(k0), ipbbt(mkbt,mkbt,nx*(nx+1)/2)
      double precision ipabt(k0, mkbt, nx), ipaa(k0,k0), h(mkbs,nx)
      double precision g(mkbt, nx), cmat(ncoef,ncoef), dmat(ncoef)
      
      start = 1
      sc1 = 1
      sc10 = 1
      sc2 = 1
      if (k0>0) then
         dmat(1:k0) = fmat
         start = start+k0
         cmat(1:k0,1:k0) = ipaa
         do i=1,nx
            ncoefi = product(kb(i,:))
            call vdprod(g(1:kb(i,1),i), 1, kb(i,1), h(1:kb(i,2),i),
     *           kb(i,2), dmat(start:(start+ncoefi-1)))
            call kron(h(1:kb(i,2),i), kb(i,2), 1,
     *           transpose(ipabt(:,1:kb(i,1),i)), kb(i,1), k0,
     *           cmat(start:(start+ncoefi-1), 1:k0))
            start = start+ncoefi
         end do
         sc1 = k0+1
         sc2 = k0+1
         sc10 = k0+1
      end if

      t = 1
      do j=1,nx
         sc22 = sc2+kb(j,1)*kb(j,2)-1         
         do i=j,nx
            sc12 = sc1+kb(i,1)*kb(i,2)-1
            call cpkron(h(1:kb(i,2),i), 1, kb(i,2), h(1:kb(j,2),j),
     *           kb(j,2), ipbbt(1:kb(i,1), 1:kb(j,1),t),kb(i,1),kb(j,1),
     *           cmat(sc1:sc12, sc2:sc22))
            t = t+1
            sc1 = sc1+kb(i,1)*kb(i,2)
         end do
         sc2 = sc2+kb(j,1)*kb(j,2)
         sc10 = sc10 + kb(j,1)*kb(j,2)
         sc1 = sc10
      end do

      end
      
      subroutine sysquadra(a,n,x,q)
      integer i, j, n
      double precision a(n,n), x(n), q, s

      q = 0.0d0
      do i=1,n
         do j=i,n
            s = a(j,i)*x(i)*x(j)
            if (i==j) then
               q = q + s
            else
               q = q+2*s
            end if
         end do
      end do

      end
      
      subroutine inprodresi(g, h, ipaa, ipabt, ipbbt, fmat,
     *     nx, k0, kb, mkbt, mkbs, ncoef, coef, ipyy, ipres)
      integer nx, k0, kb(nx,2), mkbt, mkbs, ncoef, i
      double precision fmat(k0), ipbbt(mkbt,mkbt,nx*(nx+1)/2)
      double precision ipabt(k0, mkbt, nx), ipaa(k0,k0), h(mkbs,nx)
      double precision g(mkbt, nx), cmat(ncoef,ncoef), dmat(ncoef)
      double precision ipyy, ipres, tmp, coef(ncoef)

      call prepfri(g, h, ipaa, ipabt, ipbbt, fmat,
     *     nx, k0, kb, mkbt, mkbs, ncoef, cmat, dmat)

      call sysquadra(cmat, ncoef, coef, tmp)
      ipres = ipyy-2*sum(coef*dmat)+tmp
      end
      
      

      subroutine funcregcv(g, h, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *     n, nx, k0, kb, mkbt, mkbs, ncoef, lam0, lam,
     *     ipyy, cv, info)
      integer n, nx, k0, kb(nx,2), mkbt, mkbs, ncoef, info(n), i
      double precision fmat(n, k0), ipbbt(mkbt,mkbt,nx*(nx+1)/2)
      double precision ipabt(k0, mkbt, nx), ipaa(k0,k0), ipyy(n)
      double precision rs(mkbs*mkbt, mkbs*mkbt, nx), r0(k0,k0)
      double precision rt(mkbs*mkbt, mkbs*mkbt, nx), h(n,mkbs,nx)
      double precision g(n, mkbt, nx), coef(ncoef), cvv(n), cv
      double precision lam0, lam(nx, 2) 

      do i=1,n
         call funcregi(g, h, i, r0, rt, rs, ipaa, ipabt, ipbbt, fmat,
     *        n, nx, k0, kb, mkbt, mkbs, ncoef, lam0, lam,
     *        ipyy(i), coef, cvv(i), info(i))
      end do

      cv = sum(cvv)
      end

