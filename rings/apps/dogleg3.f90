      subroutine dogleg(n, func, grad, x, f, doptions, ioptions, &
     &     nthts, nk, nnu, pow, fnu, nrdtot, verbose)

!*****************************************************************************
!     
!     DOGLEG: Powell's hybrid method for unconstrained optimization.
!     
!*****************************************************************************
!     
!     PARAMETERS
!     
!     n      :  integer
!       Problem dimension
!     
!     FUNC   :  double precision function func(N,X,IPARM,RPARM) 
!       User-supplied routine calculating f(X), the function
!       to be minimized.
!     
!     GRAD   :  subroutine GRAD(N,X,G,IPARM,RPARM)
!       User-supplied subroutine that calculates the gradient of f
!       F'(X) and returns it in G. GRAD is not called if igrad=0.
!     
!     X      :  double precision array of dimension n.
!       On input  : Starting guess.
!       On output : The minimizer found.
!     
!     FX     :  On output : Function value f(X) at the minimizer.
!     
!     doptions(1) = eps1   :  double precision
!       Termination tolerance on X.
!     
!     doptions(2) = eps2   :  double precision
!       Termination tolerance on F.
!     
!     doptions(3) = feps   :  double precision
!       (Approximately) The accuracy of the function values 
!       calculated in FUNC.
!     
!     doptions(4) = delta0 :  double precision
!       Initial radius of the trust region.
!     
!     ioptions(1) = nfmax  :  integer
!       On input  : Maximum number of function evaluations.
!       On output : The number of function evaluations used.
!     
!     ioptions(2) = igrad  :  integer
!       If igrad is 0 then F' is calculated numerically.
!       If igrad is 1 then the user supplied subroutine GRAD is called
!       to calculate F'.
!       If igrad is 2 then the user supplied subroutine GRAD is called
!       to calculate F' and initially the results from GRAD are checked
!       by comparing them to the F' obtained using finite differences.
!     
!     ioptions(3) = iextrap  :  integer
!       If iextrap is 1 then safeguarded quadratic/cubic interpolation 
!       and extrapolation is performed.
!     
!     ioptions(4) = idiag  :  integer
!       If idiag is 0 then I is used as the initial approximation to the 
!       Hessian.
!       If idiag is 1 then the diagonal of the true Hessian is computed
!       initially using a finite difference approximation.
!     
!     ioptions(5) = iprint  :  integer
!       If iprint=1 then information is printed after each iteration.
!     
!     (C) Rasmus Munk Larsen, University of Aarhus, March 1998.

      implicit none
      integer, parameter :: double=8
      integer n,ioptions(5)
      external func,grad
!      double precision func,x(*),doptions(4)
      real(kind=double) :: func,x(*),doptions(4)
!     
!     Local variables
!     
      integer MAXN,lwork,ichkgrd
      parameter(MAXN=1000,lwork=(MAXN*(MAXN+1))/2)
      integer nf,ng,istop,i,j,newton,ipiv(MAXN),info,id1,id2,k,nfail 
      integer ierr,nfmax,igrad,iprint,iextrap,idiag,verbose,nk,nthts,nnu
      integer nrdtot
      double precision f,f1,f2,f3,h(MAXN),hn(MAXN),x1(MAXN), &
          vtemp(MAXN),x2(MAXN),x3(MAXN), &
          g(MAXN),g1(MAXN),HH(lwork), B(lwork)
      double precision alpha,alpha2,aa,bb,ab,d,theta,dL,dF,rho, &
          normh,normhn,normx,normg,normg0,normg2,fdeps,sfmin, &
          normy,fdtol,delta,gh,seps,eps,zero,one,fdstep,cdstep, &
          eps1,eps2,delta0,feps
      real(kind=double),  dimension(nk,nrdtot) :: fnu
      real(kind=double),  dimension(nthts,nk,nnu) :: pow
      character*20 how
      logical extrap
      parameter(one=1d0,zero=0d0)
!     
!     External subroutine and functions
!     
      external  dnrm2,ddot,extrap3p,dlamch,cubic3
      double precision dnrm2,ddot,extrap3p,dlamch,cubic3

!************************Here begins dogleg ************************


      if (n.gt.MAXN) then
         write(*,*) 'DOGLEG only dimensioned for n <= ',MAXN
         stop 'Problem too Large'
      endif
!     
!     Initialize.
!     
      if (iprint.eq.1) then
         write(6,*) 'ENTERING DOGLEG'
         write (6,*) 'n,delta0,eps1,eps2,nfmax,igrad,iprint = ', &
     &        n,delta0,eps1,eps2,nfmax,igrad,iprint 
         call flush(6)
      endif
!     
!     Set tolerances
!     
      eps1 = doptions(1)
      eps2 = doptions(2)
      feps = doptions(3)
      delta0 = doptions(4)
!
!     Set algorithmic switches
!
      nfmax = ioptions(1)
      igrad = ioptions(2)
      if (igrad .eq. 2) then
         igrad = 1
         ichkgrd = 1
      else
         ichkgrd = 0
      endif
      iextrap = ioptions(3)
      idiag = ioptions(4)
      iprint = ioptions(5)

      eps = dlamch('e')
      sfmin = dlamch('s')
      seps = sqrt(eps)
      fdstep = sqrt(feps)
      fdtol = sqrt(fdstep)
      cdstep = feps**(1d0/3d0)
      eps1 = max(eps1,eps)
      eps2 = max(eps2,eps)

!     
!     Calculate function value and gradient in starting point
      
      f = func(n, x, nthts, nk, nnu, pow, fnu, nrdtot)
      nf = 1

      if (ichkgrd.eq.1) then
         write (*,*) 'CHECKING USER SUPPLIED GRADIENTS'
         fdeps = 0.5d0
         call grad(n,x,vtemp, nthts, nk, nnu, pow, fnu, nrdtot)
         write (6,*) '     eps         err1          err2'
         do i=1,20
            call fdgrad(func,fdeps,n,x,f,g,nthts,nk,nnu,pow, fnu, nrdtot)
            write (*,*) 'grad   = ',(vtemp(j),j=1,n)
            write (*,*) 'fdgrad = ',(g(j),j=1,n)
            write (*,*) 'err FD = ',((g(j)-vtemp(j))/vtemp(j),j=1,n)
           f1 = 0d0
            do j=1,n
               if (abs(g(j)-vtemp(j)).gt. abs(f1)) then
                  f1 = g(j)-vtemp(j)
                  id1 = j
               endif
            enddo
            call cdgrad(func,fdeps,n,x,f,g,nthts,nk,nnu,pow, fnu, nrdtot)
!            write (*,*) (g(j)-vtemp(j),j=1,n)
            f2 = 0d0
            do j=1,n
               if (abs(g(j)-vtemp(j)) .gt. abs(f2)) then
                  f2 = g(j)-vtemp(j)
                  id2 = j
               endif
            enddo
            write (6,11) fdeps,f1,f2,id1,id2
 11         format(1p3e13.4,2i4)
            fdeps = fdeps/4.0d0
         enddo
      endif
         
      if (igrad.eq.1) then
         call grad(n,x,g, nthts, nk, nnu, pow, fnu, nrdtot)
         ng = 1
      else
         call cdgrad1(func,cdstep,n,x,f,g,vtemp,nthts,nk,nnu,pow,fnu,nrdtot)
         ng = 2
         nf = nf + n
      endif
     
!     
!     Set up the Hessian (B) and its inverse HH stored in packed format. 
!     
      do i=1,(n*(n+1))/2
         HH(i) = 0D0      
         B(i) = 0D0
      enddo
      
      if (idiag.eq.1) then
         j = 1
         do i=1,n
            B(j) = vtemp(i)
            HH(j) = 1d0/vtemp(i)
            j = j+(i+1)
         enddo      
         
      else
!     HH0 =  B0 = I
         j = 1
         do i=1,n
            HH(j) = 1D0
            B(j) = 1d0
            j = j+(i+1)
         enddo      
      endif
      
      
      normg = dnrm2(n,g,1)
      normg2 = normg*normg
      normg0 = normg
      if (normg.eq.0d0) then
         write (*,*) 'DOGLEG WARNING: Zero gradient at starting point.'
         goto 90
      endif      
      normx = dnrm2(n,x,1)

!     Set initial width of trust region. 
      if (delta0 .le. 0d0)  delta0 = 0.5*normx
      if (delta0 .le. 0d0)  delta0 = 0.5
      delta = min(max(delta0,eps1*(normx+1)),100*(1+normx))


      istop = 0
      if (iprint.eq.1) then
         write (6,*) 'normx,delta = ',normx,delta
         how = ''      
         write (*,*) '  #f   #g         f       ' &
     &        // '   norm(g)      delta        df/dfold     step type'
         write (*,100) nf,ng,f,normg,delta,0.0,'       '
      endif

!     
!     Begin iteration
!     
      k = 0
      nfail = 0
      do while ((istop.eq.0).and.(nf.lt.nfmax))

!     
!     Compute Newton step h = -inv(B)*g
!     
         call dspmv('u',n,-1D0,HH,g,1,0d0,hn,1)
         normh = dnrm2(n,hn,1)
         normhn = normh

!     
!     Calculate length of steepest descent step alpha=g'*g/(g'*B*g)
!     
         call dspmv('u',n,1D0,B,g,1,0d0,vtemp,1)
         d = ddot(n,g,1,vtemp,1)
         if (d.gt.seps*normg2) then
            alpha = normg2/d
         else
            alpha = 2*delta/normg
         endif

!     Compute (h_n - h_cp)'*h_cp, where h_cp = -alpha*g
!     call dcopy(n,hn,1,vtemp,1)
!     call daxpy(n,alpha,g,1,vtemp,1)
!     aa = -alpha*ddot(n,vtemp,1,g,1)
!     c         write (*,*) 'aa = ',aa
!     
!     Try Newton step
!     
         newton = 0
!     if (normh.le.delta .and. aa.gt.seps) then
         if (normh.le.delta) then
            call dcopy(n,hn,1,h,1)
            how = 'Newton'
            newton = 1
         elseif (alpha*normg.lt.delta) then
            aa = alpha*alpha * normg2
            bb = normhn*normhn
            ab = -alpha*ddot(n,g,1,hn,1)
            d = delta*delta
            theta = ab - aa + dsqrt( (ab-d) *  &
     &           (ab-d) + (bb-d) * (d-aa))
!     write (*,*) 'theta = ',theta
            if (abs(theta).lt.abs(d-aa)*seps) then
               theta = 1d0
            else
               theta = (d-aa) / theta
            endif
            do i=1,n
               h(i) = -(1.0 - theta)*alpha*g(i)+theta*hn(i)
            enddo
            normh = dnrm2(n,h,1)
            how = 'mixed'
         else
            do i=1,n
               h(i) = -(delta/normg)*g(i)
            enddo
            normh = delta
            how = 'gradient'
         endif

!     Evaluate function at new point
         do i=1,n
            x1(i) = x(i)+h(i)
         enddo
!     Calculate F(x1)
         f1 = func(n,x1, nthts, nk, nnu, pow, fnu, nrdtot)
         nf = nf + 1
         dF = f1-f
         gh = ddot(n,g,1,h,1)

         if (dF.ge.zero) then
            delta = min(normhn/1.01d0,delta/(3d0*(1+nfail)))
         endif
         if (delta.lt.max(eps,eps1/10*normx)) then
            istop=1
         endif
!     
!     Try parabolic and cubic extrapolation/interpolation steps
!     if they look promising.
!     
         extrap = .false.
         if (iextrap.eq.1 .and. ng.gt.5 .and. f1.ne.(f+gh) &
     &        .and. f1.ne.f) then
            alpha = -0.5*gh/(f1-f-gh)
            if ((alpha.gt.0.1 .and. alpha.lt.0.90) .or. &
     &           (alpha.gt.1.5 .and. alpha.lt.20) ) then
               do i=1,n
                  x2(i) = x(i)+alpha*h(i)
               enddo
               f2 = func(n,x2,nthts,nk,nnu,pow, fnu, nrdtot)
               nf = nf + 1
               if (iprint.eq.1) then
                  write (6,*) 'ratio2,alpha = ',(f-f2)/(f-f1), &
     &                 alpha
                  write (6,*) 'f2 = ',f2,(x2(i),i=1,n)
               endif
               if (f2.lt.f) then   
                  alpha2 = cubic3(zero,f,gh,one,f1,alpha,f2,ierr)

                  if (ierr.eq.1) then
                     alpha2 = extrap3p(zero,one,alpha,f,f1,f2)
                     if (iprint.eq.1)  write (*,*) 'Quadratic3'
                  else
                     if (iprint.eq.1)  write (*,*) 'Cubic3'
                  endif
                  if (abs(alpha-alpha2) .gt. 0.3*abs(alpha-1.0) .and. &
     &                 alpha2.gt.0.1 .and. alpha2 .lt.20) then
                     do i=1,n
                        x3(i) = x(i)+alpha2*h(i)
                     enddo
                     f3 = func(n,x3,nthts,nk,nnu,pow, fnu, nrdtot)
                     nf = nf + 1
                     if (iprint.eq.1) then
                        write (6,*) 'ratio3,alpha2 = ',(f-f3)/(f-f1), &
     &                       alpha2
                        write (6,*) 'f3 = ',f3,(x3(i),i=1,n)
                     endif
                     if (f3.lt.min(f1,f2)) then
                        f2=f3
                        call dcopy(n,x3,1,x2,1)
                        alpha = alpha2
                     endif
                  endif
                  if (f2.lt.f1) then
                     f1 = f2
                     call dcopy(n,x2,1,x1,1)
                     dF = f1-f
                     do i=1,n
                        h(i) = x1(i)-x(i)
                     enddo
                     gh = alpha*gh
                     delta = max(0.9*alpha*normh,delta)
                     normh = alpha*normh
                     extrap = .true.
                  endif
               endif
            endif
         endif

!     
!     Update trust region and Hessian approximation
!     
         if (istop.ne.1 .and. dF.lt.0d0) then

            if (.not.extrap) then
               if (normhn.gt.delta) then
!     compare with linear model
                  dL = gh
               else
!     compare with quadratic model
                  call dspmv('u',n,1D0,B,h,1,0d0,vtemp,1)
                  dL = gh - 0.5*ddot(n,h,1,vtemp,1)
               endif
               if (dL.eq.0.0) then
!     avoid division by zero
                  rho = 0.2
               else
                  rho = dF/dL
               endif
               if (rho.gt.0.75) then
                  delta = min(1.5*delta,100*(normx+1))
               else if (rho.le.0.25) then
                  delta = max(delta/2,2*seps*(normx+1))
               endif
            endif
            
!     Update stopping criterion
            istop = 1
            if (2*max(abs(f-f1),-gh).gt.eps2*(abs(f1)+abs(f))) then
               istop = 0        
            endif
            do i=1,n
               if (2*abs(h(i)).gt.eps1*(eps1+abs(x(i)))) then
                  istop = 0
               endif
            enddo


            
!     j = 1
!     do i=1,n
!     write (*,*) 'B(',i,',',i,') = ',B(j) 
!     write (*,*) 'HH(',i,',',i,') = ',HH(j) 
!     j = j+(i+1)
!     enddo  


!     
!     Calculate gradient in x1
!     
            if (istop.ne.1) then
               if (igrad.eq.1) then
                  call grad(n,x1,g1,nthts,nk,nnu,pow,fnu,nrdtot)
                  ng = ng + 1
               else
                  if (normh.lt.fdtol*normx) then
!     write (*,*) 'CENTRAL'
                     call cdgrad(func,cdstep,n,x1,f1,g1,nthts,nk,nnu,pow, fnu, nrdtot)
                     ng = ng + 1
                     nf = nf + 2*n
                  else
!     write (*,*) 'FORWARD'
                     call fdgrad(func,fdstep,n,x1,f1,g1,nthts,nk,nnu,pow, fnu, nrdtot)
                     ng = ng + 1
                     nf = nf + n
                  endif
               endif
               
               do i=1,n
                  g(i) = g1(i) - g(i)
               enddo
               normy = dnrm2(n,g,1)
               alpha = ddot(n,g,1,h,1)

!     If the curvature condition is fulfilled then
!     update Hessian B and its inverse HH
               if (alpha.gt.seps*normy*normh) then
!     
!     Update B via the BFGS updating formula.
                  alpha = 1d0/alpha
                  call dspmv('u',n,1D0,B,h,1,0d0,vtemp,1)
                  call dspr('u',n,alpha,g,1,B)
                  aa = 1D0/ddot(n,h,1,vtemp,1)
                  call dspr('u',n, -aa ,vtemp,1,B)
!     
!     Update HH via the INVERSE BFGS updating formula.
                  call dspmv('u',n,1D0,HH,g,1,0d0,vtemp,1)
                  call dspr2('u',n,-alpha,vtemp,1,h,1,HH)
                  aa = alpha*(1.0+alpha*ddot(n,g,1,vtemp,1))
                  call dspr('u',n,aa,h,1,HH)
               endif
               
               if (delta.lt.eps1/10*normx) then
                  if (iprint.eq.1) then
                     write (6,*) 'Stop due to small delta:'
                     write (6,*) '  delta = ',delta
                  endif
                  istop = 1
               endif
            endif           
            call dcopy(n,x1,1,x,1)
            f = f1
            call dcopy(n,g1,1,g,1)
            
            normg = dnrm2(n,g,1)
            normg2 = normg*normg
            normx = dnrm2(n,x,1)
         elseif (nfail.eq.3) then

            nfail = 0
!     Step was uphill - restart with Hessian = I.
            if (iprint.eq.1)    write (*,*) 'RESTARTING'
            do i=1,(n*(n+1))/2
               HH(i) = 0D0      
               B(i) = 0D0
            enddo
            
!     HH0 = I / normg, B0 = normg * I
!     HH0 = I, B0 = I
            j = 1
            do i=1,n
               HH(j) = 1D0
               B(j) = 1D0
               j = j+(i+1)
            enddo  
         else
            nfail = nfail+1
         endif
         if (iprint.eq.1) write (*,100) nf,ng,f,normg,delta, &
     &        dF/f,how
         if (iprint.eq.1) write (*,*) 'X = ',(x(i),i=1,n)
         if (iprint.eq.1) write (*,*) 'G = ',(g(i),i=1,n)
         k = k+1
      enddo
      if (nf.ge.nfmax) then
        if(verbose == 1)then
         write (6,*) 'Warning: Maximum # of function evaluations' &
     &        // ' exceeded'
        endif
      endif
      
 90   if (iprint.eq.1)      write(6,*) 'EXITING DOGLEG'
      if (iprint.eq.1)      call flush(6)
      ioptions(1) = nf
      return
 100  format(2i5,1p4e15.6,'  ',a10)
      end


      subroutine fdgrad(func,fdstep,n,x,f,g,nthts,nk,nnu,pow, fnu, nrdtot)
!     
!     Gradient calculation using finite differences (Forward difference).
!     
      implicit none
      integer n
      double precision x(*),f,g(*),f1,x0
      integer i, nthts, nk, nnu, nrdtot
      external func
      double precision func,temp,delta,fdstep
      integer, parameter :: double=8
      real(kind=double),  dimension(nk,nrdtot) :: fnu
      real(kind=double),  dimension(nthts,nk,nnu) :: pow

      
!     fdstep = 7.45058059D-9
!     fdstep = 1d-4
!     write(*,*) 'ENTER fdgrad'
      do i=1,n
         x0 = x(i)
         delta = min(fdstep*max(abs(x0),1d00),0.1d0)
         temp = x0 + delta
         delta = temp - x0
         x(i) = x0 + delta
         f1 = func(n,x,nthts,nk,nnu,pow,fnu,nrdtot)
         g(i) = (f1 - f)/delta
         x(i) = x0
      enddo
!     write(*,*) 'EXIT fdgrad'
      return
      end


      subroutine cdgrad1(func,cdstep,n,x,f,g,diagB,nthts,nk,nnu,pow, fnu, nrdtot)
!     
!     Gradient calculation using finite differences (Central difference).
!     
      implicit none
      integer n
      double precision x(*),f,g(*),f1,f2,diagB(*)
      integer i, nthts,nk,nnu,nrdtot
      integer, parameter :: double=8
      real(kind=double),  dimension(nk,nrdtot) :: fnu
      real(kind=double),  dimension(nthts,nk,nnu) :: pow

      external func
      double precision func,temp,delta, delta1,delta2,cdstep,q,q2,x0

!     cdstep = 6.05545445D-6
!     cdstep = 5d-3
!     write(*,*) 'ENTER cdgrad'
      do i=1,n
         x0 = x(i)
         delta = min(cdstep*max(abs(x0),1.0d0),0.1d0)
         temp = x0 + delta
         delta1 = temp - x0
         x(i) = x0 + delta1
         f1 = func(n,x,nthts,nk,nnu,pow,fnu,nrdtot)

         temp = x0 - delta
         delta2 = x0 - temp
         x(i) = x0 - delta2
         f2 = func(n,x,nthts,nk,nnu,pow,fnu,nrdtot)
         x(i) = x0
         
         q = delta1/delta2
         q2 = q*q

         g(i) = ((f1-f2) - (1d0-q2)*(f-f2))/(delta1*(1d0+q))
         diagB(i) = (f1 -2*f + f2)/delta**2
      enddo
!     write(*,*) 'EXIT cdgrad'
      return
      end


      subroutine cdgrad(func,cdstep,n,x,f,g,nthts,nk,nnu,pow,fnu,nrdtot)
!     
!     Gradient calculation using finite differences (Central difference).
!     
      implicit none
      integer n
      double precision x(*),f,g(*),f1,f2
      integer i,nthts,nk,nnu,nrdtot
      external func
      double precision func,temp,delta, delta1,delta2,cdstep,q,q2,x0
      integer, parameter :: double=8
      real(kind=double),  dimension(nk,nrdtot) :: fnu
      real(kind=double),  dimension(nthts,nk,nnu) :: pow


!     cdstep = 6.05545445D-6
!     cdstep = 5d-3
!     write(*,*) 'ENTER cdgrad'
      do i=1,n
         x0 = x(i)
         delta = min(cdstep*max(abs(x0),1.0d0),0.1d0)
         temp = x0 + delta
         delta1 = temp - x0
         x(i) = x0 + delta1
         f1 = func(n,x,nthts,nk,nnu,pow,fnu,nrdtot)

         temp = x0 - delta
         delta2 = x0 - temp
         x(i) = x0 - delta2
         f2 = func(n,x,nthts,nk,nnu,pow,fnu,nrdtot)
         x(i) = x0
         
         q = delta1/delta2
         q2 = q*q

         g(i) = ((f1-f2) - (1d0-q2)*(f-f2))/(delta1*(1d0+q))
      enddo
!     write(*,*) 'EXIT cdgrad'
      return
      end




      double precision function extrap3p(xa,xb,xc,fxa,fxb,fxc)
!     
!     Find the minimum of the parabola going through the points
!     (a,fa), (b,fb) and (c,fc)
!     
      implicit none
      double precision a,b,c,fa,fb,fc
      double precision xa,xb,xc,fxa,fxb,fxc
      integer i
      double precision tmp,alpha,r,q

      a=xa
      b=xb
      c=xc
      fa=fxa
      fb=fxb
      fc=fxc
!     write (6,*) 'a,b,c = ',a,b,c
!     write (6,*) 'fb-fa,fc-fa = ',fb-fa,fc-fa
      do i=1,2
         if (a.gt.b) then
            tmp=a
            a=b
            b=tmp
            tmp=fa
            fa=fb
            fb=tmp
         endif
         if (b.gt.c) then
            tmp=c
            c=b
            b=tmp
            tmp=fc
            fc=fb
            fb=tmp
         endif
      enddo
!     write (6,*) 'a,b,c = ',a,b,c
!     write (6,*) 'fb-fa,fc-fa = ',fb-fa,fc-fa

      if (fb .lt. fa+(b-a)*(fc-fa)/(c-a)) then
         r = (b-a)*(fb-fc)
         q = (b-c)*(fb-fa)
         alpha = b - 0.5d0*((b-c)*q-(b-a)*r)/sign(max(abs(q-r),1d-16), &
     &        q-r)
      else
         alpha = 0
      endif
      extrap3p = alpha
      end



      double precision function cubic3(a,fa,dfa,b,fb,c,fc,ierr)

!     
!     Find the minimum of the cubic polynomial p(x) which satifies
!     p(a)=fa, p'(a)=dfa, p(b)=fb and p(c)=fc
!     
!     Points must satisfy a.ne.b and a.ne.c.

      implicit none
      double precision a,fa,dfa,b,fb,c,fc
      integer ierr
      double precision b2,b3,c2,c3,det,rb,rc,z1,z2,d
      
      if (b.eq.c) then
!     Only two distinct points: Use quadratic instead of cubic 
!     interpolation.
         if (fb.ne.(fa+dfa*(b-a))) then
            cubic3 = -0.5*dfa*(b-a)/(fa-fb+dfa*(b-a))
            ierr = 0
         else
            cubic3 = 0d0
            ierr = 1
            return
         endif
      endif

      if (b.lt.c) then
         b2 = b*b
         b3 = b2*b
         rb = fb-fa-dfa*b
         c2 = c*c
         c3 = c2*c
         rc = fc-fa-dfa*c
      else
         c2 = b*b
         c3 = b2*b
         rc = fb-fa-dfa*b
         b2 = c*c
         b3 = c2*c
         rb = fc-fa-dfa*c
      endif

!     Solve the 2-by-2 linear system 
!     
!     ( b2  b3 ) ( z1 )  =  ( rb )
!     ( c2  c3 ) ( z2 )  =  ( rc )

      if (det .eq. 0d0) then
         cubic3 = 0d0
         ierr = 1
         return
      endif
      det = b2*c3-b3*c2
      z1 = (c3*rb-b3*rc)/det
      z2 = (b2*rc-c2*rb)/det

      d = z1*z1 - 3*z2*dfa
      if (d .lt. 0d0) then
         cubic3 = 0d0
         ierr = 1
         return
      endif
      cubic3 =  (-z1 + sqrt(d))/(3*z2)
      ierr = 0
      return
      end

