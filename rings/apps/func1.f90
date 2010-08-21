    module func1
      integer, parameter :: d_1=8
      contains
      REAL(kind=d_1) FUNCTION FUNK1(npar,pr,nthts,nk,nnu,pow,fnu,nrdtot)
      use ring_pass
      implicit none
      real(kind=d_1), dimension(npar) :: pr
      integer ik,iendnu,istrtnu,inu,itht,npar,ier,ierr,i
      integer nk,nnu,nthts,nrdtot
      real(kind=d_1) :: frng,nu1,chunk,akt,x,den,pmd,bkgt1,rchi
      real(kind=d_1) :: dtht,theta
      real(kind=d_1) :: twopi,xchi,ukx,uky,ga,tmp,den1,Enum,Denom
      real(kind=d_1), dimension(nnu) :: Prod 
      real(kind=d_1), dimension(nthts) :: uk, cosine, sine
      real(kind=d_1), dimension(nthts,nk,nnu) :: pow
      real(kind=d_1), dimension(nk,nrdtot) :: fnu
      real(kind=d_1) :: safemin,safemax,fmode,gsun
      external dlamch
      real(kind=d_1) :: dlamch
 
!     Determine underflow and overflow limits and allow a safety margin
!     (4d16) to guarantee accurate results even when the computer does 
!     not implement gradual underflow. This is e.g. the case on Alpha and 
!     Mips based machines. 
!      write(*,'("in func1: ni,nj,npj,ntot,nnu",5(i6,2x))')ni,nj,npj,ntot,nnu 
!       write(*,'("in func1: dk, dnu, pr",8(1pd13.5,2x))')dk,dnu,(pr(i),i=1,6)
!       write(*,'("in func1: istk,iendk",2(i5,2x))')istk, iendk
!       write(*,'("in func1: fnu",8(1pd13.5,2x))')((fnu(ik,i),ik=istk,iendk), &
!                  &i=ni,nj)
!       write(*,'("in func1: pow ",8(1pd13.5,2x))')(((pow(itht,ik,inu), &
!                  &itht=1,8),ik=15,20),inu=300,302)
      safemin = dlamch('s')*4d16
      safemax = 1d0/safemin

      Prod = 0d0
!     Set various constants
      ga = exp(pr(5))
      den1 =  ga*ga/4.0d0
      pmd = exp(pr(3))*ga/2.0D0
      twopi = 2d0*dacos(-1d0)
      chunk = 0d0
      xchi = 0d0
      npj = 0
      dtht=twopi/dble(nthts)

      do itht=1,nthts
         theta=(itht-1)*dtht
         cosine(itht) = dcos(theta)
         sine(itht) = dsin(theta)
      enddo
!
! Do the actual calculation of the likelihood function  by
!  summing up the contributions for each k, omega and theta.
!  More than one k is only needed when ikm > 0, i.e. you're 
!  smoothing in k.
! 
      do ik = istk,iendk
         akt = dble(ik)*dk          
         bkgt1 = exp(pr(4))/(akt**3)
         if ((nj+1 <= nrdtot) .and. (fnu(ik,nj+1) > 5d0)) then
           frng = 0.40d0 * (fnu(ik,nj+1) - fnu(ik,nj)) / dnu
         else if (nj+1 > nrdtot) then
           frng = 0.35d0 * (fnu(ik,nj) - fnu(ik,nj-1)) / dnu
         endif
         if ((ni /= 0) .and. (fnu(ik,nj+1) < 5d0)) then
           frng = 0.35d0 * (fnu(ik,nj) - fnu(ik,nj-1)) / dnu
         endif
         istrtnu = nint(fnu(istk,nj) / dnu - frng)
         iendnu = nint(fnu(iendk,nj) / dnu + frng)

         if (istrtnu > iendnu) then
           frng = 0.35d0 * (fnu(ik,nj) - fnu(ik,nj-1)) / dnu
           istrtnu = nint(fnu(istk,nj) / dnu - frng)
           iendnu = nint(fnu(iendk,nj) / dnu + frng)
         endif

!
!if starting and stopping end values still inappropriate due to
!extension of guess table, use asymptotic relationship to give
!it a try
!
          gsun=2.74d8 !(Mm/microsec^2)
          if (istrtnu > iendnu) then
            fmode=dble(sqrt(gsun*akt)) / twopi
            frng = 0.4d0*fmode*dble(sqrt(real(nj)+1.5)-sqrt(real(nj)+0.5)) / dnu
            istrtnu = nint(fnu(istk,nj) / dnu - frng)
            iendnu = nint(fnu(iendk,nj) / dnu + frng)
            write(*,'("bad strt and stop frequency bins")')
            write(*,'("recalculated using asymptotic formula",2(i6,2x))') &
                   &   istrtnu,iendnu
          endif

         if (iendnu > nnu) iendnu = nnu
         do itht = 1,nthts
            ukx = pr(1)*akt/twopi*cosine(itht)
            uky = pr(2)*akt/twopi*sine(itht)
            uk(itht) = ukx+uky-fnu(ik,nj)-pr(6)
         enddo

         do inu = istrtnu,iendnu
            nu1 = (inu-1)*dnu
            Enum = 0d0
            Denom = 1d0
            do  itht = 1,nthts
!
!     Only use points with positive power. 
!     Spots of negative power are an artifact from the 
!     cubic interpolation and filtering process.
!
               if (pow(itht,ik,inu) .gt. 0d0) then
                  den = nu1 + uk(itht)
                  den = den*den + den1
                  x = pmd/den + bkgt1
                  x = x*pow(itht,ik,inu)
                  tmp = (1d0 - x)
                  xchi = xchi + tmp*tmp
                  Enum = x*Enum + Denom
                  Denom = x*Denom
               endif
            enddo
 
            if (Denom.lt.safemin .or. Denom.gt.safemax .or. &
     &          Enum.lt.safemin .or. Enum.gt.safemax  ) then
!
!     Over- or underflow probably occured. Use old reliable
!     and slow method.
!
               do itht = 1,nthts 
                  if (pow(itht,ik,inu) .gt. 0d0) then
                     den = nu1 + uk(itht)
                     den = den*den + den1
                     x = pmd/den + bkgt1
                     x = x*pow(itht,ik,inu)
                     tmp = (1d0 - x)
                     xchi = xchi + tmp*tmp
                     chunk = chunk + log(x) + 1d0/x
                  endif
               enddo
               Prod(inu) = 1d0
            else
               Prod(inu) = Denom
               chunk = chunk +  Enum/Denom
            endif
         enddo
!         chunk = sum(log(Prod(istrtnu:iendnu)))
         do inu = istrtnu,iendnu
            chunk = chunk + log(Prod(inu))
         enddo
         npj = npj + (iendnu-istrtnu+1)*nthts
 120  end do 
      funk1 = chunk
      rchi = xchi
 9999 return
      end function
    end module
