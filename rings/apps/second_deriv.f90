    module second_deriv 
      contains
      SUBROUTINE DF2(npar,pr,flag,nthts,nk,nnu,pow,fnu,nrdtot,s2)
 
!     Changed to be able to use
!     first derivative x first derivative hessian (no second derivs)
!     on Jan. 10, 2007, if flag = 0 use 2nd derivative hessian, if flag = 1 just
!     use first derivatives
!
      use ring_pass
      implicit none
      integer, parameter :: double = 8
      integer :: ik,iendnu,istrtnu,inu,itht,i,j,flag,npar,ier,ierr
      integer :: nnu, nthts, nk, nrdtot
      real(kind=double) :: dpd2, ppsqinv,dtht,theta
      real(kind=double), dimension(npar) :: dpr
      real(kind=double), dimension(npar) :: pr
      real(kind=double), dimension(npar,npar) :: s2
      real(kind=double) :: den,den1,den2,dent,dent2,dpdg,dpdnu,dch2
      real(kind=double) :: frng,nu1,akt,pmd,bkgt1,wg2,wg3, gsun, fmode
      real(kind=double) :: twopi,ga,pp,ch,dch,dpdu,dpdv,dpdp,dpdb
      real(kind=double), dimension(nthts) :: uk,c1,c2
      real(kind=double), dimension(nthts) :: cosine, sine
      real(kind=double), dimension(nthts,nk,nnu) :: pow
      real(kind=double), dimension(nk,nrdtot) :: fnu
      s2=0d0
      twopi = 2d0*dacos(-1d0)
      dtht=twopi/dble(nthts)
      ga = exp(pr(5))
      den1 = ga*ga/4d0
      pmd = exp(pr(3))*ga/2D0
      do itht=1,nthts
         theta=(itht-1)*dtht
         cosine(itht) = dcos(theta)
         sine(itht) = dsin(theta)
      enddo
      if (v == 1) then
        write(*,'(" df2, nj = ",i5)')nj
      endif
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
!           write(*,'("in second_deriv: istrtnu,iendnu",2(i6,2x))')istrtnu,iendnu

!      if (v == 1) then
!        write(*,'(" df2, frng = ",es15.5)')frng
!        write(*,'(" df2, dk = ",es15.5)')dk
!        write(*,'(" df2, dnu = ",es15.5)')dnu
!        write(*,'(" df2, istrtnu, iendnu ",2(i5,3x))')istrtnu,iendnu
!        write(*,'(" df2, flag = ",i5)')flag
!      endif
        if ( iendnu > nnu) iendnu = nnu
        do itht = 1,nthts
          c1(itht) = (akt/twopi)*cosine(itht)
          c2(itht) = (akt/twopi)*sine(itht)
          uk(itht) = -fnu(ik,nj) - pr(6) + pr(1)*c1(itht) &
     &           + pr(2)*c2(itht)
        enddo
 
         do inu = istrtnu,iendnu
            nu1 = (inu-1)*dnu
            do  itht = 1,nthts
               if (pow(itht,ik,inu) .gt. 0d0) then
                  den = nu1 + uk(itht)
                  den2 = den*den
                  dent = den2 + den1
                  dent2 = dent*dent
                  pp = pmd/dent + bkgt1
                  ch = 1d0/(pp*pow(itht,ik,inu))
                  dch = (1d0-ch)/pp
                  dch2 = (1d0/pp - 2*dch)/pp
                  dpdu = -2*c1(itht)*pmd*den/dent2
                  dpdv = -2*c2(itht)*pmd*den/dent2
                  dpdp = pmd/dent
                  dpdb = bkgt1
                  dpdg = pmd/dent - 2*pmd*den1/dent2
                  dpdnu = 2*pmd*den/dent2
                  dpd2 = (2*den*den)/(dent2*dent) - 2d0/dent2
                  dpr(1) = dpdu
                  dpr(2) = dpdv
                  dpr(3) = dpdp
                  dpr(4) = dpdb
                  dpr(5) = dpdg
                  dpr(6) = dpdnu
                  wg2 = pmd/dent2
                  wg3 = (pmd*den)/(dent2*dent)
                  ppsqinv = 1d0/pp*pp
                  if (flag .eq. 0) then 
                   s2(1,1) = s2(1,1) + dch2*dpdu*dpdu + &
     &                 dch*dpd2*pmd*c1(itht)*c1(itht)
                   s2(2,2) = s2(2,2) + dch2*dpdv*dpdv + &
     &                 dch*dpd2*pmd*c2(itht)*c2(itht)
                   s2(3,3) = s2(3,3) + dch2*dpdp*dpdp + dch*dpdp
                   s2(4,4) = s2(4,4) + dch2*dpdb*dpdb + dch*dpdb
                   s2(5,5) = s2(5,5) + dch2*dpdg*dpdg &
     &               + dch*(dpdp-2.*wg2*ga**2+.5*wg3*ga**4/den)
                   s2(6,6) = s2(6,6) + dch2*dpdnu*dpdnu &
     &               + dch*2.*(-wg2+4.*den*wg3)
                   s2(2,1) = s2(2,1) + dch2*dpdu*dpdv + &
     &                 dch*dpd2*pmd*c1(itht)*c2(itht)
                   s2(1,2) = s2(2,1)
                   s2(3,1) = s2(3,1) + dch2*dpdu*dpdp + dch*dpdu
                   s2(1,3) = s2(3,1)
                   s2(4,1) = s2(4,1) + dch2*dpdu*dpdb
                   s2(1,4) = s2(4,1)
                   s2(5,1) = s2(5,1) + dch2*dpdu*dpdg &
     &                + dch*2*c1(itht)*((wg3*ga**2)-(wg2*den))
                   s2(1,5) = s2(5,1)
                   s2(6,1) = s2(6,1) + dch2*dpdu*dpdnu &
     &                + dch*(2*c1(itht))*(wg2-4.*wg3*den)
                   s2(1,6) = s2(6,1)
                   s2(3,2) = s2(3,2) + dch2*dpdv*dpdp + dch*dpdv
                   s2(2,3) = s2(3,2)
                   s2(4,2) = s2(4,2) + dch2*dpdv*dpdb
                   s2(2,4) = s2(4,2)
                   s2(5,2) = s2(5,2) + dch2*dpdv*dpdg &
     &               + dch*2*c2(itht)*((wg3*ga**2)-(wg2*den))
                   s2(2,5) = s2(5,2)
                   s2(6,2) = s2(6,2) + dch2*dpdv*dpdnu &
     &                + dch*(2*c2(itht))*(wg2-4.*wg3*den)
                   s2(2,6) = s2(6,2)
                   s2(4,3) = s2(4,3) + dch2*dpdp*dpdb
                   s2(3,4) = s2(4,3)
                   s2(5,3) = s2(5,3) + dch2*dpdp*dpdg + dch*dpdg
                   s2(3,5) = s2(5,3)
                   s2(6,3) = s2(6,3) + dch2*dpdp*dpdnu + dch*dpdnu
                   s2(3,6) = s2(6,3)
                   s2(5,4) = s2(5,4) + dch2*dpdb*dpdg
                   s2(4,5) = s2(5,4)
                   s2(6,4) = s2(6,4) + dch2*dpdb*dpdnu
                   s2(4,6) = s2(6,4)
                   s2(6,5) = s2(6,5) + dch2*dpdg*dpdnu &
     &                + dch*2.*(wg2*den - wg3*ga**2)
                   s2(5,6) = s2(6,5)
                  else 
                   do i = 1,npar
                     do j = 1,npar                    
                      s2(j,i) = s2(j,i) + dpr(j)*dpr(i)*ppsqinv
                     enddo
                   enddo
                  endif
               endif
            enddo
         enddo
      enddo 
     
 9999 return
      end subroutine
     end module
