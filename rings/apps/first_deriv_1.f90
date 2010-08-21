    module first_deriv_1
      contains
      SUBROUTINE DF1(nparm,pr,G,nthts,nk,nnu,pow,fnu,nrdtot)
      use ring_pass
      implicit none
      integer, parameter :: double=8
      real(kind=double), dimension(nparm) :: pr, G
      integer :: ik,iendnu,istrtnu,inu,itht,nparm,ier,ierr,i
      integer :: nthts, nk, nnu, nrdtot
      real(kind=double) :: den,den1,den2,dent,dent2,frng,nu1,akt,pmd,bkgt1
      real(kind=double) :: twopi,ga,pp,ch,dch,theta,dtht,gsun,fmode
      real(kind=double), dimension(nthts) :: uk,c1,c2, cosine, sine
      real(kind=double), dimension(nthts,nk,nnu) :: pow
      real(kind=double), dimension(nk,nrdtot) :: fnu

      G = 0d0
!      if (v == 1) then
!        write(*,'(" df1,g = ",6es15.4)')(G(i),i=1,6)
!      endif 
 
      twopi = 2d0*dacos(-1d0)      
      ga = exp(pr(5))
      den1 = ga*ga/4d0
      pmd = exp(pr(3))*ga/2.D0
      dtht=twopi/dble(nthts)
      do itht=1,nthts
         theta=(itht-1)*dtht
         cosine(itht) = dcos(theta)
         sine(itht) = dsin(theta)
      enddo
      do ik = istk,iendk
         akt  =  dble(ik)*dk          
         bkgt1  =  exp(pr(4))/(akt**3)
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
!           write(*,'("in first_deriv_1: istrtnu,iendnu",2(i6,2x))')istrtnu,iendnu

         do itht  =  1,nthts
            c1(itht)  =  (akt/twopi)*cosine(itht)
            c2(itht)  =  (akt/twopi)*sine(itht)
            uk(itht)  =  -fnu(ik,nj) - pr(6) + pr(1)*c1(itht) &
     &           + pr(2)*c2(itht) 
         enddo
!
!  Calculate first derivatives
!
!         if(v == 1) then
!           write(*,'("df1, frng = ", es15.4)') frng
!           write(*,'("df1, istrtnu, iendnu ", i6,3x,i6)')istrtnu, iendnu
!           write(*,'("df1, dnu, dk", 2es15.4)')dk, dnu
!         endif
         do inu = istrtnu,iendnu
            nu1 = (inu-1)*dnu
            do itht = 1,nthts
               if (pow(itht,ik,inu) .gt. 0d0) then
                  den = nu1 + uk(itht)
                  den2 = den*den
                  dent = den2 + den1
                  dent2 = dent*dent
                  pp = pmd/dent + bkgt1
                  ch = 1d0/(pp*pow(itht,ik,inu))
                  dch = (1d0-ch)/pp
                  G(1) = G(1) + dch*( -2*c1(itht)*pmd*den/dent2 )
                  G(2) = G(2) + dch*( -2*c2(itht)*pmd*den/dent2 )
                  G(3) = G(3) + dch*( pmd/dent )
                  G(4) = G(4) + dch*( bkgt1 )
                  G(5) = G(5) + dch*( pmd/dent - 2*pmd*den1/dent2 )
                  G(6) = G(6) + dch*( 2*pmd*den/dent2 )               
               endif
            enddo
         enddo
      enddo
!      if (v == 1) then
!        write(*,'(" df1, g = ",6es15.4)')(g(i),i=1,6)
!      endif 
 9999 return
      end subroutine
    end module
