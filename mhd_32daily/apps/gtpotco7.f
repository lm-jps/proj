**************************************************
* GTPOTCO7.FOR
* 
* INPUT  : stnRRRR.dat  : raw Stanford data
*          stncrk.lst   : list of file containing brank of data
* OUTPUT : potRRRRc.LLL : potential coefficients
*                                                 by K.Hayashi
**************************************************
       program gtpotco
       implicit none
       integer nrots, nrote, nrot, i, nng
       integer lons, lone, dlon, lon
       integer ngrt(100)
       character*30 cdummy
       logical loldname


       do i = 1, 100
         ngrt(i) = -1000
       enddo
* for JSOC-purpose, this is not needed
!       open(unit = 3, file = 'stncrk3.lst', status = 'old')
!       nng = 0
! 200   continue
!         nng = nng + 1
!         read(3,'(i4,1x,a)', END = 299) ngrt(nng), cdummy
!         goto 200
! 299   continue

       write(*,'('' input start CRN and Longitude : '',$)')
       read(*,*) nrots, lons
       write(*,'('' input end   CRN and Longitude : '',$)')
       read(*,*) nrote, lone
       write(*,'('' input step of  Longitude : '')')
       write(*,'(''       0 or negative ... old_style : '',$)')
       read(*,*) dlon

       if (dlon .LE. 0) then
         lons = 180
         lone = 180
         loldname = .true.
       else
         loldname = .false.
       endif

       nrot = nrots
       lon  = lons
       do while (nrot .LE. nrote)
         write(*,'('' Now CRN Lon = '',2i5)') nrot, lon
         call getpotcoef(nrot,lon,ngrt,loldname)
         if (dlon .GT. 0) then
           lon = lon - dlon
           if (lon .LE. 0) then
             lon = lon + 360
             nrot = nrot + 1
           endif
         else
           lon = 180
           nrot = nrot + 1
         endif
         if ((nrot .EQ. nrote) .AND. (lon .LT. lone)) nrot = nrot + 1
       enddo

       stop
       end

************** END of MAIN PROGRAM **************************  

************** SUB PROGRAMS *********************************
*
       subroutine getpotcoef(nrot,lon,ngrt,loldname)
       implicit none
*
       integer nrot, lon
       integer ngrt(100)
       logical loldname
*
       real    br(72, 30), brb(216,30)
       real    pota(11,11), potb(11,11)
       real    pnmx
       real    rw, pi, domeg
       logical lfilled
       integer i, j, n, m, k, k1, k2, dncr
       character*48 flname1, flname2
       real    sumbr, sumbr2, pnmx2
       real    mphi, th, xx
       real    suma, sumb, dsuma, dsumb, nmfactor, an, cn
       real    theta(30), phi(72), zz(30)
       parameter(rw = 2.5, pi = 3.1415926535)
       parameter(domeg = 4. * pi / 30. / 72.)

** preparing : location of data (i, j)
* theta and z ! latitude
       do i = 1, 30
         zz(i) = float(31 - 2 * i) / 30.0
         theta(i) = acos(zz(i))
       enddo
* phi         ! longitude
       do i = 1, 72
         phi(i) = (360.0 - 5.0 * float(i - 1)) * pi / 180.00
       enddo

** open raw data file and read synoptic chart of solar magnetic field
       do dncr = 0, 2
         write(flname1,'(''stn'',i4,''.dat'')') nrot + dncr - 1
         open(unit = 1, file = flname1 ,status = 'old', ERR=999)
         write(*,*) ' Input Data = ', flname1
         do i = 1, 72 
           read(1,10) (brb(i+dncr*72, j), j =  1,  6)
           read(1,20) (brb(i+dncr*72, j), j =  7, 14)
           read(1,20) (brb(i+dncr*72, j), j = 15, 22)
           read(1,20) (brb(i+dncr*72, j), j = 23, 30)
         enddo
         close(1)
 999     continue
       enddo
 10    format(18x,6f9.3)
 20    format(8f9.3)

* make sure stnXXXX.dat were correctly loaded...
       open(unit=11,file='stndummy.dat',status='unknown')
         dncr = 1
         do i = 1, 72 
           write(11,10) (brb(i+dncr*72, j), j =  1,  6)
           write(11,20) (brb(i+dncr*72, j), j =  7, 14)
           write(11,20) (brb(i+dncr*72, j), j = 15, 22)
           write(11,20) (brb(i+dncr*72, j), j = 23, 30)
         enddo
       close(11)

* make a combined(?) map
       k1 = (365 - lon) / 5 +  36
       do k = k1, k1 + 71
         k2 = mod(k,72)
         if (k2 .EQ. 0) k2 = 72
         do j = 1, 30
           br(k2,j) = brb(k,j)
         enddo
       enddo

* save new map
       if (.NOT. loldname) then
         write(flname2,'(''stn'',i4,''.'',i3.3)') nrot, lon
         open(unit = 2, file = flname2 ,status = 'unknown')
         write(*,*) ' OutputData = ', flname2
         do i = 1, 72
           write(2,'(''CT'',i4,'':'',i3.3,''        '',6f9.3)') 
     +       nrot, 365 - i*5, (br(i,j),j=1,6)
           write(2,'(8f9.3)') (br(i,j),j= 7,14)
           write(2,'(8f9.3)') (br(i,j),j=15,22)
           write(2,'(8f9.3)') (br(i,j),j=23,30)
         enddo
         close(2)
       endif

* let sum(magr) = 0
       sumbr = 0.000000
       do j = 1, 30
         th = theta(j)
         do i = 1, 72
           br(i, j) = br(i, j) / sin(th)
           sumbr = sumbr + br(i, j) * domeg
         enddo
       enddo
* confirm sum(magr) = 0
       sumbr2 = 0.000000000
       do j = 1, 30
         do i = 1, 72
           br(i, j) = br(i, j) - sumbr / (4. * pi)
           sumbr2 = sumbr2 + br(i, j) * domeg
         enddo
       enddo
       write(*,'('' Sum Br : '',f11.6,'' => '', f11.6)') sumbr, sumbr2

** integrate on surface and estimate harmonic coefficients
       do n = 1, 10
         do m = 0, n
           suma = 0.000000000
           sumb = 0.000000000
           do j = 1, 30       ! integrate on surface (i-j loop)
             xx = zz(j)
             pnmx2 = pnmx(n, m, xx)
             do i = 1, 72
               mphi = float(m) * phi(i)
               dsuma = br(i, j) * pnmx2 * cos(mphi) * domeg
               dsumb = br(i, j) * pnmx2 * sin(mphi) * domeg
               suma = suma + dsuma
               sumb = sumb + dsumb
             enddo   ! End of j-loop : latitude
           enddo     ! End of i-loop : longitude = End of integrating
           cn = -1.0 / ((rw ** (2 * n + 1)) - 1.)
           an = -1.0 / (cn * float(n) - (1.0 - cn) * float(n + 1))
*           an = 1.0 / (float(n + 1) + float(n) / (rw ** (2 * n + 1)))
           nmfactor = an * float(2 * n + 1) / (4.0 * pi) 
           pota(n + 1, m + 1) = suma * nmfactor
           potb(n + 1, m + 1) = sumb * nmfactor
         enddo       ! End of m-loop
       enddo         ! End of n-loop

** save harmonic coefficients into file
       lfilled = .true.
       do i = 1, 100
         if (nrot .EQ. ngrt(i)) lfilled = .false.
       enddo
       if (lfilled) then
         if (loldname) then
           write(flname2,'(''pot'',i4,''c.dat'')') nrot
         else
           write(flname2,'(''pot'',i4,''c.'',i3.3)') nrot, lon
         endif
       else
         if (loldname) then
           write(flname2,'(''pot'',i4,''d.dat'')') nrot
         else
           write(flname2,'(''pot'',i4,''d.'',i3.3)') nrot, lon
         endif
       endif
       open(unit = 2, file = flname2, status = 'unknown')
       write(*,*) ' OUTPUT DATA = ', flname2

       write(2,'(A)') ' pot_a '
       write(2,'('' m / n '',i7,4i14)') 1, 2, 3, 4, 5
       do i = 1, 11 
         write(2,33) (i - 1), (pota(j, i), j = 2, 6)
       enddo
       write(2,'('' m / n '',i7,4i14)') 6, 7, 8, 9, 10
       do i = 1, 11 
         write(2,33) (i - 1), (pota(j, i), j = 7, 11)
       enddo
       write(2,'(A)') ' pot_b '
       write(2,'('' m / n '',i7,4i14)') 1, 2, 3, 4, 5
       do i = 1, 11 
         write(2,33) (i - 1), (potb(j, i), j = 2, 6)
       enddo
       write(2,'('' m / n '',i7,4i14)') 6, 7, 8, 9, 10
       do i = 1, 11 
         write(2,33) (i - 1), (potb(j, i), j = 7, 11)
       enddo
* 33    format(i3,5e14.7)
 33    format(i3,5(1x,e13.7))
       close(2)
       return
       end


**  The normalized Legendre Polynomials
       real function pnmx(l, m, x)
       implicit none
       integer l, m
       real    x
       real    pmm, pll, pmmp1, somx2, ans, fact
       real    aan, bbn
       integer i, ll
       pmm = 1.
       if (m .GT. 0) then
         somx2 = sqrt((1. - x) * (1. + x))
         fact = 1.
         do i = 1, m
           pmm = pmm * fact * somx2
           fact = fact + 2.
         enddo
       endif
       if (l .EQ. m) then
         ans = pmm
       else
         pmmp1 = x * float(2 * m + 1) * pmm
         if (l .EQ. (m + 1)) then
           ans = pmmp1
         else
           do ll = (m + 2), l
             pll = (x * float(ll * 2 - 1) * pmmp1 - 
     +               float(ll + m - 1) * pmm) / float(ll - m)
             pmm = pmmp1
             pmmp1 = pll
           enddo
           ans = pll
         endif
       endif
* normalize
       if (m .GT. 0) then
         aan = 1.
         bbn = 1.
         do i = (l - m + 1), (l + m)
           aan = aan * float(i)
         enddo
         ans = ans * sqrt(2. * bbn / aan)
       endif
       pnmx = ans
       return
       end

********** end of sub programs ********************

********** end of this program ********************

