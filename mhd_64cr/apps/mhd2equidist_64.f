        program mhdati
        implicit none
* number of grid for each three direction.
        integer ipickup, jpickup, kpickup, i, j, k, imenu
        integer ncr, nstep, ncr2, nstep2
        real*8  var8(8)
        logical lmhd
        real*8  rrpickup, thpickup, phpickup
        real*8  gamma,omega,v0,b0,n0,tmp0,r0
        real*8  pi
        parameter(pi = 3.14159265358979D+00)
        real*8  rupper, rbottom
        parameter(rupper = 5.0E+00, rbottom = 1.0E+00) ! in Rs
        integer ii, jj, kk
        parameter(ii = 72, jj = 64, kk = 128)
        real*8  ro3(0:ii+1,0:jj,0:kk+1),pg3(0:ii+1,0:jj,0:kk+1) ! in
        real*8  ur3(0:ii+1,0:jj,0:kk+1),ut3(0:ii+1,0:jj,0:kk+1)
        real*8  up3(0:ii+1,0:jj,0:kk+1),br3(0:ii+1,0:jj,0:kk+1)
        real*8  bt3(0:ii+1,0:jj,0:kk+1),bp3(0:ii+1,0:jj,0:kk+1)
        integer iout, jout, kout
        parameter(iout = 80 - 1, jout = 60 - 1, kout = 120 - 1)
        real*8  roo(0:iout,0:jout,0:kout),tto(0:iout,0:jout,0:kout) ! output
        real*8  uro(0:iout,0:jout,0:kout),uto(0:iout,0:jout,0:kout)
        real*8  upo(0:iout,0:jout,0:kout),bro(0:iout,0:jout,0:kout)
        real*8  bto(0:iout,0:jout,0:kout),bpo(0:iout,0:jout,0:kout)
        real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
        common  /coord/  rr, theta, phi
        common  /mhdvar/ ro3,pg3,ur3,ut3,up3,br3,bt3,bp3
        common  /mhdout/ roo,tto,uro,uto,upo,bro,bto,bpo
       
        call readinit(gamma,omega,v0,b0,n0,tmp0,r0,rr,theta,phi)

        ncr2   = -1
        nstep2 = -1
        do 111 k = 0, kout
        do 111 j = 0, jout
        do 111 i = 0, iout

          ncr   = 1000
          nstep = 4000

          rrpickup = rbottom
     &             +(float(i) + 0.5E+00) / float(iout+1)
     &             *(rupper-rbottom)
          rrpickup = rrpickup / r0 * 6.96D+10 ! in unit
          thpickup =(float(j) + 0.5E+00) / float(jout+1) * pi
          phpickup =(float(k) + 0.5E+00) / float(kout+1) * pi * 2.0E+00

          if ((ncr .NE. ncr2) .OR. (nstep .NE. nstep2))
     +      call readmhd(ncr,nstep,lmhd)

          call varatrtp(rrpickup, thpickup, phpickup,var8)
          var8(2) = var8(2) / var8(1) * tmp0 / 1.0D+06 ! in MK
          var8(1) = var8(1) * n0 ! in /cc
          var8(3) = var8(3) * v0 / 1.0D+05 ! km/s
          var8(4) = var8(4) * v0 / 1.0D+05 
          var8(5) = var8(5) * v0 / 1.0D+05
          var8(6) = var8(6) * b0 * 1.0D+05 ! BR       =  Bx in nT
          var8(7) = var8(7) * b0 * 1.0D+05 ! Btheta   = -Bz
          var8(8) = var8(8) * b0 * 1.0D+05 ! Bphi     = -By


          roo(i,j,k) = var8(1)
          tto(i,j,k) = var8(2)
          uro(i,j,k) = var8(3)
          uto(i,j,k) = var8(4)
          upo(i,j,k) = var8(5)
          bro(i,j,k) = var8(6)
          bto(i,j,k) = var8(7)
          bpo(i,j,k) = var8(8)

          ncr2 = ncr
          nstep2 = nstep
 111    continue

        call outtext()

        stop
        end

* --------------------------------------------------------------------
*
        subroutine outtext()
        implicit none
        integer iout, jout, kout
        parameter(iout = 80 - 1, jout = 60 - 1, kout = 120 - 1)
        real*8  roo(0:iout,0:jout,0:kout),tto(0:iout,0:jout,0:kout) ! output
        real*8  uro(0:iout,0:jout,0:kout),uto(0:iout,0:jout,0:kout)
        real*8  upo(0:iout,0:jout,0:kout),bro(0:iout,0:jout,0:kout)
        real*8  bto(0:iout,0:jout,0:kout),bpo(0:iout,0:jout,0:kout)
        integer i,j,k
        common  /mhdout/ roo,tto,uro,uto,upo,bro,bto,bpo
       
        open(unit=2,file='d3equidist.dat',status='unknown')
        do k = 0, kout
        do j = 0, jout
        do i = 0, iout
          write(2,'(3i5,5(1x,e11.5),1x,3i5,3(1x,e11.5))')
     &      i+1000, j+2000, k+3000,
     &      roo(i,j,k),tto(i,j,k),
     &      uro(i,j,k),uto(i,j,k),upo(i,j,k),
     &      i+1000, j+2000, k+3000,
     &      bro(i,j,k),bto(i,j,k),bpo(i,j,k)
        enddo
        enddo
        enddo
        close(2)
        return
        end


*
* --------------------------------------------------------------------
*
       subroutine varatrtp(ra,th,ph,var8)
       implicit none
* interface
       real*8  ra, th, ph, var8(8)
* local
       real*8  dra1, dra2, dth1, dth2, dph1, dph2
       integer i1, i2, j1, j2, k1, k2, ijk(3)
       real*8  ans, rtp(3)
*
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
*
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  ro3(0:ii+1,0:jj,0:kk+1),pg3(0:ii+1,0:jj,0:kk+1)
       real*8  ur3(0:ii+1,0:jj,0:kk+1),ut3(0:ii+1,0:jj,0:kk+1)
       real*8  up3(0:ii+1,0:jj,0:kk+1),br3(0:ii+1,0:jj,0:kk+1)
       real*8  bt3(0:ii+1,0:jj,0:kk+1),bp3(0:ii+1,0:jj,0:kk+1)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       common  /coord/  rr, theta, phi
       common  /mhdvar/ ro3,pg3,ur3,ut3,up3,br3,bt3,bp3

       rtp(1) = ra
       rtp(2) = th
       rtp(3) = ph
       if (ph .LT. 0.0D+00) then
         ph = ph + 2.0D+00 * pi
       else if (ph .GT. 2.0D+00 * pi) then
         ph = dmod(ph, 2.0D+00 * pi)
       endif
       rtp(3) = ph
       call rtp2ijk(rtp,ijk)
       i1 = ijk(1)
       j1 = ijk(2)
       k1 = ijk(3)
       i2 = i1 + 1
       j2 = j1 + 1
       k2 = k1 + 1
       dra1 =(ra - rr(i1))/(rr(i2)-rr(i1))
       dra2 =(rr(i2) - ra)/(rr(i2)-rr(i1))
       dth1 =(th - theta(j1))/(theta(j2)-theta(j1))
       dth2 =(theta(j2) - th)/(theta(j2)-theta(j1))
       dph1 =(ph - phi(k1))/(phi(k2)-phi(k1))
       dph2 =(phi(k2) - ph)/(phi(k2)-phi(k1))

       var8(1)=(ro3(i1,j1,k1)*dth2*dph2+ro3(i1,j2,k1)*dth1*dph2
     +        + ro3(i1,j1,k2)*dth2*dph1+ro3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(ro3(i2,j1,k1)*dth2*dph2+ro3(i2,j2,k1)*dth1*dph2
     +        + ro3(i2,j1,k2)*dth2*dph1+ro3(i2,j2,k2)*dth1*dph1)*dra1
       var8(2)=(pg3(i1,j1,k1)*dth2*dph2+pg3(i1,j2,k1)*dth1*dph2
     +        + pg3(i1,j1,k2)*dth2*dph1+pg3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(pg3(i2,j1,k1)*dth2*dph2+pg3(i2,j2,k1)*dth1*dph2
     +        + pg3(i2,j1,k2)*dth2*dph1+pg3(i2,j2,k2)*dth1*dph1)*dra1
       var8(3)=(ur3(i1,j1,k1)*dth2*dph2+ur3(i1,j2,k1)*dth1*dph2
     +        + ur3(i1,j1,k2)*dth2*dph1+ur3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(ur3(i2,j1,k1)*dth2*dph2+ur3(i2,j2,k1)*dth1*dph2
     +        + ur3(i2,j1,k2)*dth2*dph1+ur3(i2,j2,k2)*dth1*dph1)*dra1
       var8(4)=(ut3(i1,j1,k1)*dth2*dph2+ut3(i1,j2,k1)*dth1*dph2
     +        + ut3(i1,j1,k2)*dth2*dph1+ut3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(ut3(i2,j1,k1)*dth2*dph2+ut3(i2,j2,k1)*dth1*dph2
     +        + ut3(i2,j1,k2)*dth2*dph1+ut3(i2,j2,k2)*dth1*dph1)*dra1
       var8(5)=(up3(i1,j1,k1)*dth2*dph2+up3(i1,j2,k1)*dth1*dph2
     +        + up3(i1,j1,k2)*dth2*dph1+up3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(up3(i2,j1,k1)*dth2*dph2+up3(i2,j2,k1)*dth1*dph2
     +        + up3(i2,j1,k2)*dth2*dph1+up3(i2,j2,k2)*dth1*dph1)*dra1
       var8(6)=(br3(i1,j1,k1)*dth2*dph2+br3(i1,j2,k1)*dth1*dph2
     +        + br3(i1,j1,k2)*dth2*dph1+br3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(br3(i2,j1,k1)*dth2*dph2+br3(i2,j2,k1)*dth1*dph2
     +        + br3(i2,j1,k2)*dth2*dph1+br3(i2,j2,k2)*dth1*dph1)*dra1
       var8(7)=(bt3(i1,j1,k1)*dth2*dph2+bt3(i1,j2,k1)*dth1*dph2
     +        + bt3(i1,j1,k2)*dth2*dph1+bt3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(bt3(i2,j1,k1)*dth2*dph2+bt3(i2,j2,k1)*dth1*dph2
     +        + bt3(i2,j1,k2)*dth2*dph1+bt3(i2,j2,k2)*dth1*dph1)*dra1
       var8(8)=(bp3(i1,j1,k1)*dth2*dph2+bp3(i1,j2,k1)*dth1*dph2
     +        + bp3(i1,j1,k2)*dth2*dph1+bp3(i1,j2,k2)*dth1*dph1)*dra2
     +        +(bp3(i2,j1,k1)*dth2*dph2+bp3(i2,j2,k1)*dth1*dph2
     +        + bp3(i2,j1,k2)*dth2*dph1+bp3(i2,j2,k2)*dth1*dph1)*dra1

       return
       end




* -----------------------------------------
*
       subroutine rtp2ijk(rtp,ijk)
       implicit none
* interface
       real*8  rtp(3)
       integer ijk(3)
* local
       integer m
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       common  /coord/  rr, theta, phi

       if (rtp(1) .LT. rr(0)) then
         ijk(1) = 0
       else if (rtp(1) .GT. rr(ii)) then
         ijk(1) = ii
       else
         do m = 0, ii
           if ((rtp(1).GE.rr(m)).AND.(rtp(1).LE.rr(m+1))) ijk(1) = m
         enddo
       endif

       if (rtp(2) .LT. 0.0D+00) then
         ijk(2) = 0
       else if (rtp(2) .GT. theta(jj-1)) then
         ijk(2) = jj - 1
       else
         do m = 0, jj - 1
           if ((rtp(2).GE.theta(m  )).AND.
     +         (rtp(2).LE.theta(m+1))) ijk(2) = m
         enddo
       endif

       if (rtp(3) .LE. 0.0D+00)      rtp(3) = rtp(3) + 2.0D+00 * pi
       if (rtp(3) .GT. 2.0D+00 * pi) rtp(3) = rtp(3) - 2.0D+00 * pi
       ijk(3) = 0
       do m = 0, kk
         if ((rtp(3).GE.phi(m)).AND.
     +       (rtp(3).LE.phi(m+1))) ijk(3) = m
       enddo

       return
       end

*
* --------------------------------------------------------------------
*
       subroutine readinit(gamma,omega,v0,b0,n0,tmp0,r0,rr,theta,phi)
       implicit none
*
       real*8  gamma,omega,v0,b0,n0,tmp0,r0
*
       real*8  a1, a2, t0
       integer idummy
       character*50 strdummy
       character*50 strdumm2
       integer i1, i2, i, j, k
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       character*50  cmnt(30)
       integer imax, jmax, kmax, cmax

       open(unit=1,file='init.dat',status='old')

* read coordinate ... old version does not work for ii > 1000
*       imax = 0
*       jmax = 0
*       kmax = 0
* 200   continue
*         read(1,*,END=299) idummy, a1, a2
*         if (idummy .EQ. 0) goto 299
*         aa = dfloat(idummy) / 1.0D+03 + 1.0D-05
*         iorien = int(aa)
*         igrid = mod(idummy,1000)
*         if (iorien .EQ. 1) then
*           rr(igrid) = a1
*           if (imax .LT. igrid) imax = igrid
*         endif
*         if (iorien .EQ. 2) then
*           theta(igrid) = a1
*           if (jmax .LT. igrid) jmax = igrid
*         endif
*         if (iorien .EQ. 3) then
*           phi(igrid) = a1
*           if (kmax .LT. igrid) kmax = igrid
*         endif
*         goto 200
* 299   continue

* read coordinate
       imax = ii
       jmax = jj
       kmax = kk
       do i = 0, ii + 1
         read(1,*) idummy, a1,a2
         rr(i) = a1
       enddo
       do j = 0, jj
         read(1,*) idummy, a1,a2
         theta(j) = a1
       enddo
       do k = 0, kk + 1
         read(1,*) idummy, a1,a2
         phi(k) = a1
       enddo
       read(1,*) idummy, a1,a2
       write(*,'('' rr   (   0) : '',f10.5)') rr(0)
       write(*,'('' rr   ('',i4,'') : '',f10.5)') imax, rr(imax)
       write(*,'('' theta(   0) : '',f10.5)') theta(0)
       write(*,'('' theta('',i4,'') : '',f10.5)') jmax, theta(jmax)
       write(*,'('' phi  (   0) : '',f10.5)') phi(0)
       write(*,'('' phi  ('',i4,'') : '',f10.5)') kmax, phi(kmax)
* read comment
       cmax = 0
 300   continue
         read(1,'(A)',END=399) strdummy
         cmax = cmax + 1
         cmnt(cmax) = strdummy
         i1 = index(strdummy, 'gamma')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') gamma
           write(*,'('' gamma = '',e11.5)') gamma
         endif
         i1 = index(strdummy, 'omega')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') omega
           write(*,'('' omega = '',e11.5)') omega
         endif
         i1 = index(strdummy, 't0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') t0
           write(*,'('' t0    = '',e11.5)') t0
         endif
         i1 = index(strdummy, 'r0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') r0
           write(*,'('' r0    = '',e11.5)') r0
         endif
         i1 = index(strdummy, 'n0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') n0
           write(*,'('' n0    = '',e11.5)') n0
         endif
         i1 = index(strdummy, 'b0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') b0
           write(*,'('' b0    = '',e11.5)') b0
         endif
         i1 = index(strdummy, 'T_c')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') tmp0
           write(*,'('' T_c   = '',e11.5)') tmp0
         endif

         i1 = index(strdummy, 'v0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') v0
           write(*,'('' v0    = '',e11.5)') v0
         endif
         goto 300
 399   continue
       write(*,'('' Num. of Comment = '',i2)') cmax
       close(1)

       omega = omega * t0 / (3.6D+03 * 2.4D+01 * 1.80D+02) * 3.14

       return
       end
*
* --------------------------------------------------------------------
*
       subroutine readmhd(rt,nnn,lmhd)
       implicit none
*
       integer rt, nnn
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  ro3(0:ii+1,0:jj,0:kk+1),pg3(0:ii+1,0:jj,0:kk+1)
       real*8  ur3(0:ii+1,0:jj,0:kk+1),br3(0:ii+1,0:jj,0:kk+1)
       real*8  ut3(0:ii+1,0:jj,0:kk+1),bt3(0:ii+1,0:jj,0:kk+1)
       real*8  up3(0:ii+1,0:jj,0:kk+1),bp3(0:ii+1,0:jj,0:kk+1)
       logical lmhd
*
       real*8  vv(8)
       integer i, j, k, l
       integer idummy1, jdummy1, kdummy1
       integer idummy2, jdummy2, kdummy2
       logical ldummy
       character*64 flname
       common  /mhdvar/ ro3,pg3,ur3,ut3,up3,br3,bt3,bp3

       write(flname,'(''d'',i6.6,''.'',i4)') nnn + 300000, rt
       inquire(file = flname, exist = ldummy)
       if (ldummy) then
         open(unit=1,file=flname,status='old')
 100     continue
           read(1,'(3i5,5(1x,e11.5),1x,3i5,3(1x,e11.5))',END=99)
     +       idummy1, jdummy1, kdummy1, (vv(l),l=1,5),
     +       idummy2, jdummy2, kdummy2, (vv(l),l=6,8)
           i = idummy1 - 1000
           j = jdummy1 - 2000
           k = kdummy1 - 3000
           ro3(i,j,k) = vv(1)
           pg3(i,j,k) = vv(2)
           ur3(i,j,k) = vv(3)
           ut3(i,j,k) = vv(4)
           up3(i,j,k) = vv(5)
           br3(i,j,k) = vv(6)
           bt3(i,j,k) = vv(7)
           bp3(i,j,k) = vv(8)
           goto 100
 99      continue
         close(1)
         lmhd = .true.
       else
         write(flname,'(''u'',i6.6,''.'',i4)') nnn + 300000, rt
         inquire(file = flname, exist = ldummy)
         if (ldummy) then
           open(unit=1,file=flname,status='old',form='unformatted')
           read(1) ro3
           read(1) pg3
           read(1) ur3
           read(1) ut3
           read(1) up3
           read(1) br3
           read(1) bt3
           read(1) bp3
           close(1)
           lmhd = .true.
         else
           write(*,*) ' File not found ',flname
           lmhd = .false.
         endif
       endif

       return
       end

*
* --------------------------------------------------------------------
*
