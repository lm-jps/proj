*** ---------------------------------------------------
*
* showdaily <== slcalt15.for <=== papslal4.for <=== papslalX.for <====== slcalt14.for
*
* plot at r = r_i
* 8 variables and rot(v) etc.
*
* modified contour tracing
*
*** ---------------------------------------------------
*
       program slcalt
       implicit none
*
       integer numver
       parameter(numver = 15) ! version num.
*
       character*48 flname
       real*8  omega
       real*8  fx1, fy1
       integer rt, ivar, ir, iarr, nnn
       integer rt2(1000), ivar2(1000), ir2(1000)
       integer iarr2(1000), nnn2(1000), icolor2(1000)
       real*8  dcrtrn2(1000)
       integer i, j, k
       integer npage
       integer nframe, nframend, nframe2
       integer frmppg, ifmax, jfmax
       integer icolor, isinlat, irest, ititle, imagneu
       logical ldummy, lmhdfile, ldummy1, ldummy2
       integer idummy
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  xx(0:ix), yy(0:jy)
       real*8  var8(0:ix,0:jy,8)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  ut2(0:ii+1,0:jj,0:kk+1),bt2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)
       real*8  zz(0:ix,0:jy)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  crtrn, dcrtrn
       real*8  r0, t0, v0, b0, d0
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common /vecvec/ vcr, vct, vcp
       common /crt/    crtrn, dcrtrn
       common /crdnt/  xx, yy
       common /var2d/  var8
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0
       common /var3do/ ro2, pg2, ur2, ut2, up2, br2, bt2, bp2
       common /varplt/ zz
       common /sphcrd/ rr, theta, phi
       common /nrmfct/ r0, t0, v0, b0, d0

* Open LOG file
       write(flname,'(''slcalt'',i2.2,''.log'')') numver
       open(unit=13,file=flname,status='unknown')

       write(*,*) 'INPUT Num of Fig.'
       read(*,*) nframend
       write(*,*) 'INPUT Iframe Jframe Fig.perPage'
       read(*,*) ifmax, jfmax, frmppg
       if (ifmax  .LT. 1) ifmax = 1
       if (jfmax  .LT. 1) jfmax = 1
       if (frmppg .GT. ifmax * jfmax) frmppg = ifmax * jfmax

       write(*,*) 'INPUT 4 integers'
       write(*,*) ' Isin(lat) Irest Ititle(for Paper) Imag_neutral'
       write(*,*) '      : 1=Yes / others=No'
       read(*,*)  isinlat, irest, ititle, imagneu

* set coordinate
       call readinit(omega,isinlat)
       if (irest .NE. 1) omega = 0.0D+00

       write(*,*) 'INPUT 6 integers and 1 float'
       write(*,*) '    for Rot NNN Ir Ivar Iarrow Icolor and Dcrtrn'
       write(*,*) ' Ivar :  1-8,39-41  => variables, rhoV'
       write(*,*) ' Ivar : 43          => V_r : special color table'
       write(*,*) ' Ivar : 42,20,21,44 => T, |B|, |B_h|, |V_h|'
       write(*,*) ' Ivar : 45          => beta-ratio^(-1)'
       write(*,*) ' Ivar : 46          => Mag hight/close'
       write(*,*) ' Ivar : 9--11,18,19 => rot(V), |rot(V)|, rot(V)_h'
       write(*,*) ' Ivar : 12--14,26   => rot(B), |rot(B)|'
       write(*,*) ' Ivar : 22--24,25,28=> VxB, |VxB|, |VxB|/|V|/|B|'
       write(*,*) ' Ivar : 15--17,27 => rot(VxB),|rot(VxB)|'
       write(*,*) ' Ivar : 29--31,32 => rot(B)xB,|rot(B)xB|'
       write(*,*) ' Ivar : 33        => rot(B)xB_h'
       write(*,*) ' Ivar : 34-36     => rot(rot(B)xB/rho)'
       write(*,*) ' Ivar : 37        =>|rot(rot(B)xB/rho)|'
       write(*,*) ' Ivar : 38        => rot(rot(B)xB/rho)_h'
       write(*,*) ' Icolor : 1=color 2=grey / others=contour only'
       write(*,*) ' Iarrow : 1, 2, 3 /other = flow, vortex, mag / No'
       write(*,*) ' Dcrtrn in km/s cc^1, or nT '

       do i = 1, nframend
         write(*,'('' # '',i4,'' : '',$)') i
         read(*,*) rt2(i),nnn2(i),ir2(i),ivar2(i),icolor2(i),
     &             iarr2(i),dcrtrn2(i)
         ivar = ivar2(i)
         if (dcrtrn2(i) .GT. 0.0D+00) then
           if  (ivar .EQ. 1) then
             dcrtrn2(i) = dcrtrn2(i) / d0 * 1.67D-24
           else if ((ivar .EQ. 3) .OR. (ivar .EQ. 4) .OR.
     &              (ivar .EQ. 5) .OR. (ivar .EQ. 43).OR.
     &              (ivar .EQ.44)) then
             dcrtrn2(i) = dcrtrn2(i) / v0 * 1.0D+05
           else if (((ivar .GE.  6) .AND. (ivar .LE.  8)) .OR.
     &               (ivar .EQ. 20) .OR.  (ivar .EQ. 21)) then
             dcrtrn2(i) = dcrtrn2(i) / b0 / 1.0D+05
           else if ( (ivar .EQ. 45) .OR. (ivar .EQ. 46) .OR.
     &               (ivar .EQ. 42) ) then
             idummy = 1 ! keep it.
           else
             dcrtrn2(i) = -1.0 ! automatic
           endif
         endif
       enddo

       write(*,'(A)') '--'

       npage = 0
       do nframe = 1, nframend

* Open PS file and Write PS header
         i = mod(nframe, frmppg)
         if (nframend .EQ. 1) i = 1
         if (i .EQ. 1) then
           npage = npage + 1

           ldummy = .true.
           idummy = 0
           do while (ldummy) ! find vacant name
             idummy = idummy + 1
             write(flname,'(''slal'',2i2.2,''.ps.gz'')') numver,idummy
             inquire(file = flname, exist = ldummy2)
             write(flname,'(''slal'',2i2.2,''.ps'')')    numver,idummy
             inquire(file = flname, exist = ldummy1)
             ldummy = (ldummy1 .OR. ldummy2)
           enddo
           open(unit= 2,file=flname,status='unknown')
           write(*,*) 'Now open file : ', flname

           call pshead()
           write(2,'(A)') '/AC { 2 0 360 A C} def'
           write(2,'(A)')
     &       '/AB { M 2 2 rM -4 0 rL 0 -4 rL 4 0 rL C} def'
           write(2,'(A)') ' NP 1 SL'

           fx1 = 270.0
           fy1 = 180.0

           if ((ifmax .GT. 2) .OR. (jfmax .GT. 4))
     &       write(2,'(A)') 'gsave 0 300 translate 0.5 0.5 scale'

         endif

         rt    = rt2(nframe)
         nnn   = nnn2(nframe)
         ir    = ir2(nframe)
         ivar  = ivar2(nframe)
         iarr  = iarr2(nframe)
         dcrtrn = dcrtrn2(nframe)
         icolor = icolor2(nframe)

* restrict the value of i in the case derivation(s) will be calculated.
         if (((ivar .GT.  8) .AND. (ivar .LT. 20)) .OR.
     &       ((ivar .GT. 25) .AND. (ivar .LT. 28)) .OR.
     &       ((ivar .GT. 28) .AND. (ivar .LT. 39)) .OR.
     &        (ivar .EQ. 41)) then
           if (ir .LT.  1) ir = 1
           if (ir .GT. ii) ir = ii
         endif
         if (ivar .EQ. 46) ir = 0

         write(*,'('' Now  Rot NNN Ir Ivar Iflow : '',i4,i7,i4,2i3)')
     &       rt, nnn, ir, ivar, iarr

         write(2,'(A)') '%'
         write(2,'(''% Ivar Ir Iflow : '',4i3)') ivar, ir, iarr
         write(2,'(A)') '%'
         write(2,'(A)') ' gsave'
         i = mod(nframe, frmppg)
         if (i .EQ. 0) i = frmppg
*         j = mod(i+1,ifmax)
*         k = jfmax - (i - 1) / ifmax
         j = (i-1) /jfmax
         k = jfmax - mod(i-1,jfmax)
         write(2,'(2f8.2,'' translate '')')
     &     fx1 * dfloat(j), fy1 * (dfloat(k) - 0.7)

* read data
         if ((nframe .EQ. 1) .OR.
     &       (rt2(nframe) .NE. rt2(nframe-1)) .OR.
     &       (nnn2(nframe) .NE. nnn2(nframe-1)))
     &          call readmhd(rt,nnn,lmhdfile)

         call mhd2alt(ir,lmhdfile)

         nframe2 = nframe
         call onebox(rt,ivar,ir,iarr,nnn,
     &               nframe2,
     &               icolor,ititle,imagneu,omega,i,lmhdfile)

         write(2,'(A)') ' grestore '

* PS Trailer
         i = mod(nframe, frmppg)
         if ((i .EQ. 0) .OR. (nframe .EQ.  nframend)) then
           if ((ifmax .GT. 2) .OR. (jfmax .GT. 4))
     &       write(2,'(A)') 'grestore'
           write(2,'(A)') 'showpage'
           write(2,'(A)') '%%Trailer'
           close(2)
           write(*,*) 'Now close file : ', flname
         endif

       enddo

       close(13)

       stop
       end

**** -------------------------
*
       subroutine onebox(rt,ivar,ir,iarr,nnn,
     &                   nframe,
     &                   icolor,ititle,imagneu,omega,ifrm,lmhdfile)
       implicit none
* interface
       integer rt,ivar,ir,iarr,nnn,icolor,ititle,ifrm,imagneu
       integer nframe
       real*8  omega
       logical lmhdfile
* local
       real*8  fx ,fy, xg0, yg0
       real*8  maxvv, minvv, maxvv2, minvv2
       integer i, j, idummy
       character*26 cbotlab /'abcdefghijklmnopqrstuvwxyz'/
       character*1 cdummy
       real*8  xg, yg, ul !         graphic coordinates & unit length
       real*8  xt1, xt2, xt3, xt4 ! location of title
       real*8  yt1, yt2, yt3, yt4
       real*8  rdummy, pi
       parameter(pi = 3.14159265358979D+00)
       character*9 num2eee
       character*9 sdummy1, sdummy2, sdummy3
       real*8  rdummy1, rdummy2, rdummy3
       integer ncar, loncar, lonleft, lonright
* some for CR instead of date
       integer ncr2
       real*8  rlon, rlat
       real*8  fncr
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  xx(0:ix), yy(0:jy)
       real*8  var8(0:ix,0:jy,8)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  zz(0:ix,0:jy)
       real*8  crtrn, dcrtrn
       real*8  r0, t0, v0, b0, d0
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  ut2(0:ii+1,0:jj,0:kk+1),bt2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common /vecvec/ vcr, vct, vcp
       common /var3do/ ro2, pg2, ur2, ut2, up2, br2, bt2, bp2
       common /crt/    crtrn, dcrtrn
       common /crdnt/  xx, yy
       common /var2d/  var8
       common /varplt/ zz
       common /sphcrd/ rr, theta, phi
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0
       common /nrmfct/ r0, t0, v0, b0, d0

       fx  = 200.0
       fy  = 100.0
       xg0 =  50.0
       yg0 =  50.0

       ul = fx / abs(xx(ix) - xx(0))

* draw
       if (lmhdfile) then
         write(2,'(A)') '% color-or-grey box'
         call cboxdraw(xg0,yg0,fx,fy,ivar,omega,ir,icolor,
     &                  maxvv2,minvv2,nframe)
         write(2,'(A)') '% contour.'
         call cntrdraw(xg0,yg0,fx,fy,ivar,omega,ir,maxvv,minvv,imagneu,
     &                 nframe)
         write(2,'(A)') '% arrow.'
         call arrowdrw(xg0,yg0,fx,fy,omega,ir,iarr)
       endif

       write(2,'(A)') '/Times-Roman findfont 14 scalefont setfont'

* frame
       write(2,'(A)') '% frame etc.'
       write(2,'(A)') ' NP 1 SL 0 SG [] 0 SD'
       write(2,'(2f8.2,'' M '')') xg0, yg0
       write(2,'(2f8.2,'' L '')') xg0 + fx, yg0
       write(2,'(2f8.2,'' L '')') xg0 + fx, yg0 + fy
       write(2,'(2f8.2,'' L C S '')') xg0,  yg0 + fy

* measure right
       if ((icolor .EQ. 1) .OR. (icolor .EQ. 2)) then
         write(*, '('' Var : '',i3,1x,e12.5,'' to '',e12.5)')
     &     ivar, minvv, maxvv
         write(13,'('' Var : '',i3,1x,e12.5,'' to '',e12.5)')
     &     ivar, minvv, maxvv
         write(2,'('' NP '',f9.2,'' SL'')') fy / 20.0 * 1.1 ! 1.1 for overlapping margin

         do j = 0, 19
           rdummy1 = dfloat(j) / 20.0 * (maxvv2 - minvv2) + minvv2
           call defpscol(rdummy1,ivar,icolor,v0,maxvv2,minvv2)
           rdummy1 =(dfloat(j) + 0.5) / 20.0 * fy + yg0
           write(2,'(2f8.2,'' M 15 0 rL S'')') xg0 + fx + 5.0, rdummy1
         enddo
         write(2,'(A)') ' NP 0.2 SL 0 SG [] 0 SD'
         write(2,'(2f8.2,'' M '')')     xg0 + fx +  5.0, yg0
         write(2,'(2f8.2,'' L '')')     xg0 + fx + 20.0, yg0
         write(2,'(2f8.2,'' L '')')     xg0 + fx + 20.0, yg0 + fy
         write(2,'(2f8.2,'' L C S '')') xg0 + fx +  5.0, yg0 + fy
       endif

       write(2,'(A)') ' NP 0.4 SL 0 SG [] 0 SD'
       write(2,'(A)') '/Times-Roman findfont 10 scalefont setfont'
* label : left
       write(2,'(2f8.2,'' M '')')   xg0 - 4.0, yg0
       write(2,'(2f8.2,'' L S '')') xg0 - 4.0, yg0 + fy
       do j = 0, jy
         xg = xg0 - 4.0
         yg = yy(j) * ul + yg0 + fy/2.0D+00
         if ((j .EQ. 0) .OR. (j .EQ. jy) .OR. (j .EQ. jy/2)) then
           write(2,'(2f8.2,'' M -4 0 rL S'')') xg, yg
           write(2,'(2f8.2,'' M '')') xg - 21.0, yg - 4.0
           rdummy = 90.0 - theta(j) / pi * 180.0 + 0.0001
           idummy = int(rdummy)
           if (idummy .LT. 0) then
             write(2,'(''(S'',i2.2,'') show'')') abs(idummy)+1
           else if (idummy .GT. 0) then
             write(2,'(''(N'',i2.2,'') show'')')     idummy
           else
             write(2,'(''( '',i2,  '') show'')')     idummy
           endif
         else
           write(2,'(2f8.2,'' M -2 0 rL S'')') xg, yg
         endif
       enddo
       write(2,'(A)') ' NP'
* label : bottom
       write(2,'(2f8.2,'' M '')')   xg0,      yg0 - 4.0
       write(2,'(2f8.2,'' L S '')') xg0 + fx, yg0 - 4.0
       do i = 0, ix
         xg = xx(i) * ul + xg0
         yg = yg0 - 4.0
         if ((i .EQ. 0) .OR. (i .EQ. ix) .OR. (i .EQ. ix/2)) then
           write(2,'(2f8.2,'' M 0 -4 rL S'')') xg, yg
           write(2,'(2f8.2,'' M '')') xg - 10.0, yg - 16.0
           rdummy = phi(i) / pi * 180.0 + 0.0001
           idummy = int(rdummy)
           write(2,'(''( '',i3.3,  '') show'')')     idummy
         else
           write(2,'(2f8.2,'' M 0 -2 rL S'')') xg, yg
         endif
       enddo
*
       write(2,'(A)') '/Times-Roman findfont 12 scalefont setfont'
       xg = xg0 + fx / 5.0
       yg = yg0 - 29.0
       write(2,'(2f8.2,'' M (Carrington longitude (deg)) show'')')
     &   xg, yg
       write(2,'('' gsave '',2f8.2,'' translate 90 rotate'')')
     &   xg0 - 30.0, yg0 - 10.0
       write(2,'(A)') ' 0 0 M (Heliographic latitude (deg)) show'
       write(2,'(A)') ' grestore'
       write(2,'(A)') ' NP 1 SL 0 SG [] 0 SD'

* title
       xt1 = xg0 + fx - 34.0 ! location of "CR number"
       yt1 = yg0 + fy + 11.0
       xt2 = xg0 !             location of "name of variables"
       yt2 = yg0 + fy + 11.0
       xt3 = xg0 + fx - 60.0 ! location of "heliocentric distance"
       yt3 = yg0 + fy + 22.0
       xt4 = xg0 !             location of "Max Min & step of contour level"
       yt4 = yg0 + fy +  2.0

       cdummy = cbotlab(ifrm:ifrm)
       write(2,'(2f8.2,'' M (('',a1,'')'',$)') xt2, yt2, cdummy

       if (ivar .EQ.  1) write(2,*)'Number density) show'
       if (ivar .EQ.  2) write(2,*)'Pg) show'
       if (ivar .EQ.  3) write(2,*)'V_r) show'
       if (ivar .EQ.  4) write(2,*)'V_theta) show'
       if (ivar .EQ.  5) write(2,*)'V_phi) show'
       if (ivar .EQ.  6) write(2,*)'B_r) show'
       if (ivar .EQ.  7) write(2,*)'B_theta) show'
       if (ivar .EQ.  8) write(2,*)'B_phi) show'
       if (ivar .EQ.  9) write(2,*)'rot(V)_r) show'
       if (ivar .EQ. 10) write(2,*)'rot(V)_theta) show'
       if (ivar .EQ. 11) write(2,*)'rot(V)_phi) show'
       if (ivar .EQ. 12) write(2,*)'rot(B)_r) show'
       if (ivar .EQ. 13) write(2,*)'rot(B)_theta) show'
       if (ivar .EQ. 14) write(2,*)'rot(B)_phi) show'
       if (ivar .EQ. 15) write(2,*)'rot(VxB)_r) show'
       if (ivar .EQ. 16) write(2,*)'rot(VxB)_theta) show'
       if (ivar .EQ. 17) write(2,*)'rot(VxB)_phi) show'
       if (ivar .EQ. 18) write(2,*)'|rot(V)|) show'
       if (ivar .EQ. 19) write(2,*)'|rot(V)_h|) show'
       if (ivar .EQ. 20) write(2,*)'|B|) show'
       if (ivar .EQ. 21) write(2,*)'|B_h|) show'
       if (ivar .EQ. 22) write(2,*)'VxB_r) show'
       if (ivar .EQ. 23) write(2,*)'VxB_theta) show'
       if (ivar .EQ. 24) write(2,*)'VxB_phi) show'
       if (ivar .EQ. 25) write(2,*)'|VxB|) show'
       if (ivar .EQ. 26) write(2,*)'|rot(B)|) show'
       if (ivar .EQ. 27) write(2,*)'|rot(VxB)|) show'
       if (ivar .EQ. 28) write(2,*)'|VxB|/|V|/|B|) show'
       if (ivar .EQ. 29) write(2,*)'rot(B)xB_r) show'
       if (ivar .EQ. 30) write(2,*)'rot(B)xB_theta) show'
       if (ivar .EQ. 31) write(2,*)'rot(B)xB_phi) show'
       if (ivar .EQ. 32) write(2,*)'|rot(B)xB|) show'
       if (ivar .EQ. 33) write(2,*)'rot(B)xB_h) show'
       if (ivar .EQ. 34) write(2,*)'rot(rot(B)xB/rho)_r) show'
       if (ivar .EQ. 35) write(2,*)'rot(rot(B)xB/rho)_theta) show'
       if (ivar .EQ. 36) write(2,*)'rot(rot(B)xB/rho)_phi) show'
       if (ivar .EQ. 37) write(2,*)'|rot(rot(B)xB/rho)|) show'
       if (ivar .EQ. 38) write(2,*)'rot(rot(B)xB/rho)_h) show'
       if (ivar .EQ. 39) write(2,*)'rho V_r) show'
       if (ivar .EQ. 40) write(2,*)'rho V_theta) show'
       if (ivar .EQ. 41) write(2,*)'rho V_phi) show'
       if (ivar .EQ. 42) write(2,*)'Temperature) show'
       if (ivar .EQ. 43) write(2,*)'V_r) show'
       if (ivar .EQ. 44) write(2,*)'|V_h|) show'
       if (ivar .EQ. 45) write(2,*)'1/Beta) show'
       if (ivar .EQ. 46) write(2,*)'Mag Hight) show'


       if (ititle .NE. 1) then
         write(2,'(2f8.2,'' M '',$)') xt1 - 50.0, yt1

         open(12,file='carr.txt',status='old')
         read(12,*) ncr2
         read(12,*) rlon ! rlat, rlon
         close(12)
         fncr = float(ncr2) + (360.0 - rlon) / 360.0
         write(sdummy1,'(''CR'',f7.2)') fncr
** substitute calling date() at Stan.
*         call date(sdummy1)
* substitute calling date()
*         open(12,file='date.txt',status='old')
*         read(12,*) sdummy1
*         close(12)
*
         write(2,'(''('',a9,'') show '')') sdummy1
*         write(2,'(''(CR'',i4.4,'':'',i7,'') show '')') rt, nnn
         write(2,'(2f8.2,'' M '',$)') xt3, yt3
         write(2,'(''(R_'',i2,''='',f7.4,'') show '')') ir, rr(ir)
         write(2,'(A)') '/Times-Roman findfont 10 scalefont setfont'
         rdummy1 = maxvv
         rdummy2 = minvv
         rdummy3 = dcrtrn
         sdummy1 = num2eee(rdummy1)
         sdummy2 = num2eee(rdummy2)
         sdummy3 = num2eee(rdummy3)
         write(2,'(2f8.2,'' M '',$)') xt4, yt4
         write(2,'(''(Max,Min,D = '',2(a9,'',''),a9,'') show'')')
     &        sdummy1, sdummy2, sdummy3
       else
*
         open(12,file='carr.txt',status='old')
         read(12,*) ncr2
         read(12,*) rlon ! rlat, rlon
         close(12)
         fncr = float(ncr2) + (360.0 - rlon) / 360.0
         write(sdummy1,'(''CR'',f7.2)') fncr
**
*         call date(sdummy1)
* substitute calling date()
*         open(12,file='date.txt',status='old')
*         read(12,*) sdummy1
*         close(12)
**
*         write(2,'(2f8.2,'' M '',$)') xt3 - 42.0, yt3
*         write(2,'(''(1 day before '',a9,'') show '')') sdummy1
*         write(2,'(''(Simulated on '',a9,'') show '')') sdummy1
         write(2,'(2f8.2,'' M '',$)') xt3+10.0, yt3
         write(2,'(''('',a9,'') show '')') sdummy1

*         write(2,'(2f8.2,'' M '',$)') xt3 + 10.0, yt3
*         write(2,'(''('',a9,'') show '')') sdummy1
*         write(2,'(''(CR'',i4.4,'':'',i7,'') show '')') rt, nnn

         write(2,'(2f8.2,'' M '',$)') xt1 - 50.0, yt1
         write(2,'(''(R='',f8.3,'' R_sun) show '')')
     &                             rr(ir) * r0 / 6.96E+10
         write(2,'(''%'',2f8.2,'' M '',$)') xt1 - 50.0, yt1
         write(2,'(''(%R='',f8.3,'' AU) show '')')
     &                             rr(ir) * r0 / 1.496E+13
         write(2,'(A)') '/Times-Roman findfont 10 scalefont setfont'
         write(2,'(2f8.2,'' M '',$)') xt4, yt4
         if  (ivar .EQ. 1) then
*           rdummy1 = maxvv  * d0 / 1.67E-24 !<==  mass of hydrogen
*           rdummy2 = minvv  * d0 / 1.67E-24
*           rdummy3 = dcrtrn * d0 / 1.67E-24
*           sdummy1 = num2eee(rdummy1)
*           sdummy2 = num2eee(rdummy2)
*           sdummy3 = num2eee(rdummy3)
*           write(2,'(''(Max,Min,D = '',2(a9,'',''),a9,$)')
*     &        sdummy1, sdummy2, sdummy3
*           write(2,'(A)') ' (/cc)) show '

           rdummy1 = minvv  * d0 / 1.67E-24 !<==  mass of hydrogen
           rdummy2 =(minvv +dcrtrn) * d0 / 1.67E-24
           rdummy3 = maxvv * d0 / 1.67E-24
           sdummy1 = num2eee(rdummy1)
           sdummy2 = num2eee(rdummy2)
           sdummy3 = num2eee(rdummy3)
           write(2,'(''(Contour= '',2(a9,'',''),'' , ,'',a9,$)')
     &        sdummy1, sdummy2, sdummy3
           write(2,'(A)') ' /cc) show '
         else if ((ivar .EQ. 3) .OR. (ivar .EQ. 4) .OR.
     &            (ivar .EQ. 5) .OR. (ivar .EQ. 43).OR.
     &            (ivar .EQ.44)) then
*           rdummy1 = maxvv  * v0 / 1.00E+05 ! km/sec
*           rdummy2 = minvv  * v0 / 1.00E+05
*           rdummy3 = dcrtrn * v0 / 1.00E+05
*           sdummy1 = num2eee(rdummy1)
*           sdummy2 = num2eee(rdummy2)
*           sdummy3 = num2eee(rdummy3)
*           write(2,'(''(Max,Min,D = '',2(a9,'',''),a9,$)')
*     &        sdummy1, sdummy2, sdummy3
*           write(2,'(A)') ' (km/s)) show '

           rdummy1 = minvv  * v0 / 1.00E+05 ! km/sec
           rdummy2 =(minvv +dcrtrn)  * v0 / 1.00E+05
           rdummy3 = maxvv * v0 / 1.00E+05
           sdummy1 = num2eee(rdummy1)
           sdummy2 = num2eee(rdummy2)
           sdummy3 = num2eee(rdummy3)
           write(2,'(''(Contour= '',2(a9,'',''),'' , ,'',a9,$)')
     &        sdummy1, sdummy2, sdummy3
           write(2,'(A)') ' km/s) show '
         else if (((ivar .GE.  6) .AND. (ivar .LE.  8)) .OR.
     &             (ivar .EQ. 20) .OR.  (ivar .EQ. 21)) then
*           rdummy1 = maxvv  * b0 * 1.00E+05 ! nT
*           rdummy2 = minvv  * b0 * 1.00E+05
*           rdummy3 = dcrtrn * b0 * 1.00E+05
*           sdummy1 = num2eee(rdummy1)
*           sdummy2 = num2eee(rdummy2)
*           sdummy3 = num2eee(rdummy3)
*           write(2,'(''(Max,Min,D = '',2(a9,'',''),a9,$)')
*     &        sdummy1, sdummy2, sdummy3
*           write(2,'(A)') ' (nT)) show '

           rdummy1 = minvv * b0 * 1.00E+05 ! nT
           rdummy2 =(minvv + dcrtrn)  * b0 * 1.00E+05
           rdummy3 = maxvv * b0 * 1.00E+05
           sdummy1 = num2eee(rdummy1)
           sdummy2 = num2eee(rdummy2)
           sdummy3 = num2eee(rdummy3)
           write(2,'(''(Contour= '',2(a9,'',''),'' , ,+'',a9,$)')
     &        sdummy1, sdummy2, sdummy3
           write(2,'(A)') ' nT) show '
         else if (ivar .GE. 42) then
           rdummy1 = maxvv  ! / 1.00E+06
           rdummy2 = minvv  ! / 1.00E+06
           rdummy3 = dcrtrn ! / 1.00E+06
           sdummy1 = num2eee(rdummy1)
           sdummy2 = num2eee(rdummy2)
           sdummy3 = num2eee(rdummy3)
           write(2,'(''(Max,Min,D = '',2(a9,'',''),a9,$)')
     &        sdummy1, sdummy2, sdummy3
           write(2,'(A)') ' (MK)) show '
         else
           rdummy1 = maxvv
           rdummy2 = minvv
           rdummy3 = dcrtrn
           sdummy1 = num2eee(rdummy1)
           sdummy2 = num2eee(rdummy2)
           sdummy3 = num2eee(rdummy3)
           write(2,'(''(Max,Min,D = '',2(a9,'',''),a9,'') show'')')
     &        sdummy1, sdummy2, sdummy3
         endif

         if (.false.) then
* draw updated part
*         if ((nframe .EQ. 1) .OR. (nframe .EQ. 3)) then
           open(12,file='carr.txt',status='old')
           read(12,*) ncar
           read(12,*) loncar
           close(12)

           if (mod(nframe,3) .EQ. 1) then ! Vr
             write(2,'(A)') ' 1 1.0 1.0 SR 2 SL [4 2] 0 setdash' ! line
           else if (mod(nframe,3) .EQ. 2) then ! den.
             write(2,'(A)') ' 1 0.5 0.5 SR 2 SL [4 2] 0 setdash' ! line
           else ! mag
             write(2,'(A)') ' 1 1.0 0.0 SR 2 SL [4 2] 0 setdash' ! line
           endif
           if (loncar .LT. 60) then
             lonleft = 0
             lonright= loncar + 60
             write(2,'(2f8.2,'' M '',2f8.2,'' L '',$)')
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.05,
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.05
             write(2,'(2f8.2,'' L '',2f8.2,'' L '',$)')
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.95,
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.95
             write(2,'(A)') ' S NP'

             lonleft = loncar - 60 + 360
             lonright= 360
             write(2,'(2f8.2,'' M '',2f8.2,'' L '',$)')
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.05,
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.05
             write(2,'(2f8.2,'' L '',2f8.2,'' L '',$)')
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.95,
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.95
             write(2,'(A)') ' S NP'
           else if (loncar .GT. 300) then
             lonleft = loncar - 60
             lonright= 360
             write(2,'(2f8.2,'' M '',2f8.2,'' L '',$)')
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.05,
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.05
             write(2,'(2f8.2,'' L '',2f8.2,'' L '',$)')
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.95,
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.95
             write(2,'(A)') ' S NP'

             lonleft = 0
             lonright= loncar + 60 - 360
             write(2,'(2f8.2,'' M '',2f8.2,'' L '',$)')
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.05,
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.05
             write(2,'(2f8.2,'' L '',2f8.2,'' L '',$)')
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.95,
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.95
             write(2,'(A)') ' S NP'
           else
             lonleft = loncar - 60
             lonright= loncar + 60
             write(2,'(2f8.2,'' M '',2f8.2,'' L '',$)')
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.05,
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.05
             write(2,'(2f8.2,'' L '',2f8.2,'' L '',$)')
     &         lonright / 360.0 * fx + xg0, yg0 + fy*0.95,
     &         lonleft  / 360.0 * fx + xg0, yg0 + fy*0.95
             write(2,'(A)') ' C S NP'
           endif

           write(2,'(A)') ' 0 SG 1 SL [] 0 SD ' ! reset line property
*         endif
         endif ! endif false
       
       endif

       return
       end


*
** ---------------------
*
       character*9 function num2eee(aaa)
       implicit none
       real*8  aaa
*
       integer shisuu
       real*8  rdummy, seibun
       character*9 sdummy

       rdummy = abs(aaa)
       if (rdummy .GT. 1.0D-40) then
         rdummy = dlog10(rdummy) + 0.000001
         shisuu = int(rdummy)
         if (shisuu .LT. 0) shisuu = shisuu - 1
         seibun = aaa / 10.0**shisuu
       else
         shisuu = 0
         seibun = 0.0
       endif

       if (aaa .GE. 1.0D-40) then
         if (shisuu .GE. 0) then
           write(sdummy,'('' '',f4.2,''E+'',i2.2)')
     &             seibun,      shisuu
         else
           write(sdummy,'('' '',f4.2,''E-'',i2.2)')
     &             seibun,  abs(shisuu)
         endif
       else if (aaa .LT. -1.0D-40) then
         aaa = dabs(aaa)
         if (shisuu .GE. 0) then
           write(sdummy,'(''-'',f4.2,''E+'',i2.2)')
     &        dabs(seibun),     shisuu
         else
           write(sdummy,'(''-'',f4.2,''E-'',i2.2)')
     &        dabs(seibun), abs(shisuu)
         endif
       else
         sdummy = ' 0.00E+00'
       endif

       num2eee = sdummy

       return
       end



* --------------------------------------------------------------
*
       subroutine arrowdrw(xg0,yg0,fx,fy,omega,ir,iarr)
       implicit none
* interface
       real*8  xg0,yg0,fx,fy
       real*8  omega
       integer ir, iarr
* local
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  ut(0:ix,0:jy), up(0:ix,0:jy)
       integer ixs, ixe, jys, jye
       parameter(ixs = 0, ixe = ix, jys = 1, jye = jy - 1)
       integer i, j, k
       real*8  xg, yg, vx, vy, vv, la, ul
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       real*8  vvave, aa, ra
       real*8  valatjk
       integer ivar2
*
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  xx(0:ix), yy(0:jy)
       real*8  var8(0:ix,0:jy,8)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       common /sphcrd/ rr, theta, phi
       common /crdnt/  xx, yy
       common /var2d/  var8
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0

       ra = rr(ir)

       ul = fx / abs(xx(ix) - xx(0))
       la = 3.0D+00 ! arrow length
       i = 0

* set data to be plotted with contour map
       if (iarr .EQ. 1) then
         vvave = 0.0D+00
         do j = 0, jy
           do k = 0, ix
             ivar2 = 4
             ut(k,j) = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 5
             up(k,j) = valatjk(i,j,k,omega,ivar2,ir)
             aa = ut(k,j)**2 + up(k,j)**2
             aa = dsqrt(aa)
             vvave = vvave + aa
           enddo
         enddo
         vvave = vvave / dfloat(ix+1) / dfloat(jy+1)
       else if (iarr .EQ. 3) then
         vvave = 0.0D+00
         do j = 0, jy
           do k = 0, ix
             ivar2 = 7
             ut(k,j) = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 8
             up(k,j) = valatjk(i,j,k,omega,ivar2,ir)
             aa = ut(k,j)**2 + up(k,j)**2
             aa = dsqrt(aa)
             vvave = vvave + aa
           enddo
         enddo
         vvave = vvave / dfloat(ix+1) / dfloat(jy+1)
       else if (iarr .EQ. 2) then
         vvave = 0.0D+00
         do j = 1, jy - 1
           do k = 1, ix
             ivar2 = 10
             ut(k,j) = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 11
             up(k,j) = valatjk(i,j,k,omega,ivar2,ir)
             aa = ut(k,j)**2 + up(k,j)**2
             aa = dsqrt(aa)
             vvave = vvave + aa
           enddo
           ut(0,j) = ut(ix,j)
           up(0,j) = up(ix,j)
         enddo
         do k = 0, ix
           ut(k, 0) = 0.0D+00
           up(k, 0) = 0.0D+00
           ut(k,jy) = 0.0D+00
           up(k,jy) = 0.0D+00
         enddo
         vvave = vvave / dfloat(ix) / dfloat(jy-1)
       else
         goto 300
       endif

* draw arrow
       write(2,'(A)') '0.2 SL [] 0 SD 0 SG'
       do i = ixs, ixe
         do j = jys, jye

           xg = xx(i) * ul + xg0
           yg = yy(j) * ul + yg0 + fy/2.0D+00
           write(2,'(2f8.2,'' 1 CR F'')') xg, yg
           vy = -ut(i,j)
           vx =  up(i,j)
           vv = vx**2 + vy**2
           if ((vvave .LT. 1.0D-10) .OR. (vv .LT. 1.0D-10)) then
             vx = 0.0D+00
             vy = 0.0D+00
           else
             vv = dsqrt(vv)
             vx = vx * la / vv ! / vvave !  / vv**0.5D+00
             vy = vy * la / vv ! / vvave !  / vv**0.5D+00
           endif
           write(2,'(2f8.2,'' M '',2f8.2,'' rL S'')')
     &         xg, yg, vx, vy

         enddo
       enddo

 300   continue

       return
       end


*
* ---------------------------------------------------------------
*
       subroutine cboxdraw(xg0,yg0,fx,fy,ivar,omega,ir,icolor,
     &                     maxvv,minvv,nframe)
       implicit none
* interface
       integer ivar, ir, icolor
       integer nframe
       real*8  xg0, yg0, fx, fy, omega
* local
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       integer ixs, ixe, jys, jye
       parameter(ixs = 0, ixe = ix, jys = 1, jye = jy - 1) ! range of contour map
       real*8  aa, bb, cc
       real*8  xb1, xb2, yb1, yb2
       integer i2, j2
       integer ivar2
       real*8  ul
       real*8  maxvv, minvv, maxabs
       real*8  r0, t0, v0, b0, d0
* common
       real*8  zz(0:ix,0:jy)
       real*8  xx(0:ix), yy(0:jy)
       real*8  var8(0:ix,0:jy,8)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  crtrn, dcrtrn
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  ut2(0:ii+1,0:jj,0:kk+1),bt2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common /vecvec/ vcr, vct, vcp
       common /var3do/ ro2, pg2, ur2, ut2, up2, br2, bt2, bp2
       common /varplt/ zz
       common /crdnt/  xx, yy
       common /var2d/  var8
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0
       common /sphcrd/ rr, theta, phi
       common /crt/    crtrn, dcrtrn
       common /nrmfct/ r0, t0, v0, b0, d0

       ul = fx / abs(xx(ix) - xx(0))

       ivar2 = ivar
       if (ivar .EQ. 43) ivar2 = 3

* set data to be plotted with contour map
       call setzz(ivar2,omega,ir)

* get statistic properties to determine the contour level
* often not-used.
       maxvv = -1.0E+20
       minvv =  1.0E+20
       do i2 = ixs, ixe
         do j2 = jys, jye
           if (zz(i2,j2) .GT. maxvv) maxvv = zz(i2,j2)
           if (zz(i2,j2) .LT. minvv) minvv = zz(i2,j2)
         enddo
       enddo

* fix maxvv and min for dailyMHD
       if (nframe .EQ. 1) then !             V at i = 0
         maxvv = 15.0E+00 * 1.0E+05 / v0
         minvv = -1.0E+00 * 1.0E+05 / v0
       else if (nframe .EQ. 2) then !        N
         maxvv = 2.2E+08 / (d0 / 1.67E-24)
         minvv = 2.0E+08 / (d0 / 1.67E-24)
       else if (nframe .EQ. 3) then !        B
         maxvv =  5.0E+00 / b0
         minvv = -5.0E+00 / b0
       else if (nframe .EQ. 4) then !        i = 35
         maxvv = 80.0E+00 * 1.0E+05 / v0
         minvv = -1.0E+00 * 1.0E+05 / v0
       else if (nframe .EQ. 5) then
         maxvv = 1.40E+06 / (d0 / 1.67E-24)
         minvv = 0.95E+06 / (d0 / 1.67E-24)
       else if (nframe .EQ. 6) then
         maxvv =  1.5E-01 / b0
         minvv = -1.5E-01 / b0
       else if (nframe .EQ. 7) then !         i = 48
         maxvv = 220.0E+00 * 1.0E+05 / v0
         minvv = 190.0E+00 * 1.0E+05 / v0
       else if (nframe .EQ. 8) then
         maxvv = 1.5E+04 / (d0 / 1.67E-24)
         minvv = 1.0E+04 / (d0 / 1.67E-24)
       else if (nframe .EQ. 9) then
         maxvv =  4.0E-03 / b0
         minvv = -4.0E-03 / b0
       else
         write(*,*) 'unexpected num of frame'
         stop
       endif

       aa = dabs(maxvv)
       bb = dabs(minvv)
       maxabs = dmax1(aa,bb)

       cc = maxvv - minvv
       cc = dabs(cc)
       if (ivar .EQ. 43) cc = cc + 1.0 ! make criteron invalid

* draw color box
       if (((icolor .EQ. 1) .OR. (icolor .EQ. 2)) .AND.
     &     (cc .GT. 1.0D-30)) then
         if (icolor .EQ. 1) write(2,'(A)') '% color box'
         if (icolor .EQ. 2) write(2,'(A)') '% grey  box'
         do i2 = ixs, ixe
           do j2 = jys, jye
* determine location of color box
             if (i2 .EQ. 0) then
               xb1 = xg0 + ul * xx(i2)
             else
               xb1 = xg0 + ul *(xx(i2)+xx(i2-1)) / 2.0D+00
               xb1 = xb1 - 0.6
             endif
             if (i2 .EQ. ix) then
               xb2 = xg0 + ul * xx(i2)
             else
               xb2 = xg0 + ul *(xx(i2)+xx(i2+1)) / 2.0D+00
               xb2 = xb2 + 0.6
             endif
             if (j2 .EQ. 0) then
               yb1 = yg0 + fy/2.0D+00 + ul * yy(j2)
             else
               yb1 = yg0 + fy/2.0D+00 + ul *(yy(j2)+yy(j2-1))/2.0D+00
               yb1 = yb1 + 0.6
             endif
             if (j2 .EQ. jy) then
               yb2 = yg0 + fy/2.0D+00 + ul * yy(j2)
             else
               yb2 = yg0 + fy/2.0D+00 + ul *(yy(j2)+yy(j2+1))/2.0D+00
               yb2 = yb2 - 0.6
             endif
* determine color
             cc = zz(i2,j2)
             call defpscol(cc,ivar,icolor,v0,maxvv,minvv)
* write color box
             write(2,'(2f8.2,'' M '',2f8.2,'' L'',$)')
     &         xb1, yb1, xb2, yb1
             write(2,'(2f8.2,'' L '',2f8.2,'' L C F NP'')')
     &         xb2, yb2, xb1, yb2
           enddo
         enddo
       endif
       write(2,'(A)') ' 0 SG' ! reset color : black

       return
       end


*
* ---------------------------------------------------------------
*
       subroutine defpscol(zzval,ivar,icolor,v0,maxvv,minvv)
       implicit none
       real*8  zzval, v0, maxvv, minvv
       integer ivar,icolor
*
       real*8  aa, bb, cc, red, gre, blu
       real*8  maxabs

       aa = dabs(maxvv)
       bb = dabs(minvv)
       maxabs = dmax1(aa,bb)

       if (ivar .EQ. 43) then
         cc = zzval * v0 / 1.0D+05
         if (icolor .EQ. 1) then
           if (cc .LT. 250.0) then
             red = 0.5
             gre = 0.0
             blu = 0.0
           else if (cc .LT. 300.0) then
             red = (cc - 200.0) / 100.0
             gre = 0.0
             blu = 0.0
           else if (cc .LT. 400.0) then
             red = 1.0
             gre = (cc - 300.0) / 100.0
             blu = 0.0
           else if (cc .LT. 500.0) then
             red = (500.0 - cc) / 100.0
             gre = 1.0
             blu = 0.0
           else if (cc .LT. 600.0) then
             red = 0.0
             gre = 1.0
             blu = (cc - 500.0) / 100.0
           else if (cc .LT. 700.0) then
             red = 0.0
             gre = (700.0 - cc) / 100.0
             blu = 1.0
           else if (cc .LT. 800.0) then
             red = 0.0
             gre = 0.0
             blu = (900.00 - cc) / 200.0
           else
             red = 0.0
             gre = 0.0
             blu = 0.5
           endif
           if (red .GT. 1.0) red = 1.0
           if (red .LT. 0.0) red = 0.0
           if (gre .GT. 1.0) gre = 1.0
           if (gre .LT. 0.0) gre = 0.0
           if (blu .GT. 1.0) blu = 1.0
           if (blu .LT. 0.0) blu = 0.0
           write(2,'(3f5.2,'' SR'',$)') red, gre, blu
         else
           if (cc .GT. 800.0) cc =  800.0
           if (cc .LT. 250.0) cc =  250.0
           cc = (cc - 250.0) / 550.0
           write(2,'(f4.2,'' SG'',$)') 0.5 + 0.5 * cc
         endif
       else if (ivar .EQ. 3) then ! special for Vr...
         if (minvv .LT. 0.0) then
           cc = zzval / maxabs
         else
           cc =(zzval - minvv) / (maxabs - minvv)
         endif
         if (cc .LT. 0.0) cc = 0.0

           if (cc .LT. 0.0) then
             red = 0.5
             gre = 0.0
             blu = 0.0
           else if (cc .LT. 0.2) then
             red = (cc - 0.2) / 0.4 + 1.0
             gre = 0.0
             blu = 0.0
           else if (cc .LT. 0.4) then
             red = 1.0
             gre = (cc - 0.2) / 0.2
             blu = 0.0
           else if (cc .LT. 0.6) then
             red = (0.6 - cc) / 0.2
             gre = 1.0
             blu = 0.0
           else if (cc .LT. 0.8) then
             red = 0.0
             gre = 1.0
             blu = (cc - 0.6) / 0.2
           else if (cc .LT. 1.0) then
             red = 0.0
             gre = (1.0 - cc) / 0.2
             blu = 1.0
           else 
             red = 0.0
             gre = 0.0
             blu = 1.0
           endif
           if (red .GT. 1.0) red = 1.0
           if (red .LT. 0.0) red = 0.0
           if (gre .GT. 1.0) gre = 1.0
           if (gre .LT. 0.0) gre = 0.0
           if (blu .GT. 1.0) blu = 1.0
           if (blu .LT. 0.0) blu = 0.0
           write(2,'(3f5.2,'' SR'',$)') red, gre, blu
       else
         if (maxvv * minvv .LT. 0.0D+00) then
           aa = zzval / maxabs
           aa = dabs(aa)
           if (aa .GT. 1.0) aa = 1.0
           if (icolor .EQ. 1) then
             if (zzval .GE. 0.0D+00) then
               write(2,'(2f5.2,'' 1 SR'',$)')
     &           (1.0D+00 - aa), (1.0D+00 - aa) ! blue for positive
             else
               write(2,'('' 1 '',2f5.2,'' SR'',$)')
     &            (1.0D+00 - aa), (1.0D+00 - aa) ! red for negative
             endif
           else
             write(2,'(f4.2,'' SG'',$)') 1.0 - aa * 0.5
           endif
         else
           cc =  maxvv - minvv
           if (cc .GT. 1.0D-10) then
             cc = (zzval - minvv) / (maxvv - minvv)
             if (cc .GT. 1.0) cc = 1.0
             if (cc .LT. 0.0) cc = 0.0
           else
             cc = 0.0
           endif
           if (icolor .EQ. 1) then
             write(2,'(2f5.2,'' 1 SR'',$)')
     &            (1.0D+00 - cc), (1.0D+00 - cc)
           else
             write(2,'(f4.2,'' SG'',$)') 1.0 - cc * 0.5
           endif
         endif
       endif

       return
       end

*
* ---------------------------------------------------------------
*
       subroutine cntrdraw(xg0,yg0,fx,fy,ivar,omega,ir,
     &                     maxvv0,minvv0,imagneu,nframe)
       implicit none
* interface
       integer ivar, ir, imagneu
       integer nframe
       real*8  xg0, yg0, fx, fy, omega, maxvv0, minvv0
* local
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       integer ixs, ixe, jys, jye
       parameter(ixs = 0, ixe = ix, jys = 1, jye = jy - 1) ! range of contour map
       real*8  cc
       integer i, j, l, iloopmax
       integer ivar2, m2, lc, ls2, le2
       integer iscndrct
       integer crossx(0:ix,0:jy)
       integer crossy(0:ix,0:jy)
       real*8  ul
       real*8  maxvv, minvv
* common
       real*8  zz(0:ix,0:jy)
       real*8  xx(0:ix), yy(0:jy)
       real*8  var8(0:ix,0:jy,8)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  crtrn, dcrtrn
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  ut2(0:ii+1,0:jj,0:kk+1),bt2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common /vecvec/ vcr, vct, vcp
       common /var3do/ ro2, pg2, ur2, ut2, up2, br2, bt2, bp2
       common /varplt/ zz
       common /crdnt/  xx, yy
       common /var2d/  var8
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0
       common /sphcrd/ rr, theta, phi
       common /crt/    crtrn, dcrtrn

       ul = fx / dabs(xx(ix) - xx(0))

       if (imagneu .EQ. 1) then
         iloopmax = 2
       else
         iloopmax = 1
       endif
       do 500 m2 = 1, iloopmax

       if (m2 .EQ. 2) then
         ivar2 = 6
       else
         ivar2 = ivar
         if (ivar .EQ. 43) ivar2 = 3
       endif
* set data to be plotted with contour map
       call setzz(ivar2,omega,ir)

* get statistic properties to determine the contour level
       maxvv = -1.0E+20
       minvv =  1.0E+20
       do i = ixs, ixe
         do j = jys, jye
           if (zz(i,j) .GT. maxvv) maxvv = zz(i,j)
           if (zz(i,j) .LT. minvv) minvv = zz(i,j)
         enddo
       enddo

       write(*, '('' Var : '',i3,1x,e12.5,'' to '',e12.5)')
     &   ivar2, minvv, maxvv
       write(13,'('' Var : '',i3,1x,e12.5,'' to '',e12.5)')
     &   ivar2, minvv, maxvv

* define contour level
       if (m2 .EQ. 2) then
         ls2 = 0
         le2 = 0
       else
         maxvv0  = maxvv
         minvv0  = minvv
         if (dcrtrn .LE. 0.0D+00) then
           dcrtrn = (maxvv - minvv) / 10.0
           ls2 = int(minvv / dcrtrn) - 2
           le2 = int(maxvv / dcrtrn) + 2
         else
*           ls2 = int(minvv / dcrtrn) - 2
*           le2 = int(maxvv / dcrtrn) + 2
           if (minvv * maxvv .LT. 0.0) then
             ls2 =-int(abs(minvv)*0.99999 / dcrtrn)
             le2 = int(maxvv*0.99999 / dcrtrn)
           else if (minvv .GT. 0.0) then
             ls2 = int(minvv / dcrtrn) + 1
             le2 = int(maxvv / dcrtrn)
           else
             ls2 =-int(abs(minvv) / dcrtrn)
             le2 =-int(abs(maxvv) / dcrtrn)
           endif
*           if ((le2 - ls2) .GT. 20) then
*             le2 = ls2 + 20
*             write(*,*) ' Num.of.Contour is too large '
*           endif
         endif
         maxvv0  = dcrtrn * dfloat(le2)
         minvv0  = dcrtrn * dfloat(ls2)
         if (((ivar .EQ. 3) .OR. (ivar .EQ. 43)) .AND.
     &        (minvv0 .LT. 0.0)) minvv0 = 0.0
       endif

       cc = dabs(maxvv - minvv)
       if (cc .LT. 1.0D-20)  goto 666

* Contour
       do 777 l =  ls2, le2

* determine contour level
       if (m2 .EQ. 1) then
         lc = 0
         crtrn = dcrtrn * dfloat(l)
         write(2,'(A,$)') '% contour Ivar Loop.num.  Level = '
         write(2,'(2i5,1x,e11.4)') ivar, l, crtrn
         if (l .LT. lc) write(2,'(A)') ' 0 SG [2 2] 0 SD 0.2 SL' ! broken line
         if (l .EQ. lc) write(2,'(A)') ' 0 SG    [] 0 SD 0.5 SL' ! solid and thick line
         if (l .GT. lc) write(2,'(A)') ' 0 SG    [] 0 SD 0.2 SL' ! solid line

         if ((nframe .EQ. 9) .AND. (l .NE. lc))
     &             write(2,'(A)') ' 0.4 SL' ! thicker ...
*     &             write(2,'(A)') ' 1 SG' ! white line.....

       else
         crtrn = 0.0
         write(2,'(A)') '% contour lines for neutral line'
         write(2,'(A)') ' NP 1 SL 0 SG [4 4] 0 SD'
       endif

* initializing ....
       do 90 i = ixs, ixe
       do 90 j = jys, jye
         crossx(i,j) = 0
         crossy(i,j) = 0
 90    continue

* seek start point : edge
       i = ixs
       do j = jys, jye - 1
         iscndrct = 1
         call tracecnt(i,j,iscndrct,xg0,yg0,ul,fy,
     &                 crossx,crossy)
       enddo

       i = ixe
       do j = jys, jye - 1
         iscndrct = 1
         call tracecnt(i,j,iscndrct,xg0,yg0,ul,fy,
     &                 crossx,crossy)
       enddo

       j = jys
       do i = ixs, ixe - 1
         iscndrct = 2
         call tracecnt(i,j,iscndrct,xg0,yg0,ul,fy,
     &                 crossx,crossy)
       enddo

       j = jye
       do i = ixs, ixe - 1
         iscndrct = 2
         call tracecnt(i,j,iscndrct,xg0,yg0,ul,fy,
     &                 crossx,crossy)
       enddo

* seek start point : domain

       do i = ixs + 1, ixe - 1
         do j = jys + 1, jye - 1
           iscndrct = 1
           call tracecnt(i,j,iscndrct,xg0,yg0,ul,fy,
     &                   crossx,crossy)
         enddo
       enddo

       do i = ixs + 1, ixe - 1
         do j = jys + 1, jye - 1
           iscndrct = 2
           call tracecnt(i,j,iscndrct,xg0,yg0,ul,fy,
     &                   crossx,crossy)
         enddo
       enddo

 777   continue ! contour loop

 666   continue

 500   continue ! m2-loop

       return
       end



*
* ----------------------------------------------------------------
*
       subroutine tracecnt(i,j,iscndrct,xg0,yg0,ul,fy,
     &                     crossx,crossy)
       implicit none
* interface
       integer i, j, iscndrct
       real*8  xg0, yg0, ul, fy
*
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       integer ixs, ixe, jys, jye
       parameter(ixs = 0, ixe = ix, jys = 1, jye = jy - 1) ! range of contour map
       integer crossx(0:ix,0:jy)
       integer crossy(0:ix,0:jy)
* local
       real*8  aa, bb
       logical ltrace
       integer i1, i2, j1, j2
       real*8  x1, x2, y1, y2
       integer istart,jstart,iend,jend
       integer idrctx, idrcty
       integer nline, ifind
       integer ihashis, ihashie
* common
       real*8  zz(0:ix,0:jy)
       real*8  xx(0:ix), yy(0:jy)
       real*8  crtrn, dcrtrn
       common /varplt/ zz
       common /crdnt/  xx, yy
       common /crt/    crtrn, dcrtrn

       ltrace = .false.
       if (iscndrct .EQ. 1) then ! Y direction
         if (crossy(i,j) .EQ. 0)
     &     aa = (zz(i,j) - crtrn) * (zz(i,j+1) - crtrn)
         if (aa .LT. 0.0) ltrace = .true.
         if (zz(i,j) .GT. crtrn) then
           idrctx = -1
           idrcty =  0
           if (i .EQ. ixs) goto 777
         else
           idrctx =  1
           idrcty =  0
           if (i .EQ. ixe) goto 777
         endif
         if (ltrace) then
           bb = (crtrn-zz(i,j)) /  (zz(i,j+1)-zz(i,j))
           x1 = xx(i)
           y1 = yy(j) + (yy(j+1) - yy(j)) * bb
         endif
       else if (iscndrct .EQ. 2) then ! X-direction
         if (crossx(i,j) .EQ. 0)
     &     aa = (zz(i,j) - crtrn) * (zz(i+1,j) - crtrn)
         if (aa .LT. 0.0) ltrace = .true.
         if (zz(i,j) .GT. crtrn) then
           idrctx =  0
           idrcty =  1
           if (j .EQ. jye) goto 777
         else
           idrctx =  0
           idrcty = -1
           if (j .EQ. jys) goto 777
         endif
         if (ltrace) then
           bb = (crtrn-zz(i,j)) / (zz(i+1,j)-zz(i,j))
           x1 = xx(i) + (xx(i+1) - xx(i)) * bb
           y1 = yy(j)
         endif
       else
         write(*,*) ' Wrong value of Iscndrct'
       endif

       nline = 0
       ifind = 1

       if (ltrace) then ! found the candidate for start point

         i1 = i
         j1 = j

         istart = i1
         jstart = j1

         if ((i1 .EQ. ixs) .OR. (i1 .EQ. ixe) .OR.
     &       (j1 .EQ. jys) .OR. (j1 .EQ. jye)) then
           ihashis = 1
         else
           ihashis = 0
         endif

         write(2,'(2f8.2,'' M %  '',2i4)')
     &     x1*ul+xg0, y1*ul+yg0+fy/2.0,i1,j1

         do while ((nline .LT. 200) .AND. (ifind .EQ. 1))
           call seeknext(i1,j1,i2,j2,x2,y2,idrctx,idrcty,ifind)
           if (ifind .EQ. 1) then
             nline = nline + 1
             if (idrctx .EQ. 0) then
               if (crossx(i2,j2) .EQ. 1) then
                 write(2,'(2f8.2,'' L S %'',2i4)')
     &             x1*ul+xg0,y1*ul+yg0+fy/2.0,i1,j1
                 ifind  = -2
               else
                 write(2,'(2f8.2,'' L'')')
     &             x2*ul+xg0,y2*ul+yg0+fy/2.0
                 crossx(i2,j2) = 1
               endif
             endif
             if (idrcty .EQ. 0) then
               if (crossy(i2,j2) .EQ. 1) then
                 write(2,'(2f8.2,'' L S %'',2i4)')
     &             x1*ul+xg0,y1*ul+yg0+fy/2.0,i1,j1
                 ifind  = -2
               else
                 write(2,'(2f8.2,'' L'')')
     &             x2*ul+xg0,y2*ul+yg0+fy/2.0
                 crossy(i2,j2) = 1
               endif
             endif
             x1 = x2
             y1 = y2
             i1 = i2
             j1 = j2
           endif
           if (ifind .EQ. -1) then
             write(2,'(2f8.2,'' L S %'',4i4)')
     &         x1*ul+xg0,y1*ul+yg0+fy/2.0,istart,jstart,i1,j1
             ifind = -2
           endif
         enddo

         if (ifind .NE. -2) write(*,*) ' Wrong at 111a'
         iend = i1
         jend = j1
         if ((i1 .EQ. ixs) .OR. (i1 .EQ. ixe) .OR.
     &       (j1 .EQ. jys) .OR. (j1 .EQ. jye)) then
           ihashie = 1
         else
           ihashie = 0
         endif
         if (ihashis .NE. ihashie) then
           write(13,'(''Wrong 111b at'',4i4)') istart,jstart,iend,jend
           write(*,'(''Wrong 111b at'',4i4)') istart,jstart,iend,jend
         endif

       endif

 777   continue

       return
       end


*
* ----------------------------------------------------------------
*
       subroutine seeknext(i1,j1,i2,j2,x2,y2,idrctx,idrcty,ifind)
       implicit none
* interface
       integer i1,j1,i2,j2
       real*8  x2,y2
       integer idrctx,idrcty,ifind
* local
       real*8  aa, bb
       integer ia, ib, ja, jb
       integer l
       integer ix2(4), jy2(4)
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       integer ixs, ixe, jys, jye
       parameter(ixs = 0, ixe = ix, jys = 1, jye = jy - 1)
       real*8  zz(0:ix,0:jy)
       real*8  xx(0:ix), yy(0:jy)
       real*8  crtrn, dcrtrn
       common /crdnt/  xx, yy
       common /varplt/ zz
       common /crt/    crtrn, dcrtrn

       ifind = 0
       if (idrcty .EQ. 0) then
         if (idrctx .GT. 0) then
           if (i1 .EQ. ixe) then
             ifind = -1
             goto 450
           endif
           ix2(1) = i1
           jy2(1) = j1 + 1
           ix2(2) = i1 + 1
           jy2(2) = j1 + 1
           ix2(3) = i1 + 1
           jy2(3) = j1
           ix2(4) = i1
           jy2(4) = j1
         endif
         if (idrctx .LT. 0) then
           if (i1 .EQ. ixs) then
             ifind = -1
             goto 450
           endif
           ix2(1) = i1
           jy2(1) = j1
           ix2(2) = i1 - 1
           jy2(2) = j1
           ix2(3) = i1 - 1
           jy2(3) = j1 + 1
           ix2(4) = i1
           jy2(4) = j1 + 1
         endif
       endif
       if (idrctx .EQ. 0) then
         if (idrcty .GT. 0) then
           if (j1 .EQ. jye) then
             ifind = -1
             goto 450
           endif
           ix2(1) = i1
           jy2(1) = j1
           ix2(2) = i1
           jy2(2) = j1 + 1
           ix2(3) = i1 + 1
           jy2(3) = j1 + 1
           ix2(4) = i1 + 1
           jy2(4) = j1
         endif
         if (idrcty .LT. 0) then
           if (j1 .EQ. jys) then
             ifind = -1
             goto 450
           endif
           ix2(1) = i1 + 1
           jy2(1) = j1
           ix2(2) = i1 + 1
           jy2(2) = j1 - 1
           ix2(3) = i1
           jy2(3) = j1 - 1
           ix2(4) = i1
           jy2(4) = j1
         endif
       endif
       do l = 1, 3 ! try three candidate(s)
         ia = ix2(l)
         ja = jy2(l)
         ib = ix2(l+1)
         jb = jy2(l+1)
         aa = (zz(ia,ja) - crtrn) * (zz(ib,jb) - crtrn)
         if (aa .LE. 0.0) then ! found the next point
           ifind = 1
           bb = (crtrn-zz(ia,ja)) / (zz(ib,jb)-zz(ia,ja))
           x2 = xx(ia) + (xx(ib) - xx(ia)) * bb
           y2 = yy(ja) + (yy(jb) - yy(ja)) * bb
           i2 = min(ia,ib)
           j2 = min(ja,jb)
           if (ia .EQ. ib) then
             idrcty = 0
             if (i1 .GE. i2) idrctx = -1
             if (i1 .LT. i2) idrctx =  1
           endif
           if (ja .EQ. jb) then
             idrctx = 0
             if (j1 .GE. j2) idrcty = -1
             if (j1 .LT. j2) idrcty =  1
           endif
           goto 450
         endif
       enddo
       write(*,*) ' can not find the next point !! at ', i1,j1
       write(*,'(''idrctX idrctY '',2i3)') idrctx, idrcty
       write(*,*) (ix2(l),l=1,4)
       write(*,*) (jy2(l),l=1,4)
       write(13,*) ' can not find the next point !! at ', i1,j1
       write(13,'(''idrctX idrctY '',2i3)') idrctx, idrcty
       write(13,*) (ix2(l),l=1,4)
       write(13,*) (jy2(l),l=1,4)

 450   continue

       return
       end



* --------------------------------------------------------------
*
       subroutine readinit(omega,isinlat)
       implicit none
*
       real*8  omega
       integer isinlat
*
       real*8  a1, a2, aa
       integer idummy
       integer iorien, igrid
       character*50 strdummy
       character*50 strdumm2
       integer i1, i2, i, j
       integer imax, jmax, kmax, cmax
       real*8  r0, v0, pi, t0, b0, d0
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  xx(0:ix), yy(0:jy)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       common /crdnt/  xx, yy
       common /sphcrd/ rr, theta, phi
       common /nrmfct/ r0, t0, v0, b0, d0

       open(unit=15,file='init.dat',status='old')
* read coordinate
       imax = 0
       jmax = 0
       kmax = 0
 200   continue
         read(15,*,END=299) idummy, a1, a2
         if (idummy .EQ. 0) goto 299
         aa = float(idummy) / 1.0E+03 + 1.0E-05
         iorien = int(aa)
         igrid = mod(idummy,1000)
         if (iorien .EQ. 1) then
           rr(igrid) = a1
           if (imax .LT. igrid) imax = igrid
         endif
         if (iorien .EQ. 2) then
           theta(igrid) = a1
           if (jmax .LT. igrid) jmax = igrid
         endif
         if (iorien .EQ. 3) then
           phi(igrid) = a1
           if (kmax .LT. igrid) kmax = igrid
         endif
         goto 200
 299   continue

* get unit factor
 300   continue
         read(15,'(A)',END=399) strdummy
         cmax = cmax + 1
*         cmnt(cmax) = strdummy
         i1 = index(strdummy, 'r0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') r0
           write(13,'('' r0     = '',e11.5)') r0
         endif
         i1 = index(strdummy, 't0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') t0
           write(13,'('' t0     = '',e11.5)') t0
         endif
         i1 = index(strdummy, 'v0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') v0
           write(13,'('' v0     = '',e11.5)') v0
         endif
         i1 = index(strdummy, 'b0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') b0
           write(13,'('' b0     = '',e11.5)') b0
         endif
         i1 = index(strdummy, 'rhoc')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') d0
           write(13,'('' d0     = '',e11.5)') d0
         endif
         i1 = index(strdummy, 'omega')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') omega
           write(13,'('' omega0 = '',e11.5)') omega
         endif
         goto 300
 399   continue
       write(13,*) ' omega = ',omega
       write(13,*) ' t0    = ',t0
       write(13,*) ' v0    = ',v0
       write(13,*) ' r0    = ',r0
* not used in this program
       pi = 3.14159265358979
       omega = omega * pi/ 180.0 / (60.0*60.0*24.0) ! [deg/day]=>[rad/sec]
       omega = omega * t0 !                                    =>[rad/unit time]

       write(13,*) ' omega = ',omega

       close(15)

       do i = 0, ix
         xx(i) = phi(i) / pi * 180.0
       enddo

       if (isinlat .EQ. 1) then
         write(13,'(A)') ' YY are determined as cos of theta'
         do j = 0, jy
           aa = theta(j)
           yy(j) = dcos(aa) * 90.0D+00
         enddo
       else
         do j = 0, jy
           write(13,'(A)') ' YY are determined as linear of theta'
           aa = 90.0 - theta(j) / pi * 180.0
           yy(j) = aa
         enddo
       endif

       return
       end



* --------------------------------------------------------------
*
       subroutine  pshead()
       write(2,'(A)') '%!PS-Adobe-2.0'
       write(2,'(A)') '%%BoundingBox: (atend)'
       write(2,'(A)') '%%Creator: hayashi.for'
       write(2,'(A)') '%%page: (atend)'
       write(2,'(A)') '%%EndComments'
       write(2,'(A)') '/M {moveto} def   /L {lineto} def'
       write(2,'(A)') '/rM {rmoveto} def /rL {rlineto} def'
       write(2,'(A)') '/S {stroke} def   /C {closepath} def'
       write(2,'(A)') '/F {fill} def     /SL {setlinewidth} def'
       write(2,'(A)') '/A {arc} def      /SG {setgray} def'
       write(2,'(A)') '/CR { 0 360 A } def '
       write(2,'(A)') '/NP {newpath} def'
       write(2,'(A)') '/SR {setrgbcolor} def'
       write(2,'(A)') '/SD {setdash} def'
       return
       end



* ----------------------------------------------------------------------------
*
       subroutine setzz(ivar,omega,ir)
       implicit none
* interface
       integer ivar, ir
       real*8  omega
*local
       real*8  ra, aa, bb, cc, dd, ee, ff
       integer i, j, k
       integer i1, i2, i3, j1, j2, k1, k2
       real*8  dra, dth, dph
       real*8  valatjk, lrntatjk, maghight
       integer ivar2, ir2
       real*8  a1, a2, a3, alpha
       real*8  rot1, rot2
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  var8(0:ix,0:jy,8)
       real*8  zz(0:ix,0:jy)
       real*8  xx(0:ix), yy(0:jy)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  ut2(0:ii+1,0:jj,0:kk+1),bt2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common /vecvec/ vcr, vct, vcp
       common /var3do/ ro2, pg2, ur2, ut2, up2, br2, bt2, bp2
       common /var2d/  var8
       common /varplt/ zz
       common /crdnt/  xx, yy
       common /sphcrd/ rr, theta, phi
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0

* initialize
       do i = 0, ix
         do j = 0, jy
           zz(i,j) = 0.0D+00
         enddo
       enddo

       i = 0
       ra = rr(ir)

       if (ivar .LE. 8) then
         do j = 0, jy
           do k = 0, ix
             zz(k,j) = valatjk(i,j,k,omega,ivar,ir)
             dd = zz(k,j)
             if ((ivar .EQ. 3) .AND. (ir .EQ. 0) .AND.
     &           (dd .LE. 0.0D+00)) zz(k,j) = -1.0D-01 ! make vr=0 clear at ir = 0
           enddo
         enddo
       else if (ivar .EQ. 42) then
         do j = 0, jy
           do k = 0, ix
             zz(k,j) = valatjk(i,j,k,omega,ivar,ir)
           enddo
         enddo
       else if (ivar .EQ. 44) then ! |V_h|
         do j = 0, jy
           do k = 0, ix
             ivar2 = 4
             dd = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 5
             ee = valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = dsqrt(dd**2 + ee**2)
           enddo
         enddo
* abs. value of magnetic field
       else if (ivar .EQ. 45) then ! B^2 / Pg : note exact beta = B^2/(Pg/gamma)
         do j = 0, jy
           do k = 0, ix
             ivar2 = 6
             dd = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 7
             ee = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 8
             ff = valatjk(i,j,k,omega,ivar2,ir)
             aa = dd**2 + ee**2 + ff**2
             ivar2 = 2
             bb = valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = aa / bb
           enddo
         enddo
* abs. value of magnetic field
       else if ((ivar .GE. 20) .AND. (ivar .LE. 21)) then ! |B|
         do j = 0, jy
           do k = 0, ix
             ivar2 = 6
             dd = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 7
             ee = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 8
             ff = valatjk(i,j,k,omega,ivar2,ir)
             if (ivar .EQ. 20) then
               zz(k,j) = dsqrt(dd**2 + ee**2 + ff**2) !|B|
             else if (ivar .EQ. 21) then
               zz(k,j) = dsqrt(dd**2 + ff**2) !        |B_h|
             endif
           enddo
         enddo
* Lorentz term
       else if ((ivar .GE. 29) .AND. (ivar .LE. 33)) then ! rot(B)xB
         do j = 1, jy - 1
           do k = 1, ix
             ivar2 = 1
              dd = lrntatjk(i,j,k,omega,ivar2,ir) ! rot(B)xB_r
             ivar2 = 2
              ee = lrntatjk(i,j,k,omega,ivar2,ir) ! rot(B)xB_theta
             ivar2 = 3
              ff = lrntatjk(i,j,k,omega,ivar2,ir) ! rot(B)xB_phi
             if (ivar .EQ. 29) then
               zz(k,j) = dd
             else if (ivar .EQ. 30) then
               zz(k,j) = ee
             else if (ivar .EQ. 31) then
               zz(k,j) = ff
             else if (ivar .EQ. 32) then
               zz(k,j) = dsqrt(dd**2 + ee**2 + ff**2)
             else if (ivar .EQ. 33) then
               zz(k,j) = dsqrt(ee**2 + ff**2)
             endif
           enddo
           zz(0,j) = zz(ix,j)
         enddo
* rot(rot(B)xB/rho)
       else if ((ivar .GE. 34) .AND. (ivar .LE. 38)) then ! rot(rot(B)xB/rho)
         alpha =(rr(ir+1) - rr(ir)) / (rr(ir) - rr(ir-1))
         a1 =          1.0D+00 / alpha
         a2 =  alpha - 1.0D+00 / alpha
         a3 = -alpha
         i1 =   1
         i2 =   0
         i3 = - 1
         dra = rr(ir+1) - rr(ir-1)
         dph = phi(3)-phi(1)
         dth = theta(3) - theta(1)
         do j = 2, jy - 2
           j1 = j + 1
           j2 = j - 1
           do k = 1, ix
             k1 = k + 1
             k2 = k - 1
             if (k .EQ. ix) k1 =  1
             if (k .EQ.  1) k2 = ix
             ivar2 = 3
             rot1 = (lrntatjk(i,j1,k,omega,ivar2,ir)/ro0(i,j1,k)
     &              -lrntatjk(i,j2,k,omega,ivar2,ir)/ro0(i,j2,k))
     &            / dth / dsin(theta(j)) / rr(ir)
             ivar2 = 2
             rot2 = (lrntatjk(i,j,k1,omega,ivar2,ir)/ro0(i,j,k1)
     &              -lrntatjk(i,j,k2,omega,ivar2,ir)/ro0(i,j,k2))
     &            / dph / dsin(theta(j)) / rr(ir)
             dd = rot1 - rot2 ! rot(rot(B)xB/rho)_r
             ivar2 = 1
             rot1 = (lrntatjk(i,j,k1,omega,ivar2,ir)/ro0(i,j,k1)
     &              -lrntatjk(i,j,k2,omega,ivar2,ir)/ro0(i,j,k2))
     &            / dph / dsin(theta(j)) / rr(ir)
             ivar2 = 3
             rot2 = (lrntatjk(i1,j,k,omega,ivar2,ir)/ro0(i1,j,k) * a1
     &              +lrntatjk(i2,j,k,omega,ivar2,ir)/ro0(i2,j,k) * a2
     &              +lrntatjk(i3,j,k,omega,ivar2,ir)/ro0(i3,j,k) * a3)
     &            / dra / rr(ir)
             ee = rot1 - rot2 ! rot(rot(B)xB/rho)_theta

             ivar2 = 2
             rot1 = (lrntatjk(i1,j,k,omega,ivar2,ir)/ro0(i1,j,k) * a1
     &              +lrntatjk(i2,j,k,omega,ivar2,ir)/ro0(i2,j,k) * a2
     &              +lrntatjk(i3,j,k,omega,ivar2,ir)/ro0(i3,j,k) * a3)
     &            / dra / rr(ir)
             ivar2 = 1
             rot2 = (lrntatjk(i,j1,k,omega,ivar2,ir) / ro0(i,j1,k)
     &              -lrntatjk(i,j2,k,omega,ivar2,ir) / ro0(i,j2,k))
     &            / dth / rr(ir)
             ff = rot1 - rot2 ! rot(rot(B)xB/rho)_phi
             if (ivar .EQ. 34) then
               zz(k,j) = dd
             else if (ivar .EQ. 35) then
               zz(k,j) = ee
             else if (ivar .EQ. 36) then
               zz(k,j) = ff
             else if (ivar .EQ. 37) then
               zz(k,j) = dsqrt(dd**2 + ee**2 + ff**2)
             else if (ivar .EQ. 38) then
               zz(k,j) = dsqrt(ee**2 + ff**2)
             endif
           enddo
           zz(0,j) = zz(ix,j)
         enddo
* rot flow etc
       else if ((ivar .GE. 9) .AND. (ivar .LE. 17)) then ! rot(V) rot(B) rot(VxB)
         do j = 1, jy - 1
           do k = 1, ix
             zz(k,j) = valatjk(i,j,k,omega,ivar,ir)
           enddo
           zz(0,j) = zz(ix,j)
         enddo
       else if (ivar .EQ. 18) then ! |rot(V)|
         do j = 1, jy - 1
           do k = 1, ix
             ivar2 =  9
             aa =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 10
             bb =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 11
             cc =  valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = dsqrt(aa**2 + bb**2 + cc**2)
           enddo
           zz(0,j) = zz(ix,j)
         enddo
       else if (ivar .EQ. 19) then ! |rot(V_h)|
         do j = 1, jy - 1
           do k = 1, ix
             ivar2 = 10
             bb =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 11
             cc =  valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = dsqrt(bb**2 + cc**2)
           enddo
           zz(0,j) = zz(ix,j)
         enddo
       else if (ivar .EQ. 26) then ! |rot(B)|
         do j = 1, jy - 1
           do k = 1, ix
             ivar2 = 12
             aa =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 13
             bb =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 14
             cc =  valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = dsqrt(aa**2 + bb**2 + cc**2)
           enddo
           zz(0,j) = zz(ix,j)
         enddo
       else if (ivar .EQ. 27) then ! |rot(VxB)|
         do j = 1, jy - 1
           do k = 1, ix
             ivar2 = 15
             aa =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 16
             bb =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 17
             cc =  valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = dsqrt(aa**2 + bb**2 + cc**2)
           enddo
           zz(0,j) = zz(ix,j)
         enddo
       else if (ivar .EQ. 28) then ! |VxB|/|V|/|B|
         do j = 0, jy
           do k = 0, ix
             ivar2 = 22
             aa =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 23
             bb =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 24
             cc =  valatjk(i,j,k,omega,ivar2,ir)
             dd =  dsqrt(aa**2 + bb**2 + cc**2)
             ivar2 =  3
             aa =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 =  4
             bb =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 =  5
             cc =  valatjk(i,j,k,omega,ivar2,ir)
             cc =  cc - rr(ir)*dsin(theta(j))*omega
             ee =  dsqrt(aa**2 + bb**2 + cc**2)
             ivar2 =  6
             aa =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 =  7
             bb =  valatjk(i,j,k,omega,ivar2,ir)
             ivar2 =  8
             cc =  valatjk(i,j,k,omega,ivar2,ir)
             ff =  dsqrt(aa**2 + bb**2 + cc**2)
             if (abs(ee * ff) .GE. 1.0D-30) then
               zz(k,j) = dd / ee / ff
             else
               zz(k,j) = 0.0D+00
             endif
           enddo
         enddo
       else if ((ivar .GE. 22) .AND. (ivar .LE. 24)) then ! VxB
         do j = 0, jy
           do k = 0, ix
             zz(k,j) = valatjk(i,j,k,omega,ivar,ir)
           enddo
         enddo
       else if (ivar .EQ. 25) then ! |VxB|
         do j = 0, jy
           do k = 0, ix
             ivar2 = 22
             aa = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 23
             bb = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = 24
             cc = valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = dsqrt(aa**2 + bb**2 + cc**2)
           enddo
         enddo
       else if ((ivar .GE. 39) .AND. (ivar .LE. 41)) then
         do j = 0, jy
           do k = 0, ix
             ivar2 = 1
             aa = valatjk(i,j,k,omega,ivar2,ir)
             ivar2 = ivar - 36
             bb = valatjk(i,j,k,omega,ivar2,ir)
             zz(k,j) = aa * bb
           enddo
         enddo
       else if (ivar .EQ. 46) then
         do j = 0, jy
           do k = 0, ix
             ivar2 = 1
             ir2   = 0
             aa = maghight(i,j,k,ivar2,ir2)
             zz(k,j) = aa
           enddo
         enddo
       else
         write(*,*)  'Ivar is not available at setzz .. stop'
         write(13,*) 'Ivar is not available at setzz .. stop'
         stop
       endif

       return
       end


*
* ----------------------------------------------------
*
       real*8 function maghight(i,j,k,ivar2,ir2) ! max will be 10.0
       implicit none
       integer i,j,k,ivar2,ir2
* local
       integer idummy, m
       real*8   maxh
       integer linemax
       parameter(linemax = 1000)
       real*8  rtpline(3,0:linemax), avec(0:linemax)
       integer nline
       real*8  rupper2, rbottom2
       logical lfromsun

* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common /vecvec/ vcr, vct, vcp
       common /sphcrd/ rr, theta, phi

       idummy = ivar2 ! supress warnings.
       idummy = i

       idummy = 1
       lfromsun = .true.
       rupper2 = rr(ir2) * 10.0D+00
       if (rupper2 .GT. rr(ii+1)*0.999D+00)
     &        rupper2 = rr(ii+1)*0.999D+00
       rbottom2 = rr(ir2) * 1.01D+00
       rtpline(1,0) = rbottom2 * 1.01D+00
       rtpline(2,0) = theta(j)
       rtpline(3,0) = phi(k)
       call trcline(rtpline,nline,rupper2,rbottom2,
     &              idummy,avec,lfromsun)
       maxh = -1.0D+20
       do m = 0, nline
         if (rtpline(1,m) .GT. maxh) maxh = rtpline(1,m)
       enddo

       if ((rtpline(1,nline  ) .LT. 1.1*rbottom2) .OR.
     &     (rtpline(1,nline-1) .LT. 1.1*rbottom2)) maxh = 0.0
       maghight = maxh

       return
       end



*
** trace line ----------------------------------------------------------
*
       subroutine trcline(rtpline,nline,rupper,rbottom,
     &                    ivar,avec,lfromsun)
       implicit  none
* interface
       integer linemax
       parameter(linemax = 1000)
       real*8  rtpline(3,0:linemax),avec(0:linemax)
       integer nline
       real*8  rupper,rbottom
       integer ivar
       real*8  avecl
       logical lfromsun
* local
       integer idrct
       integer i
       real*8  rtp1(3), rtp2(3)
       real*8  dr, ra, th, ph
       real*8  ra0, th0, ph0
       real*8  pi, dpi
       parameter(pi = 3.14159265358979D+00, dpi = 5.0D-02)
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common  /sphcrd/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp
*
       ra0 = rtpline(1,0)
       th0 = rtpline(2,0)
       ph0 = rtpline(3,0)
       ra = ra0
       th = th0
       ph = ph0
       rtp1(1) = ra
       rtp1(2) = th
       rtp1(3) = ph
* check the direction
       idrct = 1
*       dr = ra * dsin(th) * (2.0D+00 * pi)
*     &    / dfloat(kk) * 0.5D+00 * dfloat(idrct)
       dr = (rr(1) - rr(0)) / 4.0
       call rungevec(rtp1, rtp2, dr, avecl, ivar)
       avec(0) = avecl

       if ((      lfromsun) .AND. (rtp2(1) .LT. rtp1(1))) idrct = -1
       if ((.NOT. lfromsun) .AND. (rtp2(1) .GT. rtp1(1))) idrct = -1

*       write(*,*) idrct

* estimate magnetic vector and get coordinate of line
       nline = 0
       ra = rtpline(1,0)
       th = rtpline(2,0)
       ph = rtpline(3,0)
       rtp1(1) = ra + dfloat(idrct) * 1.0D-03
       rtp1(2) = th
       rtp1(3) = ph
       do while((ra .GE. rbottom) .AND.
     &          (ra .LT. rupper) .AND.
     &          (th .GE. dpi) .AND.
     &          (th .LE. (pi - dpi)) .AND.
     &          (nline .LT. linemax))

         dr = ra * pi / dfloat(jj) / 1.0D+01 * dfloat(idrct)
         nline = nline + 1

         call rungevec(rtp1, rtp2, dr, avecl, ivar) ! 1 --> 2
         ra = rtp2(1)
         th = rtp2(2)
         ph = rtp2(3)
         if (ph .LT. 0.0D+00) then
           ph = ph + 2.0D+00 * pi
         else if (ph .GT. 2.0D+00 * pi) then
           ph = dmod(ph, 2.0D+00 * pi)
         endif
         rtp2(3) = ph

*         if (.NOT. lfromsun) goto 300

         do i = 1, 3
           rtpline(i,nline) = rtp2(i)
           rtp1(i)          = rtp2(i)
         enddo

         avec(nline) = avecl

       enddo        ! End of do-while-loop

 300   continue

*           write(11,'(8x,'' : '',i5,x,6e12.5)')
*     &        nline,  ra0, th0, ph0, ra, th, ph

       return
       end


*
* ---------------------------------------------------------------
*
       subroutine rungevec(rtp1, rtp2, dr, avecl, ivar)  ! rtp1 --> rtp2
       implicit none
* interface
       real*8  rtp1(3), rtp2(3), dr, avecl
       integer ivar
* local
       real*8  bbb, ccc
       real*8  k1(3), k2(3), k3(3), k4(3), k0(3)
       real*8  b1(3), b2(3), b3(3), b4(3), b0(3)
       real*8  rtp3(3)
       real*8  xyz1(3), xyz2(3), xyz3(3)
       integer k
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common  /sphcrd/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp

       call rtp2xyz(rtp1,xyz1)

       rtp3(1) = rtp1(1)
       rtp3(2) = rtp1(2)
       rtp3(3) = rtp1(3)
       call vecatrtp(b1,rtp3)
       ccc = b1(1)
       call vrtp2xyz(b1,rtp3)
       call normvec3(b1,bbb)

       if (ivar .EQ. 1) then
         avecl = bbb
       else if (ivar .EQ. 2) then
         avecl = ccc
       else
         avecl = 0.0D+00
       endif

       do k = 1, 3
         k1(k) = b1(k) * dr
         xyz3(k) = xyz1(k) + k1(k) * 0.5
       enddo
       call xyz2rtp(xyz3,rtp3)
       call vecatrtp(b2,rtp3)
       call vrtp2xyz(b2,rtp3)
       call normvec3(b2,bbb)
       do k = 1, 3
         k2(k) = b2(k) * dr
         xyz3(k) = xyz1(k) + k2(k) * 0.5
       enddo
       call xyz2rtp(xyz3,rtp3)
       call vecatrtp(b3,rtp3)
       call vrtp2xyz(b3,rtp3)
       call normvec3(b3,bbb)
       do k = 1, 3
         k3(k) = b3(k) * dr
         xyz3(k) = xyz1(k) + k3(k)
       enddo
       call xyz2rtp(xyz3,rtp3)
       call vecatrtp(b4,rtp3)
       call vrtp2xyz(b4,rtp3)
       call normvec3(b4,bbb)
       do k = 1, 3
         k4(k) = b4(k) * dr
         b0(k) = (b1(k) + 2.0 * b2(k) + 2.0 * b3(k) + b4(k)) / 6.0
       enddo
       call  normvec3(b0,bbb)
       do k = 1, 3
         k0(k) = b0(k) * dr
         xyz2(k) = xyz1(k) + k0(k)
       enddo

       call xyz2rtp(xyz2,rtp2)

       return
       end



*
** --------------
*
       subroutine normvec3(xyzvec, lngth)
       implicit none
* interface
       real*8  xyzvec(3), lngth
* local

       lngth = xyzvec(1)**2 + xyzvec(2)**2 + xyzvec(3)**2
       lngth = dsqrt(lngth)
       if (lngth .GT. 1.0D-40) then
         xyzvec(1) = xyzvec(1) / lngth
         xyzvec(2) = xyzvec(2) / lngth
         xyzvec(3) = xyzvec(3) / lngth
       else
         xyzvec(1) = 0.0D+00
         xyzvec(2) = 0.0D+00
         xyzvec(3) = 0.0D+00
       endif

       return
       end

*
** --------------
*
       subroutine rtp2xyz(rtp,xyz)
       implicit none
* interface
       real*8  rtp(3), xyz(3)

       xyz(1) = rtp(1) * dsin(rtp(2)) * dcos(rtp(3))
       xyz(2) = rtp(1) * dsin(rtp(2)) * dsin(rtp(3))
       xyz(3) = rtp(1) * dcos(rtp(2))

       return
       end

*
** -----------------------------------------------------------------
*
       subroutine vrtp2xyz(vec,rtp)
       implicit none
       real*8  vec(3), rtp(3)
       real*8  ra, th, ph, vec2(3)

       ra = rtp(1)
       th = rtp(2)
       ph = rtp(3)

       vec2(1) = vec(1) * dsin(th) * dcos(ph) ! UNDER
     &         + vec(2) * dcos(th) * dcos(ph)
     &         - vec(3)            * dsin(ph)
       vec2(2) = vec(1) * dsin(th) * dsin(ph)
     &         + vec(2) * dcos(th) * dsin(ph)
     &         + vec(3)            * dcos(ph)
       vec2(3) = vec(1) * dcos(th)
     &         - vec(2) * dsin(th)

       vec(1) = vec2(1)
       vec(2) = vec2(2)
       vec(3) = vec2(3)

       return
       end


*
** --------------
*
       subroutine xyz2rtp(xyza,rtp)
       implicit none
* interface
       real*8  xyza(3), rtp(3)
* local
       real*8  aa, ra, th, ph
       real*8  pi
       parameter(pi = 3.14159265358979D+00)

       ra = dsqrt(xyza(1)**2 + xyza(2)**2 + xyza(3)**2)
       th = dacos(xyza(3) / ra)
       aa = dsqrt(xyza(1)**2 + xyza(2)**2)
       if (aa .LT. 1.0D-03) aa = 1.0D-03
       ph = dacos(xyza(1) / aa)
       if (xyza(2) .LT. 0.0) ph = - ph + 2.0D+00 * pi
       if (ph      .LT. 0.0) ph =   ph + 2.0D+00 * pi
       rtp(1) = ra
       rtp(2) = th
       rtp(3) = ph

       return
       end


*
* magnetic field at (ra,th,ph) ---------------------------
*
       subroutine vecatrtp(vec,rtp)
       implicit none
* interface
       real*8  vec(3), rtp(3)
* local
       real*8  ra ,th, ph
       real*8  aa, bb, cc
       real*8  dra1, dra2, dth1, dth2, dph1, dph2
       integer i1, i2, j1, j2, k1, k2, ijk(3)
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common  /vecvec/ vcr, vct, vcp
       common  /sphcrd/  rr, theta, phi

       ra = rtp(1)
       th = rtp(2)
       ph = rtp(3)
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
       aa =(vcr(i1,j1,k1)*dth2*dph2+vcr(i1,j2,k1)*dth1*dph2
     &    + vcr(i1,j1,k2)*dth2*dph1+vcr(i1,j2,k2)*dth1*dph1)*dra2
     &    +(vcr(i2,j1,k1)*dth2*dph2+vcr(i2,j2,k1)*dth1*dph2
     &    + vcr(i2,j1,k2)*dth2*dph1+vcr(i2,j2,k2)*dth1*dph1)*dra1
       bb =(vct(i1,j1,k1)*dth2*dph2+vct(i1,j2,k1)*dth1*dph2
     &    + vct(i1,j1,k2)*dth2*dph1+vct(i1,j2,k2)*dth1*dph1)*dra2
     &    +(vct(i2,j1,k1)*dth2*dph2+vct(i2,j2,k1)*dth1*dph2
     &    + vct(i2,j1,k2)*dth2*dph1+vct(i2,j2,k2)*dth1*dph1)*dra1
       cc =(vcp(i1,j1,k1)*dth2*dph2+vcp(i1,j2,k1)*dth1*dph2
     &    + vcp(i1,j1,k2)*dth2*dph1+vcp(i1,j2,k2)*dth1*dph1)*dra2
     &    +(vcp(i2,j1,k1)*dth2*dph2+vcp(i2,j2,k1)*dth1*dph2
     &    + vcp(i2,j1,k2)*dth2*dph1+vcp(i2,j2,k2)*dth1*dph1)*dra1

       vec(1) = aa
       vec(2) = bb
       vec(3) = cc

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
       common /sphcrd/  rr, theta, phi

       if (rtp(1) .LT. rr(0)) then
         ijk(1) = 0
       else if (rtp(1) .GT. rr(ii)) then
         ijk(1) = ii
       else
         do m = 0, ii
           if ((rtp(1).GE.rr(m)).AND.(rtp(1).LE.rr(m+1))) ijk(1) = m
         enddo
       endif

       if (rtp(2) .LT. 0.0) rtp(2) = - rtp(2)
       if (rtp(2) .GT.  pi) rtp(2) = - rtp(2) + 2.0D+00 * pi
       do m = 0, jj - 1
         if ((rtp(2).GE.theta(m  )).AND.
     &       (rtp(2).LE.theta(m+1))) ijk(2) = m
       enddo

       if (rtp(3) .LE. 0.0D+00)      rtp(3) = rtp(3) + 2.0D+00 * pi
       if (rtp(3) .GT. 2.0D+00 * pi) rtp(3) = rtp(3) - 2.0D+00 * pi
       do m = 0, kk
         if ((rtp(3).GE.phi(m)).AND.
     &       (rtp(3).LE.phi(m+1))) ijk(3) = m
       enddo

       return
       end

*
* ----------------------------------------------------
*
       real*8 function lrntatjk(i,j,k,omega,ivar2,ir)
       implicit none
* interface
       integer ivar2, i, j, k, ir
       real*8  omega
* local
       real*8  ans, aa, bb, cc, dd
       real*8  a1, a2, a3, alpha
       real*8  valatjk
       integer ivar0
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  var8(0:ix,0:jy,8)
       common /sphcrd/ rr, theta, phi
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0
       common /var2d/  var8

       alpha =(rr(ir+1) - rr(ir)) / (rr(ir) - rr(ir-1))
       a1 =          1.0D+00 / alpha
       a2 =  alpha - 1.0D+00 / alpha
       a3 = -alpha

       if (ivar2 .EQ. 1) then !      rot(B)xB_r
         ivar0 = 13
         aa = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 =  8
         bb = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 = 14
         cc = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 =  7
         dd = valatjk(i,j,k,omega,ivar0,ir)
         ans = aa * bb - cc * dd
       else if (ivar2 .EQ. 2) then ! rot(B)xB_theta
         ivar0 = 14
         aa = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 =  6
         bb = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 = 12
         cc = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 =  8
         dd = valatjk(i,j,k,omega,ivar0,ir)
         ans = aa * bb - cc * dd
       else if (ivar2 .EQ. 3) then ! rot(B)xB_phi
         ivar0 = 12
         aa = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 =  7
         bb = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 = 13
         cc = valatjk(i,j,k,omega,ivar0,ir)
         ivar0 =  6
         dd = valatjk(i,j,k,omega,ivar0,ir)
         ans = aa * bb - cc * dd
       else
         write(*,*) ' wrong at Lorentz Force Func'
         write(*,*) ivar2
         stop
       endif

       lrntatjk = ans

       return
       end

*
* ----------------------------------------------------
*
       real*8 function valatjk(i,j,k,omega,ivar,ir)
       implicit none
* interface
       integer ivar, i, j, k, ir
       real*8  omega
*local
       real*8  rot1, rot2, sss
       real*8  ans
       real*8  a1, a2, a3, alpha
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  rr(0:ii+1), theta(0:jj),phi(0:kk+1)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  var8(0:ix,0:jy,8)
       common /sphcrd/ rr, theta, phi
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0
       common /var2d/  var8

       if ((i .GT. 0) .AND. (i .LE. ii)) then
         alpha =(rr(ir+1) - rr(ir)) / (rr(ir) - rr(ir-1))
         a1 =          1.0D+00 / alpha
         a2 =  alpha - 1.0D+00 / alpha
         a3 = -alpha
       else
         alpha = 1.0 ! dummy
         a1    = 1.0
         a2    = 1.0
         a3    = 1.0
       endif

*       omega = 0.0

* eight variables
       if (ivar .LE. 8) then
         ans = var8(k,j,ivar)
         if (ivar .EQ. 5) ans = ans + rr(ir)*dsin(theta(j))*omega
       else if (ivar .EQ. 42) then
         ans = var8(k,j,2) / var8(k,j,1)
* rot flow in the rest frame
       else if (ivar .EQ. 9) then ! rot(V)_r
         rot1 =(up0(i,j+1,k) * dsin(theta(j+1))
     &         -up0(i,j-1,k) * dsin(theta(j-1)))
     &        /(theta(j+1) - theta(j-1)) / dsin(theta(j)) / rr(ir)
         sss  = 2.0D+00 * omega * dcos(theta(j))
         rot2 =(ut0(i,j,k+1)
     &         -ut0(i,j,k-1))
     &        /(phi(k+1)-phi(k-1)) / dsin(theta(j)) / rr(ir)
         ans = rot1 - rot2 + sss
       else if (ivar .EQ. 10) then ! rot(V)_theta
         rot1 =(ur0(i,j,k+1)
     &         -ur0(i,j,k-1))
     &        /(phi(k+1)-phi(k-1)) / dsin(theta(j)) / rr(ir)
         rot2 =(up0(i+1,j,k) * rr(ir+1) * a1
     &         +up0(i  ,j,k) * rr(ir  ) * a2
     &         +up0(i-1,j,k) * rr(ir-1) * a3)
     &        /(rr(ir+1) - rr(ir-1)) / rr(ir)
         sss  = 2.0D+00 * omega * dsin(theta(j))
         ans = rot1 - rot2 - sss
       else if (ivar .EQ. 11) then ! rot(V)_phi
         rot1 =(ut0(i+1,j,k) * rr(ir+1)*a1
     &         +ut0(i  ,j,k) * rr(ir  )*a2
     &         +ut0(i-1,j,k) * rr(ir-1)*a3)
     &        /(rr(ir+1) - rr(ir-1)) / rr(ir)
         rot2 =(ur0(i,j+1,k)
     &         -ur0(i,j-1,k))
     &        /(theta(j+1) - theta(j-1)) / rr(ir)
         ans = rot1 - rot2
* rot mag
       else if (ivar .EQ. 12) then ! rot(B)_r
         rot1 =(bp0(i,j+1,k) * dsin(theta(j+1))
     &         -bp0(i,j-1,k) * dsin(theta(j-1)))
     &        /(theta(j+1) - theta(j-1)) / dsin(theta(j)) / rr(ir)
         rot2 =(bt0(i,j,k+1)
     &         -bt0(i,j,k-1))
     &        /(phi(k+1)-phi(k-1)) / dsin(theta(j)) / rr(ir)
         ans = rot1 - rot2
       else if (ivar .EQ. 13) then ! rot(B)_theta
         rot1 =(br0(i,j,k+1)
     &         -br0(i,j,k-1))
     &        /(phi(k+1)-phi(k-1)) / dsin(theta(j)) / rr(ir)
         rot2 =(bp0(i+1,j,k) * rr(ir+1) * a1
     &         +bp0(i  ,j,k) * rr(ir  ) * a2
     &         +bp0(i-1,j,k) * rr(ir-1) * a3)
     &        /(rr(ir+1) - rr(ir-1)) / rr(ir)
         ans = rot1 - rot2
       else if (ivar .EQ. 14) then ! rot(B)_phi
         rot1 =(bt0(i+1,j,k) * rr(ir+1) * a1
     &         +bt0(i  ,j,k) * rr(ir  ) * a2
     &         +bt0(i-1,j,k) * rr(ir-1) * a3)
     &        /(rr(ir+1) - rr(ir-1)) / rr(ir)
         rot2 =(br0(i,j+1,k)
     &         -br0(i,j-1,k))
     &        /(theta(j+1) - theta(j-1))/rr(ir)
         ans = rot1 - rot2
* rot VxB in the rotating frame
       else if (ivar .EQ. 15) then ! rot(VxB)_r
         rot1 =((ur0(i,j+1,k) * bt0(i,j+1,k)
     &          -ut0(i,j+1,k) * br0(i,j+1,k))* dsin(theta(j+1))
     &         -(ur0(i,j-1,k) * bt0(i,j-1,k)
     &          -ut0(i,j-1,k) * br0(i,j-1,k))* dsin(theta(j-1)))
     &        /(theta(j+1) - theta(j-1)) / dsin(theta(j)) / rr(ir)
         rot2 =((up0(i,j,k+1) * br0(i,j,k+1)
     &          -ur0(i,j,k+1) * bp0(i,j,k+1))
     &         -(up0(i,j,k-1) * br0(i,j,k-1)
     &          -ur0(i,j,k-1) * bp0(i,j,k-1)))
     &        /(phi(k+1)-phi(k-1)) / dsin(theta(j)) / rr(ir)
         ans = rot1 - rot2
       else if (ivar .EQ. 16) then ! rot(VxB)_theta
         rot1 =((ut0(i,j,k+1) * bp0(i,j,k+1)
     &          -up0(i,j,k+1) * bt0(i,j,k+1))
     &         -(ut0(i,j,k-1) * bp0(i,j,k-1)
     &          -up0(i,j,k-1) * bt0(i,j,k-1)))
     &        /(phi(k+1)-phi(k-1)) / dsin(theta(j)) / rr(ir)
         rot2 =((ur0(i+1,j,k) * bt0(i+1,j,k)
     &          -ut0(i+1,j,k) * br0(i+1,j,k)) * rr(ir+1) * a1
     &         +(ur0(i  ,j,k) * bt0(i  ,j,k)
     &          -ut0(i  ,j,k) * br0(i  ,j,k)) * rr(ir  ) * a2
     &         +(ur0(i-1,j,k) * bt0(i-1,j,k)
     &          -ut0(i-1,j,k) * br0(i-1,j,k)) * rr(ir-1) * a3)
     &        /(rr(ir+1) - rr(ir-1)) / rr(ir)
         ans = rot1 - rot2
       else if (ivar .EQ. 17) then ! rot(VxB)_phi
         rot1 =((up0(i+1,j,k) * br0(i+1,j,k)
     &          -ur0(i+1,j,k) * bp0(i+1,j,k)) * rr(ir+1) * a1
     &         +(up0(i  ,j,k) * br0(i  ,j,k)
     &          -ur0(i  ,j,k) * bp0(i  ,j,k)) * rr(ir  ) * a2
     &         +(up0(i-1,j,k) * br0(i-1,j,k)
     &          -ur0(i-1,j,k) * bp0(i-1,j,k)) * rr(ir-1) * a3)
     &        /(rr(ir+1) - rr(ir-1)) / rr(ir)
         rot2 =((ut0(i,j+1,k) * bp0(i,j+1,k)
     &          -up0(i,j+1,k) * bt0(i,j+1,k))
     &         -(ut0(i,j-1,k) * bp0(i,j-1,k)
     &          -up0(i,j-1,k) * bt0(i,j-1,k)))
     &        /(theta(j+1) - theta(j-1)) / rr(ir)
         ans = rot1 - rot2
* VxB
       else if (ivar .EQ. 22) then ! (VxB)_r
         ans = ut0(i,j,k)*bp0(i,j,k) - up0(i,j,k)*bt0(i,j,k)
       else if (ivar .EQ. 23) then ! (VxB)_theta
         ans = up0(i,j,k)*br0(i,j,k) - ur0(i,j,k)*bp0(i,j,k)
       else if (ivar .EQ. 24) then ! (VxB)_phi
         ans = ur0(i,j,k)*bt0(i,j,k) - ut0(i,j,k)*br0(i,j,k)
* invalid index
       else
         write(*,*)  'Ivar is not available at valatjk .. stop'
         write(13,*) 'Ivar is not available at valatjk .. stop'
         stop
       endif

       valatjk = ans

       return
       end


* ----------------
*
       subroutine readmhd(rt,nnn,lmhdfile)
       implicit none
       integer rt, nnn
       logical lmhdfile
*
       real*8  vv(8)
       integer i, j, k, l
       integer idummy1, jdummy1, kdummy1
       integer idummy2, jdummy2, kdummy2
       logical ldummy
       character*64 flname
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  ut2(0:ii+1,0:jj,0:kk+1),bt2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common /vecvec/ vcr, vct, vcp
       common /var3do/ ro2, pg2, ur2, ut2, up2, br2, bt2, bp2

       write(flname,'(''d'',i6.6,''.'',i4)') nnn + 300000, rt
       inquire(file = flname, exist = ldummy)
       if (ldummy) then
         open(unit=1,file=flname,status='old')
 100     continue
           read(1,'(3i5,5(1x,e11.5),1x,3i5,3(1x,e11.5))',END=99)
     &       idummy1, jdummy1, kdummy1, (vv(l),l=1,5),
     &       idummy2, jdummy2, kdummy2, (vv(l),l=6,8)
           i = idummy1 - 1000
           j = jdummy1 - 2000
           k = kdummy1 - 3000
           ro2(i,j,k) = vv(1)
           pg2(i,j,k) = vv(2)
           ur2(i,j,k) = vv(3)
           ut2(i,j,k) = vv(4)
           up2(i,j,k) = vv(5)
           br2(i,j,k) = vv(6)
           bt2(i,j,k) = vv(7)
           bp2(i,j,k) = vv(8)
           goto 100
 99      continue
         close(1)
         lmhdfile = .true.
       else
         write(flname,'(''u'',i6.6,''.'',i4)') nnn + 300000, rt
         inquire(file = flname, exist = ldummy)
         if (ldummy) then
           open(unit=1,file=flname,status='old',form='unformatted')
           read(1) ro2
           read(1) pg2
           read(1) ur2
           read(1) ut2
           read(1) up2
           read(1) br2
           read(1) bt2
           read(1) bp2
           close(1)
           lmhdfile = .true.
         else
           write(*,*) ' File not found ',flname
           lmhdfile = .false.
         endif
       endif

       do 101 k = 0, kk + 1
       do 101 j = 0, jj
       do 101 i = 0, ii + 1
         vcr(i,j,k) = br2(i,j,k)
         vct(i,j,k) = bt2(i,j,k)
         vcp(i,j,k) = bp2(i,j,k)
 101   continue

       return
       end


* -------------------------
*
       subroutine mhd2alt(ir,lmhdfile)
       implicit none
       integer ir
       logical lmhdfile
*
       integer i, j, k, l, i2
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 64, kk = 128)
       integer ix, jy
       parameter(ix = kk, jy = jj)
       real*8  var8(0:ix,0:jy,8)
       real*8  ro0(-1:1,0:jj,0:kk+1),pg0(-1:1,0:jj,0:kk+1)
       real*8  ur0(-1:1,0:jj,0:kk+1),br0(-1:1,0:jj,0:kk+1)
       real*8  ut0(-1:1,0:jj,0:kk+1),bt0(-1:1,0:jj,0:kk+1)
       real*8  up0(-1:1,0:jj,0:kk+1),bp0(-1:1,0:jj,0:kk+1)
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  ut2(0:ii+1,0:jj,0:kk+1),bt2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)
       common /var2d/  var8
       common /var3d/  ro0, pg0, ur0, ut0, up0, br0, bt0, bp0
       common /var3do/ ro2, pg2, ur2, ut2, up2, br2, bt2, bp2

       if (lmhdfile) then

         do 100 k = 0, kk
         do 100 j = 0, jj
           var8(k,j,1) = ro2(ir,j,k)
           var8(k,j,2) = pg2(ir,j,k)
           var8(k,j,3) = ur2(ir,j,k)
           var8(k,j,4) = ut2(ir,j,k)
           var8(k,j,5) = up2(ir,j,k)
           var8(k,j,6) = br2(ir,j,k)
           var8(k,j,7) = bt2(ir,j,k)
           var8(k,j,8) = bp2(ir,j,k)
 100     continue

         do 101 k = 0, kk + 1
         do 101 j = 0, jj
         do 101 i2 = - 1, 1
           i = ir + i2
           if (i .GT. ii + 1) i = ii + 1
           if (i .LT. 0) i = 0
           ro0(i2,j,k) = ro2(i,j,k)
           pg0(i2,j,k) = pg2(i,j,k)
           ur0(i2,j,k) = ur2(i,j,k)
           ut0(i2,j,k) = ut2(i,j,k)
           up0(i2,j,k) = up2(i,j,k)
           br0(i2,j,k) = br2(i,j,k)
           bt0(i2,j,k) = bt2(i,j,k)
           bp0(i2,j,k) = bp2(i,j,k)
 101     continue

       else

         do 200 l = 1, 8
         do 200 k = 0, kk
         do 200 j = 0, jj
           var8(k,j,l) = 0.0D+00
 200     continue

         do 201 k = 0, kk + 1
         do 201 j = 0, jj
         do 201 i2 = - 1, 1
           ro0(i2,j,k) = 0.0D+00
           pg0(i2,j,k) = 0.0D+00
           ur0(i2,j,k) = 0.0D+00
           ut0(i2,j,k) = 0.0D+00
           up0(i2,j,k) = 0.0D+00
           br0(i2,j,k) = 0.0D+00
           bt0(i2,j,k) = 0.0D+00
           bp0(i2,j,k) = 0.0D+00
 201     continue

       endif

       return
       end

