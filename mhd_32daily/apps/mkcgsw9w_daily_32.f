**********************************************************************
*
* mkcgswX <= mvflow3, mvmhdf1, mvf2de13, mkrgbcsg
*
* MaKe CG of Solar Wind
*
* fortran program for making CG
*           contour of scalar & field line of vector
*
* INPUT  file  : d3NNNNN.dat or d3NNNNN.RRRR, init.dat
* OUTPUT file  : wFFFRRRR.(rgb fal hex ppm pgm and ps), mkcgswX.log
*
* 1-point view (Itten-tohshi hou)
*
** Example of usage to make movie on NVS2000 system
*
*    % a.out < incg.lst | rgb2melt -t3 -x648 -y486 | melt2nvs - >> a.nvs
*
*
** Example of 'incg.lst'
*
* 3           : index for output device
*       negative    invalid (vacant for future use)
*          0        RGB     (std out) : for movie
*          1        RGB+FAL (file)
*          2        PPM P3  (file)    : usually not used
*          3        PPM P6  (file)    : portable format
*          4        HEX     (file)    : usually used at Fujitsu or SunMicro
*          5        PGM P5  (file)    : gray (+ negative)
*          6        PS      (file)    : rgb
*          7        PS      (file)    : gray (+ negative)
*          8        PS      (file)    : cmyk
*          9        BMP     (file)    : bgr
* 2           : 1=yes to draw grid : usually non-1
* 12          : num. of frames
* 1653 10000 40.0 0 1 1 0 0 0 60.0 180.0 0 2
*          (continue)
*
*  1) 1653        : rot. num.
*  2) 60000       : step num. of comp.
*  3) 40.0        : radii of view point (float)
*  4) compress radial distance  (1 = Yes, other = No)
*  5) show solar surface        (1 = photo.mag, 2 = yellow ball,
*                                3 = flowspeed, other = No)
*  6) show HCS                  (1 = Yes, other = No)
*  7) show Mag line             (1 = Yes, other = No)
*  8) show Flow line            (1 = Yes, other = No)
*  9) show Integ/scalar         (1--12 = Yes, other = No)
* 10) latitude of eye point (float, degree, '0' is for North pole)
* 11) longitude of eye point
* 12) string or not etc. (1=Pixel, -1=Ps, other=No)
* 13) write IPS obs line : 1 = white, 2 = color & wreight, other = no
*
** The parameter (or constant) can be changed for purposes
*
*   ifine          : degree of fineness of image
*   ismthe         : degree of gradualness between two images
*   iframe jframe  : number of images in one page
*   bkrgb(3)       : background color
*   lthdeg, lphdeg : direction of light source
*   irtube         : index for radius of tube
*
*   dproj          : if quite large value were taken,
*                       the output image approach "parallel-projection"
*
** On some changes that cannot be controlled by parameters.
*  1) seek word C2 to make larger ocult images (of LoS)
*  2) seek words "Skip near Axis" to avoid drawing line near polar axis
*  3) modify logical function lnotdrwn() (not) to draw things at particular potision
*  4) seek word "cinema-mode" to change final pixel sizes of PPM.
*  5) seek words " draw only closed line" to draw only closed lines
*
*                                 K.Hayashi    Apr, 2005
*
**********************************************************************

       program mkcgsw
       implicit none
*
       integer nrot, nnn, numframe, numpage, onepage, numfram2
       integer numend
       integer nref, rotref
       integer isunsurf, iplthcs, ipltmag, iprtcl
       integer ipltflow, ipltrho, igrid, ipltrhox, ienhintg
       integer is, ie, js, je, ks, ke
       real*8  radorb, radorb2
*
       integer i0, j0, k0
*
       real*8  rupper0, seeth0, seeph0, seeth, seeph
*
       integer bkrgb(3)  ! /0, 0, 0/  ! color of back ground
*
       real*8  gamma, t0, length0, v0, omega
       real*8  rau
       real*8  aa, bb
       integer i, j, k
*
       character*20 charout
       integer charlen, ix, jy, itrape, ixfp0, jyfp0
       integer ivar, irtube, iwristr, ihist, iips
* function
       logical findfile
       real*8  lograd
* dummy and flag
       integer idummy1, idummy2, idummy3, jdummy1
       logical ldummy
       logical leyearth
       integer itubecol
       integer nnn2, nrot2, ipltflw2
* RGB data storage
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0*ifine, jvpix2 = jvpix0*ifine)
       character*1 rgb0(3,ihpix0, jvpix0) ! RGB data for present frame
       character*1 rgb2(3,ihpix2, jvpix2) ! RGB data for fine picture
       character*1 rgbx(3,ihpix2, jvpix2) ! RGB data for present frame : sub array
       integer ihpixf, jvpixf, iframe, jframe, maxfatp
       parameter(iframe = 1, jframe = 1)
       parameter(maxfatp = iframe * jframe)
       parameter(ihpixf = ihpix0*iframe, jvpixf = jvpix0*jframe)
       character*1 rgblarge(3,ihpixf,jvpixf) ! RGB data for present large frame
       character*1 rgblargp(3,ihpixf,jvpixf) ! RGB data for previous large frame
*
       character*33 psstr(maxfatp)
       character sdummy1*18 ! , sdummy2*15
       integer numpstr, psix(maxfatp), psjy(maxfatp)
* const
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* common
       real*8  dstnt(ihpix2, jvpix2) ! distance projection plane <-> objects
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  ro0(0:ii+1,0:jj,0:kk+1)
       real*8  pg0(0:ii+1,0:jj,0:kk+1)
       real*8  ur0(0:ii+1,0:jj,0:kk+1)
       real*8  ut0(0:ii+1,0:jj,0:kk+1)
       real*8  up0(0:ii+1,0:jj,0:kk+1)
       real*8  br0(0:ii+1,0:jj,0:kk+1)
       real*8  bt0(0:ii+1,0:jj,0:kk+1)
       real*8  bp0(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       real*8  ro3(0:ii+1,0:jj,0:kk+1)
       real*8  pg3(0:ii+1,0:jj,0:kk+1)
       real*8  ur3(0:ii+1,0:jj,0:kk+1)
       real*8  ut3(0:ii+1,0:jj,0:kk+1)
       real*8  up3(0:ii+1,0:jj,0:kk+1)
       real*8  br3(0:ii+1,0:jj,0:kk+1)
       real*8  bt3(0:ii+1,0:jj,0:kk+1)
       real*8  bp3(0:ii+1,0:jj,0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       integer nfoot
       parameter(nfoot = 10000) ! approxiamte..
       logical lopen1st(nfoot)
       integer prtmax
*       parameter(prtmax = 600)
       parameter(prtmax = 200)
       real*8  prtxyz(3,prtmax)
       logical ldrwprt(prtmax)
       integer iout, irlog
       common  /lopline/ lopen1st
       common  /prtprt/ prtxyz
       common  /lprtprt/ ldrwprt
       common  /images/ dstnt, rgb2
       common  /coord/  rr, theta, phi
       common  /var01/  ro0, pg0
       common  /var02/  ur0, ut0, up0
       common  /var04/  br0, bt0, bp0
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /var31/  ro3, pg3
       common  /var32/  ur3, ut3, up3
       common  /var34/  br3, bt3, bp3
       common  /sclscl/ sclr
       common  /vecvec/ vcr, vct, vcp
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom
       common  /cgcfg/  iout, irlog
       common  /gasgas/ gamma

* initialize
       do idummy1 =  1, nfoot
         lopen1st(idummy1) = .false.
       enddo

* open log_file
       open(unit = 11, file = 'mkcgsw9w.log', status = 'unknown')

* write out parameter(s)
       write(11,'('' Ifine  = '',i4)') ifine
       write(11,'('' Iframe = '',i4)') iframe
       write(11,'('' Jframe = '',i4)') jframe

* select out put format and device
*       negative    invalid (vacant for future use)
*          0        RGB     (std out) : for movie
*          1        RGB+FAL (file)
*          2        PPM P3  (file)    : usually not used
*          3        PPM P6  (file)    : portable format
*          4        HEX     (file)    : usually not used
*          5        PGM P5  (file)    : gray (+ negative)
*          6        PS      (file)    : rgb
*          7        PS      (file)    : gray (+ negative)
*          8        PS      (file)    : cmyk
*          9        BMP     (file)    : bgr

       read(*,*) iout
       if ((iout .LT. 0) .OR. (iout .GT. 9)) then
         write(*,*) 'Output device selection is invalid'
         stop
       endif
       if (iout .NE. 0) then
         write(*,'('' Iout = '', i4)') iout
         write(*,'(A)') ' 0   = RGB (stdout), 1 = RGB and FAL'
         write(*,'(A)') ' 2,3 = PPM (P3/P6),  4 = Hex        '//
     &                  ' 5   = PGM (P5)      9 = 24bit-BMP'
         write(*,'(A)') ' 6-8 = PS (rgb, negative grey,cmyk)'
       endif

       if (iout .GT. 0)
     &   write(*,*) 'INPUT Index : Draw Grid    (1/other)=(Y/N)'
       read(*,*) igrid
       if (igrid .EQ. 1) then
         if (iout .GT. 0) then
           write(*,*) 'INPUT Grid Address'
           write(*,*) ' Is Ie Js Je Ks Ke'
         endif
         read(*,*) is, ie, js, je, ks, ke
       else
         is = 0
         ie = 0
         js = 0
         je = 0
         ks = 0
         ke = 0
       endif

       if (iout .GT. 0)
     &   write(*,*) 'INPUT Int. : Histgram Re-Color (1/other)=(Y/N)'
       read(*,*) ihist

       if (iout .GT. 0) write(*,*) 'INPUT Num of image'
       read(*,*) numend

       if (iout .GT. 0) then
         write(*,*) 'INPUT 19 variables'
         write(*,*) 'Nrot NNN Rupper(R_sun) Irlog (1/other=Y/N)'
         write(*,*) 'Isunsurf(1-7)=(Br,Yellow,Vr,T,N,P,O/C-2.5Rs,none)'
         write(*,*) 'Ihcs (1/2/other=Yes/shadowed/No)'
         write(*,*) 'Imag 7,1,--6,8,others : gray,green,STEL-Vr'
         write(*,*) '   Normal-Vr, open/close, O/C-first,eqBack,'
         write(*,*) '   std(T)_i(or Max/Min T) ,None'
         write(*,*) 'Iflow 1,2,3,4(-1,,-4 : - r sin omega),other'
         write(*,*) '   = blue,color(eqBack,ful,eqOutward),No'
         write(*,*) 'Irho : 1--12 = surf_div(rho), intg(pB)'
         write(*,*) '   surf_div(v), intg(rho/rho_prk), intg(v_t rho)'
         write(*,*) '   intg(ro^2 r^4), intg(std(1/v)_i)'
         write(*,*) '   intg(std(rho)_i), intg(1/(V*std(rho)_i)_i)'
         write(*,*) '   intg(rho), intg(Pg/Pg_prk)'
         write(*,*) '   surf(div(rho)_i=0.0), >= 13 : none '
         write(*,*) '   0 or negative = Nref of intg(rho/rho_Nref)'
         write(*,*) ' see_theta see_phi'
         write(*,*) ' Draw particle (1 = Yes, other = No)'
         write(*,*) ' IWriteString (1=Pixel, -1=Ps,   other=No)'
         write(*,*) ' IDrawIpsLine (1=White,  2=Color,other=No)'
         write(*,*) ' Iemphasize for Intg (1=Circle'
         write(*,*) '      2=DevEntire,3=Max/MinEntire'
         write(*,*) '      4=0-Dev,5=Histgram'
         write(*,*) '      6=ParkerDensity(r=linear)'
         write(*,*) '      7=Newkirk filter, other = None'
         write(*,*) ' 3 integer for Background RGB (0-255)'
         write(*,*) ' Radius(Rs) of orbit circle : negative = No'
       endif

       numpage = 0
       onepage = 0
       numpstr = 0

       nnn2  = -10
       nrot2 = -10

       do numfram2=1, numend

       numframe = numfram2

* read data
         read(*,*)
     &     nrot, nnn, rupper0, irlog,
     &     isunsurf, iplthcs, ipltmag, ipltflow, ipltrho,
     &     seeth0, seeph0, iprtcl, iwristr, iips, ienhintg,
     &     (bkrgb(idummy1),idummy1=1,3), radorb

         if (ipltrho .LE. 0) then
           ipltrhox = -1
         else
           ipltrhox = ipltrho
         endif
         if (iout .GT. 0) then
           write(*,'(A,$)') ' Nf Np   RT    NNN Rupper See_th See_ph'
           write(*,'(A)')   '/ Irad Isun Ihcs Imag Iflw Irho Iips'
           write(*,'(2i3,i5,i7,3f7.2,''/'',7i5)')
     &        mod(numframe,1000), mod(numpage+1,1000),
     &        nrot, nnn,rupper0,seeth0,seeph0,
     &        irlog,isunsurf,iplthcs,ipltmag,ipltflow,ipltrhox,iips
         endif

         open(12,file='carr.txt',status='old')
         read(12,*) i ! ncar
         read(12,*) j
         close(12)
         seeph0 = float(j)
         open(12,file='latlon.txt',status='old')
         read(12,*) aa, bb ! lat and lon
         close(12)
         seeth0 = aa

         write(11,'(A,$)') ' Nf Np   RT    NNN Rupper See_th See_ph'
         write(11,'(A)')   '/ Irad Isun Ihcs Imag Iflw Irho Iips'
         write(11,'(2i3,i5,i7,3f7.2,''/'',7i5)')
     &      numframe, numpage+1, nrot, nnn,rupper0,seeth0,seeph0,
     &      irlog,isunsurf,iplthcs,ipltmag,ipltflow,ipltrhox,iips

* initializing ......
         call initrgb(rgb2,dstnt,bkrgb)

* check data file exist or not
         ldummy = findfile(nrot,nnn)

** construct one frame
         if (ldummy) then

           if ((nrot .NE. nrot2) .OR. (nnn .NE. nnn2)) then
             call readdata(nrot,nnn) ! XX3 be fillled
             do k = 0, kk + 1
             do j = 0, jj
             do i = 0, ii + 1
               ro1(i,j,k) = ro3(i,j,k)
               pg1(i,j,k) = pg3(i,j,k)
               ur1(i,j,k) = ur3(i,j,k)
               ut1(i,j,k) = ut3(i,j,k)
               up1(i,j,k) = up3(i,j,k)
               br1(i,j,k) = br3(i,j,k)
               bt1(i,j,k) = bt3(i,j,k)
               bp1(i,j,k) = bp3(i,j,k)
             enddo
             enddo
             enddo
             if (iout .GT. 0)
     &          write(*,'('' Read .. done : '',2i8)') nrot, nnn
           endif
           if ((ipltflow .LE. -1) .AND. (ipltflow .GE. -4)) then
             call readinit(gamma,t0,length0,v0,omega)
             do k = 0, kk + 1
             do j = 0, jj
             do i = 0, ii + 1
               up1(i,j,k) = up1(i,j,k)-rr(i)*dsin(theta(j))*omega
             enddo
             enddo
             enddo
             ipltflow = - ipltflow
             if (iout .GT. 0)
     &          write(*,*) ' R sin(theta) omega correction done'
           endif
           nnn2  = nnn
           nrot2 = nrot
           ipltflw2 = ipltflow

           call readinit(gamma,t0,length0,v0,omega)

           rbottom = rr(0)
           rupper  = rupper0 / length0 ! R_sun => non-dim. length
           seeth = seeth0
           seeph = seeph0

* reset coordinate if neccesary
           if (irlog .EQ. 1) then
             call resetrad(rr, rbottom, rupper, iout)
             if (iout .GT. 0)
     &         write(*,'('' R : '',2f10.2)') rr(0), rr(ii+1)
           endif

           leyearth = .false. ! whether eye is at the Earth or not

           rau = 215.0D+00 / length0
           call setgrap(xyzcp,dproj,xyzxp,xyzyp,xyzlght,
     &                  rupper,rbottom,seeth,seeph,
     &                  ul,numframe,numend,leyearth,rau)

* set variables for eye, light source and projection plane

* surface of Sun
           if ((isunsurf .GE. 1) .AND. (isunsurf .LE. 6)) then
             call surf2rgb(isunsurf,v0)
             if (iout .GT. 0) write(*,*) 'SURF2RGB ..end'
             write(11,*) 'SURF2RGB ..end'
           else if (isunsurf .EQ. 7) then
             call ocsf2rgb()
             if (iout .GT. 0) write(*,*) 'Open/Close SURF2RGB ..end'
             write(11,*) 'Open/Close SURF2RGB ..end'
           endif

           if (radorb .GT. 0.0D+00) then
             radorb2 = radorb / length0
             if (irlog .EQ. 1) then
               aa = rr(0)
               bb = rupper0 / length0 ! R_sun => non-dim. length
               radorb2 = lograd(radorb2,aa,bb,irlog)
               if (iout .GT. 0)
     &           write(*,'('' Rorbit '',f10.3,''=>'',f10.3)')
     &             radorb, radorb2
             endif
             call orb2rgb(radorb2)
             if (iout .GT. 0) write(*,*) 'ORB2RGB ..end'
             write(11,*) 'ORB2RGB ..end'
           endif

* HCS
           if (iplthcs .EQ. 1) then
             ivar = 1
             call sclr2rgb(ivar)
             if (iout .GT. 0) write(*,*) 'HCS2RGB  ..end'
             write(11,*) 'HCS2RGB  ..end'
           else if (iplthcs .EQ. 2) then
             ivar = 1
             call scl2rgb2(ivar)
             if (iout .GT. 0) write(*,*) 'HCS2RGB shadowed ..end'
             write(11,*) 'HCS2RGB shadowed ..end'
           endif

* Field line
*           if (ipltmag .EQ. 1) then
           if ((ipltmag .GE. 1) .AND. (ipltmag .LE. 8))  then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
             endif
             irtube = 150 ! index of radius of tube ; larger irtube => slender tube
             ivar   = 1
             itubecol = ipltmag  + 10
             call tube2rgb(irtube,ivar,itubecol,v0,numframe)
             if (iout .GT. 0) write(*,*) 'TUBE2RGB 1 ..end'
             write(11,*) 'TUBE2RGB 1 ..end'
           endif

* Flow line
           if (((ipltflow .GE. 1) .AND. (ipltflow .LE. 4)) .OR.
     &         ((ipltflow .LE.-1) .AND. (ipltflow .GE.-4))) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
             endif
             ipltflow = abs(ipltflow)
*             irtube = 400 ! index of radius of tube ; larger irtube => slender tube
             irtube = 150 ! index of radius of tube ; larger irtube => slender tube
             ivar   = 2
             itubecol = ipltflow
             call tube2rgb(irtube,ivar,itubecol,v0,numframe)
             if (iout .GT. 0) write(*,*) 'TUBE2RGB 2 ..end'
             write(11,*) 'TUBE2RGB 2 ..end'
           endif

* rho / rho_ref.
           if (ipltrho .LE. 0) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 9
             nref = abs(ipltrho)
             rotref = nrot2
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) 'rho2RGB ..end'
             write(11,*) 'rho2RGB ..end'
* Sudden change of mass density
           else if (ipltrho .EQ. 1) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 2
             call sclr2rgb(ivar)
             if (iout .GT. 0) write(*,*) 'MASS2RGB ..end'
             write(11,*) 'MASS2RGB ..end'
* integ. corona pB
           else if (ipltrho .EQ. 2) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 5
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) 'pB2RGB ..end'
             write(11,*) 'INTG2RGB ..end'
* sudden change of V_r
           else if (ipltrho .EQ. 3) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 6
             call sclr2rgb(ivar)
             if (iout .GT. 0) write(*,*) 'D_VR2RGB ..end'
             write(11,*) 'D_VR2RGB ..end'
* intg rho/rho_park
           else if (ipltrho .EQ. 4) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 7
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0) ! Los of N/N_p > 0.85
*             call sclr2rgb(ivar) ! draw where rho/rho_p = 0.85,,,
             if (iout .GT. 0) write(*,*) 'N/Park ..end'
             write(11,*) 'N/Park ..end'
* Countour of dev(rho)_i
           else if (ipltrho .EQ. 12) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 17
             call sclr2rgb(ivar)
             if (iout .GT. 0) write(*,*) 'DevN_iR2RGB ..end'
             write(11,*) 'DevN_iR2RGB ..end'
* plot rho/rho_park
           else if (ipltrho .EQ. 4) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 7
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0) ! Los of N/N_p > 0.85 MIND
*             call sclr2rgb(ivar) ! draw where rho/rho_p = 0.85,,,
             if (iout .GT. 0) write(*,*) 'N/Park ..end'
             write(11,*) 'N/Park ..end'
* intg Pg/Pg_park
           else if (ipltrho .EQ. 11) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 16
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
*             call sclr2rgb(ivar) ! draw where rho/rho_p = 0.85,,,
             if (iout .GT. 0) write(*,*) 'P/Park ..end'
             write(11,*) 'P/Park ..end'
* integ. IPS
           else if (ipltrho .EQ. 5) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 8
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) 'IPS2RGB ..end'
             write(11,*) 'IPS2RGB ..end'
* integ. Density
           else if (ipltrho .EQ. 6) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 10
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) 'DenR^2 2 RGB ..end'
             write(11,*) 'DenR^2 2 RGB ..end'
* integ. 1 / V_r
           else if (ipltrho .EQ. 7) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
                call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 11
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) ' 1/Vr 2 RGB ..end'
             write(11,*) '1/Vr 2 RGB ..end'
* integ. Std(rho)_i
           else if (ipltrho .EQ. 8) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 12
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) 'Std(N)_i 2RGB ..end'
             write(11,*) 'Std(rho)_i 2 RGB ..end'
* integ. Std(rho/v)_i
           else if (ipltrho .EQ. 9) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 13
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) 'Std(N/V)_i ..end'
             write(11,*) 'Std(N/V)_i 2 RGB ..end'
* integ. rho.
           else if (ipltrho .EQ. 10) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0
               call resetrad(rr, rbottom, rupper, iout)
             endif
             ivar = 14
             call intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
             if (iout .GT. 0) write(*,*) ' N2RGB ..end'
             write(11,*) 'N2RGB ..end'
           endif

* particle movement
           if (iprtcl .EQ. 1) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             call drawprt(numframe,v0)
             if (iout .GT. 0) write(*,*) 'Particle ..end'
             write(11,*) 'Particle ..end'
           endif

* IPS obs. LoS
           if ((iips .EQ. 1) .OR. (iips .EQ. 2)) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             call ipso2rgb(nrot,v0,length0,iips)
             if (iout .GT. 0) write(*,*) 'IPSo2RGB ..end'
             write(11,*) 'IPSo2RGB ..end'
           endif

* plot grid system
           if (igrid .EQ. 1) then
             if (irlog .EQ. 1) then
               call readinit(gamma,t0,length0,v0,omega)
               rbottom = rr(0)
               rupper  = rupper0 / length0 ! R_sun => non-dim. length
               call resetrad(rr, rbottom, rupper, iout)
             endif
             call grid2rgb(is, ie, js, je, ks, ke)
             if (iout .GT. 0) write(*,*) 'GRID2RGB ..end'
             write(11,*) 'GRID2RGB ..end'
           endif

         else
           if (iout .GT. 0) write(*,*) 'Data file not found'
           write(11,*) 'Data file not found'
         endif

* make rgb data array small
         call rgbsmall(rgb2,rgb0)
         if (iout .GT. 0) write(*,*) 'RGBsmall .. done'
         write(11,*) 'RGBsmall .. done'

* re-def. color (not always to be done)
         if (ihist .EQ. 1) then
           call reclrtbl(rgb0)
           if (iout .GT. 0) write(*,*) 'Histgram re-color .. done'
           write(11,*) 'Recolor .. done'
         endif

         do 501 j0 = 1, jvpix0
         do 501 i0 = 1, ihpix0
         do 501 k0 = 1, 3
           rgbx(k0,i0,j0) = rgb0(k0,i0,j0)
 501     continue

* write & prepare character(s)
         itrape = 0 ! (2/1/0) = not_draw/transparent/non-transparent

         ix = ihpix0 - 100
         jy = jvpix0 - 15
         open(12,file='date.txt',status='old')
         read(12,*) charout
         close(12)
         charlen = 9
         call charwri(charout,charlen,ix,jy,itrape,rgbx)
**
*         ix = 1
*         jy = jvpix0 - 15
*         charout='2 days before'
*         charlen = 13
*         call charwri(charout,charlen,ix,jy,itrape,rgbx)
**

         ix = 10
         jy = jvpix0 - 20
         idummy1 = nrot + 10000
         idummy1 = mod(idummy1,10000)
         idummy2 = int(seeph0)
         idummy3 = int(seeth0)
         write(charout,'(''CR'',i4.4,x,i3.3,x,i3.3)')
     &     idummy1,idummy2,idummy3
         charlen = 14
         if (iwristr .EQ. 1)
     &     call charwri(charout,charlen,ix,jy,itrape,rgbx)

         jy = jy - 15
*         write(charout,'(''Step '',i6.6)') nnn
         write(charout,'(''Time '',i6.6)') nnn
         charlen = 11
         if (iwristr .EQ. 1)
     &     call charwri(charout,charlen,ix,jy,itrape,rgbx)

* for ps
*         write(sdummy1,'(''CR'',i4,1x,''Step'',i7)') idummy1, nnn
         write(sdummy1,'(''CR'',i4,1x,''Time'',i7)') idummy1, nnn
         numpstr = numpstr + 1
         psstr(numpstr) = sdummy1
*         psstr(numpstr) = sdummy1//sdummy2 ! here sdummy2 is not defined ....
         if (iout .GT. 0) write(*,*) 'CHRWRI etc .. done'
         write(11,*) 'CHRWRI etc .. done'

         idummy1 = (numframe-1) / jframe
         idummy1 = mod(idummy1, iframe)
         jdummy1 = jframe - mod(numframe-1,jframe) - 1
         ixfp0 = idummy1*ihpix0
         jyfp0 = jdummy1*jvpix0

         if (iwristr .EQ. -1) then
           psix(numpstr) = ixfp0
           psjy(numpstr) = jyfp0 + jvpix0 - 15
         else
           psix(numpstr) = -7777
           psjy(numpstr) = -7777
         endif

* move each frame RGB to large frame area.
         do 500 j0 = 1, jvpix0
         do 500 i0 = 1, ihpix0
         do 500 k0 = 1, 3
           rgblarge(k0,ixfp0+i0,jyfp0+j0) = rgbx(k0,i0,j0)
 500     continue

         if (iout .GT. 0) write(*,*) 'Connect Image .. done'
         write(11,*) 'Connect Image .. done'

         onepage = onepage + 1

* out put
         if (onepage .EQ. maxfatp) then
           numpage = numpage + 1

           call doout(rgblarge,rgblargp,numpage,
     &                nrot,psstr,psix,psjy,numpstr,iout)

           onepage = 0
           numpstr = 0
         endif
       enddo

       if (onepage .GT. 0) then
         numpage = numpage + 1
         call doout(rgblarge,rgblargp,numpage,
     &              nrot,psstr,psix,psjy,numpstr,iout)
       endif

       close(11) ! close log file

       stop
       end



*
***=======================================================================
*
       subroutine drawprt(numframe,v0)
       implicit none
* interface
       integer numframe
       real*8  v0
* local
       integer l, m, n, i, j, k, idummy1, idummy2, idummy3
       real*8  rad, th, ph, aa
       real*8  rtpb(3), xyzb(3), vecb(3)
       real*8  rdummy1, rdummy2
       real*8  dstnt2, xp1, yp1, dd
       integer ix1,jy1,ic,jc,ix2,jy2
       integer red,gre,blu,red2,gre2,blu2
       real*8  dev, ave
* const
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* common
       integer prtmax
       parameter(prtmax = 200)
*       parameter(prtmax = 600)
       real*8  prtxyz(3,prtmax),prtspd(prtmax)
       logical ldrwprt(prtmax)
       integer numlat(11) /3,5,7,10,15,20,15,10,7,5,3/
*       integer numlat(23) /2,3,4,5,6,7,8,10,13,15,17,20,
*     +                     17,15,13,10,8,7,6,5,4,3,2/
       integer iout, irlog
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       common  /cgcfg/  iout, irlog
       common  /coord/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp
       common  /var12/  ur1, ut1, up1
       common  /images/ dstnt, rgb2
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom
       common  /prtprt/ prtxyz
       common  /lprtprt/ ldrwprt

       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

       do 400 i = 0, ii + 1
         ave = 0.0D+00
         dev = 0.0D+00
         do 401 k = 1, kk
         do 401 j = 0, jj
           aa = ur1(i,j,k)**2 + up1(i,j,k)**2 + up1(i,j,k)**2
           ave = ave + dsqrt(aa)
           dev = dev + aa
 401     continue
         ave = ave / dfloat(jj+1) / dfloat(kk)
         dev = dev / dfloat(jj+1) / dfloat(kk)
         dev = dev - ave**2
         dev = dabs(dev)
         dev = dsqrt(dev)
         dev = 1.0D+00 ! do nothing !
         do 402 k = 0, kk + 1
         do 402 j = 0, jj
*           vcr(i,j,k) = (ur1(i,j,k) - ave) / dev
*           vct(i,j,k) = (ut1(i,j,k) - ave) / dev
*           vcp(i,j,k) = (up1(i,j,k) - ave) / dev
           vcr(i,j,k) = ur1(i,j,k) / dev
           vct(i,j,k) = ut1(i,j,k) / dev
           vcp(i,j,k) = up1(i,j,k) / dev
 402     continue
 400   continue

* loc.of.particles
       if (numframe .EQ. 1) then ! set initial location
*         rad = rbottom * 1.1D+00
         n = 0
         do i = 1, 2 ! 3
*           if (i .EQ. 1) then
*             rad = rbottom * 1.1D+00
*           else
*             rad =(rbottom + rupper) / 2.0D+00
*           endif
           rad = rbottom * 1.1D+00
     &         + (rupper-rbottom*1.1D+00)*dfloat(i-1)/3.0D+00
           do l = 1, 11 ! 23
             th = dfloat(l) / 24.0D+00 * pi
             do m = 1, numlat(l)
               n = n + 1
               ph = dfloat(m) / dfloat(numlat(l)) * pi * 2.0D+00
               prtxyz(1,n) = rad * dsin(th) * dcos(ph)
               prtxyz(2,n) = rad * dsin(th) * dsin(ph)
               prtxyz(3,n) = rad * dcos(th)
               ldrwprt(n) = .true.
             enddo
           enddo
         enddo
       else
         do n = 1, prtmax
           do i = 1, 3
             xyzb(i) = prtxyz(i,n)
           enddo
           call xyz2rtp(xyzb,rtpb)
           call vecatrtp(vecb,rtpb)
           call vrtp2xyz(vecb,rtpb)
           aa = vecb(1)**2 + vecb(2)**2 + vecb(3)**2
           aa = dsqrt(aa)
           aa = (aa/3.0D+00)**2 ! <===== MIND : speed is emphasized !!!
     &        /((rupper-rbottom)/2.0D+02)
*           aa = dfloat(15 * 60) / 30000.0D+00
*     +        /((rupper-rbottom)/2.0D+02)
*           aa = 1.0D+01 * aa
           do i = 1, 3
             prtxyz(i,n) = xyzb(i) + aa * vecb(i)
           enddo
           rad = prtxyz(1,n)**2 + prtxyz(2,n)**2 + prtxyz(3,n)**2
           rad = dsqrt(rad)
           if ((rad .GT. rupper) .OR. (rad .LT. rbottom)) then
             do i = 1, 3
               prtxyz(i,n) = prtxyz(i,n) * rbottom / rad
             enddo
           endif
* if needed, set conditions of drawing or not
*           ldrwprt(n) = .false. , ............
         enddo
       endif

* flow speed
       do n = 1, prtmax
         do i = 1, 3
           xyzb(i) = prtxyz(i,n)
         enddo
         call xyz2rtp(xyzb,rtpb)
         call vecatrtp(vecb,rtpb)
         prtspd(n) = vecb(1) * v0 / 1.0D+05
       enddo

* draw particle
       do n = 1, prtmax
        if (ldrwprt(n)) then
         do i = 1, 3
           xyzb(i) = prtxyz(i,n)
         enddo
         call prjpnt1(xyzb,xp1,yp1,dstnt2)
         call pxy2ij(xp1,yp1,ul,ix2,jy2,ic,jc)
         rdummy2 = prtspd(n)
         rdummy1 = 1.0D+00
         call ipscolor(red2,gre2,blu2,rdummy2,rdummy1)

         dd = 3.0D+00 / (dstnt2/rupper)
         do 500 ix1 = 1, ihpix2
         do 500 jy1 = 1, jvpix2
           idummy1 = (ix1 - ix2)**2 + (jy1 - jy2)**2
           rdummy1 = - dfloat(idummy1) / dd**2
           rdummy1 = dexp(rdummy1)
           rdummy2 = dfloat(red2) * rdummy1
           red = int(rdummy2)
           rdummy2 = dfloat(gre2) * rdummy1
           gre = int(rdummy2)
           rdummy2 = dfloat(blu2) * rdummy1
           blu = int(rdummy2)

           if (dstnt2 .LE. dstnt(ix1,jy1)) then
             idummy1 = ichar(rgb2(1,ix1,jy1)) + red
             idummy2 = ichar(rgb2(2,ix1,jy1)) + gre
             idummy3 = ichar(rgb2(3,ix1,jy1)) + blu
             if (idummy1 .GT. 255) idummy1 = 255
             if (idummy2 .GT. 255) idummy2 = 255
             if (idummy3 .GT. 255) idummy3 = 255
             if (idummy1 .LT.   0) idummy1 = 0
             if (idummy2 .LT.   0) idummy2 = 0
             if (idummy3 .LT.   0) idummy3 = 0
             rgb2(1,ix1,jy1) = char(idummy1)
             rgb2(2,ix1,jy1) = char(idummy2)
             rgb2(3,ix1,jy1) = char(idummy3)
           endif
 500     continue
        endif
       enddo

       return
       end

*
***=======================================================================
*
       subroutine doout(rgblarge,rgblargp,numpage,
     &                  nrot,psstr,psix,psjy,numpstr,iout)
       implicit none
* interface : RGB data storage
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpixf, jvpixf, iframe, jframe, maxfatp
       parameter(iframe = 1, jframe = 1)
       parameter(maxfatp = iframe * jframe)
       parameter(ihpixf = ihpix0*iframe, jvpixf = jvpix0*jframe)
       character*1 rgblarge(3,ihpixf,jvpixf) ! RGB data for present large frame
       character*1 rgblargp(3,ihpixf,jvpixf) ! RGB data for previous large frame
       character*1 rgblargx(3,ihpixf,jvpixf) ! RGB data for large frame output
       integer numpage
       integer nrot
       character*33 psstr(maxfatp)
       integer psix(maxfatp), psjy(maxfatp), numpstr
       integer iout
* local
       integer ismthe
       parameter(ismthe = 1)  ! index for fade in and out
       integer ismth, ismthend
       real*8  aa, bb, rdummy1
       integer idummy1, idummy2, idummy0
       integer i0, j0, k0

* smoothig frames or not.
       if (numpage .EQ. 1) then
         ismthend = 1
         do 401 j0 = 1, jvpixf
         do 401 i0 = 1, ihpixf
         do 401 k0 = 1, 3
           rgblargp(k0,i0,j0) = char(0) ! initialize
 401     continue
       else
         ismthend = ismthe
       endif

* write out rgb data to a file or standard output device
       do ismth = 1, ismthend
         aa =           dfloat(ismth) / dfloat(ismthend)
         bb = 1.0D+00 - dfloat(ismth) / dfloat(ismthend)
         do 400 j0 = 1, jvpixf
         do 400 i0 = 1, ihpixf
         do 400 k0 = 1, 3
           idummy1 = ichar(rgblargp(k0,i0,j0))
           idummy2 = ichar(rgblarge(k0,i0,j0))
           rdummy1 = dfloat(idummy1) * bb + dfloat(idummy2) * aa
           idummy0 = int(rdummy1)
           if (idummy0 .GT. 255) idummy0 = 255
           if (idummy0 .LT.   0) idummy0 =   0
           rgblargx(k0,i0,j0) = char(idummy0)
 400     continue
         if (iout .EQ. 9) then
           call bmpout(iout,rgblargx,nrot)
         else if (iout .LE. 5) then
           call rgbout(iout,rgblargx,nrot)
         else
           call psout(iout,rgblargx,nrot,
     &                psstr,psix,psjy,numpstr)
         endif
         if (iout .GT. 0) write(*,*) 'done .. IMAGEOUT'
       enddo

* copy rgb data
       do 350 j0 = 1, jvpixf
       do 350 i0 = 1, ihpixf
       do 350 k0 = 1, 3
         rgblargp(k0,i0,j0) = rgblarge(k0,i0,j0)
 350   continue
       if (iout .GT. 0) write(*,*) 'Copy RGB .. done'
       write(11,*) 'Copy RGB .. done'

       return
       end


*
***=======================================================================
*
       subroutine orb2rgb(radorb)
       implicit none
       real*8  radorb
* local
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       integer numcirc
       parameter(numcirc = 100)
       real*8  aa, xyzs(3), xyze(3)
       integer ic, jc, m
       integer red, gre, blu, itransp
* common
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       common  /images/ dstnt, rgb2
       common  /coord/  rr, theta, phi
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

       do m = 1, numcirc
         aa = 2.0D+00 * pi * dfloat(m-1) / dfloat(numcirc)
         xyzs(1) = radorb * dcos(aa)
         xyzs(2) = radorb * dsin(aa)
         xyzs(3) = 0.0D+00
         aa = 2.0D+00 * pi * dfloat(m  ) / dfloat(numcirc)
         xyze(1) = radorb * dcos(aa)
         xyze(2) = radorb * dsin(aa)
         xyze(3) = 0.0D+00

         red = 255
         gre = 255
         blu = 255
         itransp = 0 ! Non-transparent
         call drawline(xyzs,xyze,red,gre,blu,ic,jc,itransp)
       enddo

       return
       end


*
***=======================================================================
*
       subroutine reclrtbl(rgb0) ! reset color table
       implicit none
* interface
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256,          jvpix0 = 256)
       character*1 rgb0(3,ihpix0, jvpix0) ! RGB data for output
* local variables
       integer i0, j0, k, l
       integer idummy1, ilumirgb
       real*8  rdummy1, rdummy2, rdummy3, rlumirgb
       integer numsum(0:255), nsum2, numpix
       real*8  rintnum(1:254)

* stat.
       do k = 0, 255
         numsum(k) = 0
       enddo
       do 400 j0 = 1, jvpix0
       do 400 i0 = 1, ihpix0
         idummy1 = ichar(rgb0(1,i0,j0))
         rdummy1 = dfloat(idummy1) * 0.299
         idummy1 = ichar(rgb0(2,i0,j0))
         rdummy1 = dfloat(idummy1) * 0.587 + rdummy1
         idummy1 = ichar(rgb0(3,i0,j0))
         rdummy1 = dfloat(idummy1) * 0.114 + rdummy1
         idummy1 = int(rdummy1)
         numsum(idummy1) = numsum(idummy1) + 1
 400   continue

       numpix = 0
       do l = 1, 254
         numpix = numpix + numsum(l)
       enddo

       if (numpix .GT. 0) then
         do k = 1, 254
           nsum2 = 0
           do l = 1, k
             nsum2 = nsum2 + numsum(l)
           enddo
           rintnum(k) = dfloat(nsum2) / dfloat(numpix)
         enddo

         do 410 j0 = 1, jvpix0
         do 410 i0 = 1, ihpix0
           idummy1  = ichar(rgb0(1,i0,j0))
           rlumirgb = dfloat(idummy1) * 0.299
           idummy1  = ichar(rgb0(2,i0,j0))
           rlumirgb = dfloat(idummy1) * 0.587 + rlumirgb
           idummy1  = ichar(rgb0(3,i0,j0))
           rlumirgb = dfloat(idummy1) * 0.114 + rlumirgb
           ilumirgb = int(rlumirgb)

           if ((ilumirgb .NE. 0) .AND. (ilumirgb .NE. 255)) then
             rdummy2 = rintnum(ilumirgb) * 255.00
             rdummy3 = rdummy2 / rlumirgb
             do k = 1, 3
               idummy1 = ichar(rgb0(k,i0,j0))
               rdummy1 = dfloat(idummy1) * rdummy3
               idummy1 = int(rdummy1)
               if (idummy1 .LT. 0)   idummy1 = 0
               if (idummy1 .GT. 255) idummy1 = 255
               rgb0(k,i0,j0) = char(idummy1)
             enddo
           endif
 410     continue
       endif

       return
       end


*
***=======================================================================
*

       subroutine setgrap(xyzcp,dproj,xyzxp,xyzyp,xyzlght,
     &                    rupper,rbottom,seeth,seeph,
     &                    ul,numframe,numend,leyearth,rau)
       implicit none
* interface
       real*8  xyzcp(3)
       real*8  dproj
       real*8  xyzxp(3), xyzyp(3)
       real*8  xyzlght(3)
       real*8  rupper, rbottom, seeth, seeph
       real*8  ul
       integer numframe, numend
       logical leyearth
       real*8  rau
* dummy
       real*8  dummyxyz(3)
       real*8  rdummy1
       integer idummy
       real*8  absvec3
* local
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
* constant
       real*8  pi
       parameter(pi = 3.14159265358979D+00)

* surpress compiler's warnings
       idummy = numend
       idummy = numframe
       rdummy1 = rbottom

       seeth = seeth / 180.0 * pi
       seeph = seeph / 180.0 * pi

* distance between the origin of projection plane and point of view
       if (leyearth) then
         if (rupper .GT. 0.9 * rau) rupper = 0.9 * rau
         dproj    = rau - rupper ! assuming that 3 points are on the same line
       else
*         dproj    = 1.0D+01 * rupper ! standard for one-point view
         dproj    = 1.0D+03 * rupper ! almost parallel view.
       endif

* location of center (origin) of projection plane
       xyzcp(1) = 1.01D+00 * rupper * dsin(seeth) * dcos(seeph)
       xyzcp(2) = 1.01D+00 * rupper * dsin(seeth) * dsin(seeph)
       xyzcp(3) = 1.01D+00 * rupper * dcos(seeth)

* unit vectors on projection plane (perpendicular each other, normalized)
       dummyxyz(1) = 0.0D+00
       dummyxyz(2) = 0.0D+00
       dummyxyz(3) = 1.0D+00
       call gaiseki3(dummyxyz,xyzcp,xyzxp)

       rdummy1 = absvec3(xyzxp)
       rdummy1 = dabs(rdummy1)
       if (rdummy1 .LT. 1.0D-01) then
         dummyxyz(1) = 1.0D+00
         call gaiseki3(dummyxyz,xyzcp,xyzxp)
       endif

       call normvec3(xyzxp,rdummy1)
       call gaiseki3(xyzcp,xyzxp,xyzyp)
       call normvec3(xyzyp,rdummy1)
* location of light source
*       xyzlght(1) =-xyzcp(1)+(xyzxp(1)*0.5D+00+xyzyp(1))*rupper
*       xyzlght(2) =-xyzcp(2)+(xyzxp(2)*0.5D+00+xyzyp(2))*rupper
*       xyzlght(3) =-xyzcp(3)+(xyzxp(3)*0.5D+00+xyzyp(3))*rupper

*       do m = 1, 3
*         xyzlght(m) = xyzcp(m)*0.5D+00
*     +              +(xyzxp(m)*0.5D+00+xyzyp(m))*rupper
*       enddo

*       xyzlght(1) = 0.0D+00
*       xyzlght(2) = 0.0D+00
*       xyzlght(3) = rupper*3.0D+00

       xyzlght(1) = 3.0D+00 * rupper*dsin(seeth)*dcos(seeph-pi*0.33) ! - 60 degr.
       xyzlght(2) = 3.0D+00 * rupper*dsin(seeth)*dsin(seeph-pi*0.33)
       xyzlght(3) = 3.0D+00 * rupper*dcos(seeth)

* set magnification factor : unit length for one pixel size
       rdummy1 = (rupper + dproj)**2 - rupper**2
       rdummy1 = dsqrt(rdummy1)
       rdummy1 = dproj * rupper / rdummy1
       ul = 2.0D+00 * rdummy1 * 1.15D+00 / dfloat(jvpix2)

       return
       end


*
***=======================================================================
*
       subroutine ipso2rgb(rt,v0,length0,iips)
       implicit none
* interface
       integer rt
       real*8  v0, length0
       integer iips ! 1 = white, 2 = color
* local
       integer i, m
       integer ic, jc          ! pixel potition of center of project plane
       integer red, gre, blu
       integer itransp
* constant
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* dummy
       logical lipsfile
       real*8  aa, bb, cc
* IPS obs.
       integer maxobs
       parameter(maxobs = 1000)
       real*8  th1o(maxobs), th2o(maxobs)
       real*8  ph1o(maxobs), ph2o(maxobs)
       real*8  radp(maxobs), thep(maxobs), phip(maxobs)
       real*8  valips(maxobs), errips(maxobs)
       integer numobs
*
       integer n
       integer nlineseg
       parameter(nlineseg = 10)
       real*8  xyzs(3), xyze(3), xyza(3), xyzb(3), xyzc(3)
       real*8  th1, th2, ph1, ph2
       real*8  ra1au
       real*8  spdips, dpth
* common
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       common  /images/ dstnt, rgb2
       common  /coord/  rr, theta, phi
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

* pixel position on the projection plane
       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

* avoid compiler warning
       aa = v0

       ra1au = 215.0 / length0
*       write(*,*) ' ra1au = ' , ra1au

       call readips(rt,th1o,th2o,ph1o,ph2o,
     &              radp,thep,phip,valips,errips,lipsfile,numobs)

       if (lipsfile) then

         do i = 1, numobs
           th1o(i) = th1o(i) / 180.0D+00 * pi
           ph1o(i) = ph1o(i) / 180.0D+00 * pi
           th2o(i) = th2o(i) / 180.0D+00 * pi
           ph2o(i) = ph2o(i) / 180.0D+00 * pi
           valips(i) = valips(i) ! * 1.0D+05 / v0
           errips(i) = errips(i) ! * 1.0D+05 / v0
         enddo

         do i = 1, numobs
           th1 = th1o(i)
           ph1 = ph1o(i)
           th2 = th2o(i)
           ph2 = ph2o(i)
           xyzs(1) = ra1au * dsin(th1) * dcos(ph1)
           xyzs(2) = ra1au * dsin(th1) * dsin(ph1)
           xyzs(3) = ra1au * dcos(th1)
           xyze(1) = ra1au * dsin(th2) * dcos(ph2)
           xyze(2) = ra1au * dsin(th2) * dsin(ph2)
           xyze(3) = ra1au * dcos(th2)
           spdips  = valips(i)

           if (iips .EQ. 1) then
             red = 255
             gre = 255
             blu = 255
             itransp = 0
             call drawline(xyzs,xyze,red,gre,blu,ic,jc,itransp)
           endif
           if (iips .EQ. 2) then
             do n = 1, nlineseg
               aa = dfloat(n)            / dfloat(nlineseg)
               bb = dfloat(n-1)          / dfloat(nlineseg)
               cc =(dfloat(n) - 0.5D+00) / dfloat(nlineseg)
               do m = 1, 3
                 xyza(m) =  xyze(m) * aa + xyzs(m) * (1.0D+00 - aa)
                 xyzb(m) =  xyze(m) * bb + xyzs(m) * (1.0D+00 - bb)
                 xyzc(m) = (xyze(m) + xyzs(m)) / 2.0
               enddo
               dpth = cc - 0.5
               dpth = - dpth**2 * 20.0
               dpth = dexp(dpth)
               call ipscolor(red,gre,blu,spdips,dpth)
               itransp = 1
               call drawline(xyzb,xyza,red,gre,blu,ic,jc,itransp)
             enddo
           endif
         enddo
       endif

       return
       end

*
** ----------------------
*
       subroutine velcolor(ired,igre,iblu,vel,mx) ! black-red-white
       implicit none
*
       integer ired,igre,iblu
       real*8  vel, mx
*
       real*8  red,gre,blu
       real*8  aa
**
       aa = vel / mx
       aa = 1.0 - aa ! upside-down...

       if (aa .LT. 0.5) then
         aa = aa * 2.0
         gre = aa
         red = 0.0
         blu = 0.0
       else
         aa = (aa - 0.5) * 2.0
         gre = 1.0
         red = aa
         blu = aa
       endif

* IPS-vel. like
*       aa = -stdtmp + 1.0 ! red is hotter ..
*       if (aa .LT. -0.50) then
*         red = (aa + 0.25) / 0.5
*         gre = 0.0
*         blu = 0.0
*       else if (aa .LT. -0.25) then
*         red = 1.0
*         gre = (aa + 0.125) / 0.125
*         blu = 0.0
*       else if (aa .LT.  0.00) then
*         red = (0.0 - aa) / 0.25
*         gre = 1.0
*         blu = 0.0
*       else if (aa .LT.  0.25) then
*         red = 0.0
*         gre = 1.0
*         blu = (aa - 0.0) / 0.25
*       else if (aa .LT.  0.50) then
*         red = 0.0
*         gre = (0.5 - aa) / 0.25
*         blu = 1.0
*       else
*         red = 0.0
*         gre = 0.0
*         blu = (1.0 - aa) / 0.50
*       endif

       if (red .GT. 0.99) red = 0.99
       if (red .LT. 0.01) red = 0.01
       if (gre .GT. 0.99) gre = 0.99
       if (gre .LT. 0.01) gre = 0.01
       if (blu .GT. 0.99) blu = 0.99
       if (blu .LT. 0.01) blu = 0.01

       ired = int(red * 255.0)
       igre = int(gre * 255.0)
       iblu = int(blu * 255.0)

       return
       end
*
** ----------------------
*
       subroutine tmpcolor(ired,igre,iblu,stdtmp)
       implicit none
*
       integer ired,igre,iblu
       real*8  stdtmp
*
       real*8  red,gre,blu
       real*8  aa

**
*       aa = (stdtmp + 1.0D+00)
*     &    / 5.0        ! simga covered
**
       aa = stdtmp ! ........ if max/min

       if (aa .LT. 0.5) then
         aa = aa * 2.0
         red = aa ! red
         gre = 0.0
         blu = 0.0
       else
         aa = (aa - 0.5) * 2.0
         red = 1.0 !  red !!
         gre = aa
         blu = aa
       endif

* IPS-vel. like
*       aa = -stdtmp + 1.0 ! red is hotter ..
*       if (aa .LT. -0.50) then
*         red = (aa + 0.25) / 0.5
*         gre = 0.0
*         blu = 0.0
*       else if (aa .LT. -0.25) then
*         red = 1.0
*         gre = (aa + 0.125) / 0.125
*         blu = 0.0
*       else if (aa .LT.  0.00) then
*         red = (0.0 - aa) / 0.25
*         gre = 1.0
*         blu = 0.0
*       else if (aa .LT.  0.25) then
*         red = 0.0
*         gre = 1.0
*         blu = (aa - 0.0) / 0.25
*       else if (aa .LT.  0.50) then
*         red = 0.0
*         gre = (0.5 - aa) / 0.25
*         blu = 1.0
*       else
*         red = 0.0
*         gre = 0.0
*         blu = (1.0 - aa) / 0.50
*       endif

       if (red .GT. 0.99) red = 0.99
       if (red .LT. 0.01) red = 0.01
       if (gre .GT. 0.99) gre = 0.99
       if (gre .LT. 0.01) gre = 0.01
       if (blu .GT. 0.99) blu = 0.99
       if (blu .LT. 0.01) blu = 0.01

       ired = int(red * 255.0)
       igre = int(gre * 255.0)
       iblu = int(blu * 255.0)

       return
       end

*
** ----------------------
*
       subroutine ipscolor(ired,igre,iblu,ipsvel,depth)
       implicit none
       integer ired,igre,iblu
       real*8  red,gre,blu
       real*8  ipsvel, depth

       if (ipsvel .LT. 250.0) then
         red = 0.5
         gre = 0.0
         blu = 0.0
       else if (ipsvel .LT. 300.0) then
         red = (ipsvel - 200.0) / 100.0
         gre = 0.0
         blu = 0.0
       else if (ipsvel .LT. 400.0) then
         red = 1.0
         gre = (ipsvel - 300.0) / 100.0
         blu = 0.0
       else if (ipsvel .LT. 500.0) then
         red = (500.0 - ipsvel) / 100.0
         gre = 1.0
         blu = 0.0
       else if (ipsvel .LT. 600.0) then
         red = 0.0
         gre = 1.0
         blu = (ipsvel - 500.0) / 100.0
       else if (ipsvel .LT. 700.0) then
         red = 0.0
         gre = (700.0 - ipsvel) / 100.0
         blu = 1.0
       else if (ipsvel .LT. 800.0) then
         red = 0.0
         gre = 0.0
         blu = (900.00 - ipsvel) / 200.0
       else
         red = 0.0
         gre = 0.0
         blu = 0.5
       endif
       red = red * depth
       gre = gre * depth
       blu = blu * depth
       if (red .GT. 0.99) red = 0.99
       if (red .LT. 0.01) red = 0.01
       if (gre .GT. 0.99) gre = 0.99
       if (gre .LT. 0.01) gre = 0.01
       if (blu .GT. 0.99) blu = 0.99
       if (blu .LT. 0.01) blu = 0.01

       ired = int(red * 255.0)
       igre = int(gre * 255.0)
       iblu = int(blu * 255.0)

       return
       end

*
** ----------------------
*
       subroutine ipscolr2(ired,igre,iblu,ipsvel,depth)
       implicit none
       integer ired,igre,iblu
       real*8  red,gre,blu
       real*8  ipsvel, depth
*
       real*8  maxlocal
       parameter(maxlocal = 100.0D+00)
       real*8  velocal

       velocal = ipsvel / maxlocal * 500.0 + 250.0 ! need adjustment each time if used

       if (velocal .LT. 250.0) then
         red = 0.5
         gre = 0.0
         blu = 0.0
       else if (velocal .LT. 300.0) then
         red = (velocal - 200.0) / 100.0
         gre = 0.0
         blu = 0.0
       else if (velocal .LT. 400.0) then
         red = 1.0
         gre = (velocal - 300.0) / 100.0
         blu = 0.0
       else if (velocal .LT. 500.0) then
         red = (500.0 - velocal) / 100.0
         gre = 1.0
         blu = 0.0
       else if (velocal .LT. 600.0) then
         red = 0.0
         gre = 1.0
         blu = (velocal - 500.0) / 100.0
       else if (velocal .LT. 700.0) then
         red = 0.0
         gre = (700.0 - velocal) / 100.0
         blu = 1.0
       else if (velocal .LT. 800.0) then
         red = 0.0
         gre = 0.0
         blu = (900.00 - velocal) / 200.0
       else
         red = 0.0
         gre = 0.0
         blu = 0.5
       endif
       red = red * depth
       gre = gre * depth
       blu = blu * depth
       if (red .GT. 0.99) red = 0.99
       if (red .LT. 0.01) red = 0.01
       if (gre .GT. 0.99) gre = 0.99
       if (gre .LT. 0.01) gre = 0.01
       if (blu .GT. 0.99) blu = 0.99
       if (blu .LT. 0.01) blu = 0.01

       ired = int(red * 255.0)
       igre = int(gre * 255.0)
       iblu = int(blu * 255.0)

       return
       end



*
** ---------------
*
       subroutine readips(rt,th1o,th2o,ph1o,ph2o,
     &                       radp,thep,phip,valips,errips,
     &                       lipsfile,numobs)
       implicit none
* interface
       integer rt
       integer maxobs
       parameter(maxobs = 1000)
       real*8  th1o(maxobs), th2o(maxobs)
       real*8  ph1o(maxobs), ph2o(maxobs)
       real*8  radp(maxobs), thep(maxobs), phip(maxobs)
       real*8  valips(maxobs), errips(maxobs)
       logical lipsfile
       integer numobs
* local
       character*40 fipsdata
* reading interface dummy
       real*8  thdummy1, thdummy2
       real*8  phdummy1, phdummy2
       real*8  rdummy(3), aa
       integer ivv, ivverr, i
* lower bottom of IPS data
       real*8  rabot
       parameter(rabot = 0.3D+00)

       write(fipsdata,'(''ipsv'',i4.4,''.dat'')') rt
       inquire(exist=lipsfile,file=fipsdata)
       if (lipsfile) then
         open(unit=14,file=fipsdata,status='old')
         numobs = 0
 777     continue
           if (numobs .LE. maxobs) then
             read(14,*,ERR=999,END=999)
     &         thdummy1, phdummy1,
     &         thdummy2, phdummy2,
     &         (rdummy(i),i=1,3), ivv, ivverr
             aa = dfloat(ivverr) / dfloat(ivv)
             if ((rdummy(1) .GT. rabot) .AND.
     &           (rdummy(1) .LT. 1.0D+00) .AND.
     &            (aa .LT. 0.2D+00)) then
               numobs = numobs + 1
               valips(numobs) = dfloat(ivv)
               errips(numobs) = dfloat(ivverr)
               th1o(numobs) = thdummy1
               ph1o(numobs) = phdummy1
               th2o(numobs) = thdummy2
               ph2o(numobs) = phdummy2
               radp(numobs) = rdummy(1) ! radii at P-point
               thep(numobs) = rdummy(2)
               phip(numobs) = rdummy(3)
             endif
             goto 777
           endif
 999     continue
         close(14)
         write(11,*) 'numobs = ', numobs
       else
         write(11,*) ' File not found : ', fipsdata
         numobs = 0
       endif

       return
       end


*
***=======================================================================
*
       subroutine grid2rgb(is, ie, js, je, ks, ke)
       implicit none
* interface
       integer is, ie, js, je, ks, ke
* local
*       integer idummy1, idummy2, idummy3
       integer i, j, k, k2, l
       integer ic, jc          ! pixel potition of center of project plane
       real*8  xp1, yp1, xp2, yp2
       integer red, gre, blu
       integer ijk(3)
       real*8  rtp(3), xyz1(3), xyz2(3)
       real*8  dstnt2
       integer itransp
* constant
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* dummy
       logical ldummy
* common
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       common  /images/ dstnt, rgb2
       common  /coord/  rr, theta, phi
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

* pixel position on the projection plane
       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

* check the range
       rtp(1) = rupper
       rtp(2) = pi / 4.0D+00
       rtp(3) = pi / 4.0D+00
       call rtp2ijk(rtp,ijk)
       if (ie .GT. ijk(1)) ie = ijk(1)
       if (is .GT. ie) is = ie - 3

       red = 240
       gre = 240
       blu = 120

* drawing the mesh system
       do 100 i = is, ie
       do 100 j = js, je
       do 100 k = ks, ke
       do 100 l = 1, 3

         rtp(1) = rr(i)
         rtp(2) = theta(j)
         rtp(3) = phi(k)
         call rtp2xyz(rtp,xyz1)
         call prjpnt1(xyz1,xp1,yp1,dstnt2)

* segments
         ldummy = .true.
         if (l .EQ. 1) then
           rtp(1) = rr(i-1)
           if (i .EQ. is) ldummy = .false.
         endif
         if (l .EQ. 2) then
           rtp(2) = theta(j+1)
           if (j .EQ. je) ldummy = .false.
         endif
         if (l .EQ. 3) then
           k2 = mod(k+1,kk)
           rtp(3) = phi(k2)
           if (k .EQ. ke) ldummy = .false.
         endif

         if (ldummy) then
           call rtp2xyz(rtp,xyz2)
           call prjpnt1(xyz2,xp2,yp2,dstnt2)

           itransp = 0
           call drawline(xyz1,xyz2,red,gre,blu,ic,jc,itransp)
         endif

 100   continue

       return
       end


*
** ----------------------------------------------------------
*
       subroutine drawline(xyz1,xyz2,red,gre,blu,ic,jc,itransp)
       implicit none
* interface
       real*8   xyz1(3), xyz2(3)
       integer  red, gre, blu
       integer  ic, jc, itransp
* local
       real*8  aa, ee
       real*8  xp1, yp1, xp2, yp2, xp, yp
       real*8  dstnt2
       integer i, ie, m
       integer ix1, jy1
       real*8  nprj(3), xyzeye(3)
* dummy
       integer idummy1, idummy2, idummy3
       real*8  xyzdummy(3)
       real*8  rdummy
       logical ldummy
* function
       logical sameside
* common
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  dstnt(ihpix2, jvpix2)
       real*8  rupper, rbottom
       common  /images/ dstnt, rgb2
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,rdummy)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo

       call prjpnt1(xyz1,xp1,yp1,dstnt2)
       call prjpnt1(xyz2,xp2,yp2,dstnt2)

       aa = (dabs(xp2 - xp1) + dabs(yp2 - yp1)) / ul
       ie = int(aa * 4.0) + 1

       do i = 0, ie
         ee = dfloat(i) / dfloat(ie)
         xyzdummy(1) = xyz1(1) + ee * (xyz2(1) - xyz1(1))
         xyzdummy(2) = xyz1(2) + ee * (xyz2(2) - xyz1(2))
         xyzdummy(3) = xyz1(3) + ee * (xyz2(3) - xyz1(3))

         call prjpnt1(xyzdummy,xp,yp,dstnt2)
         call pxy2ij(xp,yp,ul,ix1,jy1,ic,jc)

         ldummy = sameside(xyzcp,nprj,xyzdummy,xyzeye)

         if (((ix1 * (ihpix2 - ix1) .GT. 0)) .AND.
     &       ((jy1 * (jvpix2 - jy1) .GT. 0))) then

           if (.NOT. ldummy) then
             if (dstnt2 .LE. dstnt(ix1,jy1)) then
               if (itransp .EQ. 1) then
                 idummy1 = ichar(rgb2(1,ix1,jy1)) + red
                 idummy2 = ichar(rgb2(2,ix1,jy1)) + gre
                 idummy3 = ichar(rgb2(3,ix1,jy1)) + blu
               else
                 dstnt(ix1,jy1) = dstnt2
                 idummy1 = red
                 idummy2 = gre
                 idummy3 = blu
               endif
             else
               if (itransp .EQ. 0) then ! mix color
                 idummy1 = (2 * ichar(rgb2(1,ix1,jy1)) + red) / 3
                 idummy2 = (2 * ichar(rgb2(2,ix1,jy1)) + gre) / 3
                 idummy3 = (2 * ichar(rgb2(3,ix1,jy1)) + blu) / 3
               else ! do nothing
                 idummy1 = ichar(rgb2(1,ix1,jy1))
                 idummy2 = ichar(rgb2(2,ix1,jy1))
                 idummy3 = ichar(rgb2(3,ix1,jy1))
               endif
             endif
             if (idummy1 .GT. 255) idummy1 = 255
             if (idummy2 .GT. 255) idummy2 = 255
             if (idummy3 .GT. 255) idummy3 = 255
             rgb2(1,ix1,jy1) = char(idummy1)
             rgb2(2,ix1,jy1) = char(idummy2)
             rgb2(3,ix1,jy1) = char(idummy3)
           endif
         endif
       enddo

       return
       end


*
***=======================================================================
*
       subroutine resetrad(rr, rbottom, rupper, iout)
       implicit none
* interface
       integer ii
       parameter(ii = 72)
       real*8  rr(0:ii+1)
       real*8  rbottom, rupper
       integer iout
* local
       real*8  aa, ra
       integer i, irlog
       real*8  rbot2, rtop2
       parameter(rbot2 = 1.5D+00, rtop2 = 10.0D+00)
* functions
       real*8  lograd


       irlog = 1
       do i = 0, ii + 1
         ra = rr(i)
         aa = lograd(ra,rbottom,rupper,irlog)
         rr(i) = aa
       enddo

       rbottom = rbot2
       rupper  = rtop2

       if (iout .GT. 0) write(*,'(A)') ' R : log-scale .. done'
       write(11,'(A)') ' R : log-scale .. done'

       return
       end


*
** ---------------------
*
       real*8 function lograd(ra,rbottom,rupper,irlog)
       implicit none
* interface
       real*8  ra,rbottom,rupper
       integer irlog
* local
       real*8  rbot2, rtop2
       parameter(rbot2 = 1.5D+00, rtop2 = 10.0D+00)
       real*8  aa

       if (irlog .EQ. 1) then
         if (ra .LT. rbottom) then
           aa = ra / rbottom * rbot2
         else
           aa = dlog10(ra / rbottom) / dlog10(rupper / rbottom)
           aa = aa * (rtop2 - rbot2) + rbot2
         endif
       else
         aa = ra
       endif

       lograd = aa

       return
       end


*
***=======================================================================
*
       subroutine rgbsmall(rgb2,rgb0)
       implicit none
* interface
       integer ihpix0, jvpix0
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix0 = 256,          jvpix0 = 256)
       parameter(ihpix2 = ihpix0*ifine, jvpix2 = jvpix0*ifine)
       character*1 rgb0(3,ihpix0, jvpix0) ! RGB data for output
       character*1 rgb2(3,ihpix2, jvpix2) ! RGB data
* local variables
       integer i0, j0, k0
       integer i2, j2, k2
       integer idummy(3)
       integer rdummy

       do 100 j0 = 1, jvpix0
       do 100 i0 = 1, ihpix0
         do k0 = 1, 3
           idummy(k0) = 0
         enddo
         do 110 j2 = (j0-1)*ifine+1, j0*ifine
         do 110 i2 = (i0-1)*ifine+1, i0*ifine
         do 110 k2 = 1, 3
           idummy(k2) = idummy(k2) + ichar(rgb2(k2,i2,j2))
 110     continue
         do k0 = 1, 3
           rdummy = dfloat(idummy(k0))
           rdummy = rdummy / dfloat(ifine**2)
           idummy(k0) = int(rdummy)
           if (idummy(k0) .GT. 255) idummy(k0) = 255
           if (idummy(k0) .LT.   0) idummy(k0) =   0
           rgb0(k0,i0,j0) = char(idummy(k0))
         enddo
 100   continue

       return
       end

*
*** -----------------------------------------------------------------
*
       logical function findfile(rt,nnn)
       implicit none
* interface
       integer rt, nnn
* local
       character*48 flname
       logical      ldummy1, ldummy2
       integer      rt2

       if (rt .LT. 1000)
     &   write(flname,'(''d'',i6.6,''.dat'')') 300000 + nnn
       if ((rt .GE. 1000) .AND. (rt .LT. 10000))
     &   write(flname,'(''d'',i6,''.'',i4)') 300000+nnn,rt
       if (rt .GE. 10000) then
         rt2 = rt - 10000
         write(flname,'(''e'',i7,''.'',i4)') 3000000 + nnn, rt2
       endif
       inquire(file = flname, exist = ldummy1)

       write(flname,'(''u'',i6,''.'',i4)') 300000+nnn,rt
       inquire(file = flname, exist = ldummy2)

       findfile = (ldummy1 .OR. ldummy2)

       return
       end



***=======================================================================
*
       subroutine psout(iout,rgb0,rot2,
     &                  psstr,psix,psjy,numpstr)
       implicit none
* interface
       integer iout
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpixf, jvpixf, iframe, jframe, maxfatp
       parameter(iframe = 1, jframe = 1)
       parameter(maxfatp = iframe * jframe)
       parameter(ihpixf = ihpix0 * iframe, jvpixf = jvpix0 * jframe)
       character*1 rgb0(3,ihpixf,jvpixf) ! RGB data for output (2)
       integer rot2
       character*33 psstr(maxfatp)
       integer psix(maxfatp), psjy(maxfatp), numpstr
* local
       character*1 cmyk(4,ihpixf,jvpixf) ! CYMK data for output
       character*12 flnmps
       integer idummy1, idummy2
       real*8  rdummy1, rdummy2, rdummy3
       real*8  xps0, yps0, fxps, fyps
       integer ix1, ix2, jy1, jy2
       character*2 hx
       integer hexcol
       parameter(hexcol = 30)
       logical filexist
       integer ncol, i0, j0, k0, ic, iy, im, ik
       real*8  aa
* function
       character*2 int2hex2

* determine file name
       idummy2 = 0
       filexist = .true.
       do while (filexist)
         idummy2 = idummy2 + 1
         if (idummy2 .GE. 1000) then
           write(*,*) 'Not available number at 999'
           stop
         endif
         if (rot2 .LT. 1000) then
           idummy1 = 0
         else
           idummy1 = mod(rot2,10000)
         endif
         idummy1 = idummy2 * 10000 + idummy1
         write(flnmps,'(''w'',i7.7,''.eps'')') idummy1
         inquire(file = flnmps, exist = filexist)
       enddo
       write(*,*) 'Now open PS : ', flnmps

* define location of image on paper
       rdummy1 = 740.0 / 540.0 ! 600 X 840 = > 540 X 740
       idummy1 = ihpixf
       idummy2 = jvpixf
       rdummy2 = dfloat(idummy2) / dfloat(idummy1)
       rdummy3 = rdummy2 / rdummy1
       if (rdummy3 .GT. 1.0) then ! Tate-naga
         xps0 =  30.0 + 540.0 / 2.0 * (1.0 - 1.0 / rdummy3)
         fxps = 740.0 / rdummy2
         yps0 =  50.0
         fyps = 740.0
       else
         xps0 =  30.0
         fxps = 540.0
         yps0 =  50.0 + 740.0 / 2.0 * (1.0 - rdummy3)
         fyps = 540.0 * rdummy2
       endif

       open(unit=2,file=flnmps,status='unknown')
* header
       write(2,'(A)') '%!PS-Adobe-3.0 EPSF-3.0'
       ix1 = int(xps0)
       jy1 = int(yps0)
       ix2 = int(xps0 + fxps - 0.001) + 1
       jy2 = int(yps0 + fyps - 0.001) + 1
       write(2,'(''%%BoundingBox:'',4i4)') ix1, jy1, ix2, jy2
       write(2,'(A)') '%%Creator: Keiji Hayashi'
       write(2,'(A)') '%%Title: '// flnmps
       write(2,'(A)') '%%EndComments'
       if (iout .EQ. 6) then
         write(2,'(A)') '/readstring {'
         write(2,'(A)') '  currentfile exch readhexstring pop'
         write(2,'(A)') '} bind def'
         idummy1 = ihpixf
         write(2,'(''/rpicstr'',i5,'' string def'')') idummy1
         write(2,'(''/gpicstr'',i5,'' string def'')') idummy1
         write(2,'(''/bpicstr'',i5,'' string def'')') idummy1
         write(2,'(A)') '%'
         write(2,'(''gsave '',f6.1,1x,f6.1,'' translate'')') xps0,yps0
         write(2,'(2f8.1,'' scale'')') fxps, fyps
         idummy1 = ihpixf
         idummy2 = jvpixf
         write(2,'(2i6,'' 8 '')') idummy1, idummy2
         write(2,'(''[ '',6i6,'' ]'')')
     &     idummy1, 0, 0, -idummy2, 0, idummy2
         write(2,'(A)') '{ rpicstr readstring }'
         write(2,'(A)') '{ gpicstr readstring }'
         write(2,'(A)') '{ bpicstr readstring }'
         write(2,'(A)') 'true 3 colorimage'

         ncol = 0
         do 410 j0 = jvpixf, 1, -1
         do 410 k0 = 1, 3
         do 410 i0 = 1, ihpixf
           idummy1 = ichar(rgb0(k0,i0,j0))
           hx = int2hex2(idummy1)
           ncol = ncol + 1
           if (ncol .EQ. hexcol) then
             write(2,'(A2)') hx
             ncol = 0
           else
             write(2,'(A2,$)') hx
           endif
 410     continue

       else if (iout .EQ. 7) then
         write(2,'(A)') '/readstring {'
         write(2,'(A)') 'currentfile exch readhexstring pop'
         write(2,'(A)') '} bind def'
         idummy1 = ihpixf
          write(2,'(''/picstr'',i6,'' string def'')') idummy1
         write(2,'(A)') '%'
         write(2,'(''gsave '',f6.1,1x,f6.1,'' translate'')') xps0,yps0
         write(2,'(2f8.3,'' scale'')') fxps, fyps
         idummy1 = ihpixf
         idummy2 = jvpixf
         write(2,'(2i6,'' 8 '')') idummy1, idummy2
         write(2,'(''[ '',6i6,'' ]'')')
     &     idummy1, 0, 0, -idummy2, 0, idummy2
         write(2,'(A)') '{ picstr readstring } image'

         ncol = 0
         do 400 j0 = jvpixf, 1, -1
         do 400 i0 = 1, ihpixf
           idummy1 = ichar(rgb0(1,i0,j0))
           rdummy1 = dfloat(idummy1) * 0.299
           idummy1 = ichar(rgb0(2,i0,j0))
           rdummy1 = dfloat(idummy1) * 0.587 + rdummy1
           idummy1 = ichar(rgb0(3,i0,j0))
           rdummy1 = dfloat(idummy1) * 0.114 + rdummy1
           idummy1 = int(rdummy1)
           idummy1 = 255 - idummy1 ! negative

           hx = int2hex2(idummy1)
           ncol = ncol + 1
           if (ncol .EQ. hexcol) then
             write(2,'(A2)') hx
             ncol = 0
           else
             write(2,'(A2,$)') hx
           endif
 400     continue

       else if (iout .EQ. 8) then
         write(2,'(A)') '%%Extensions: CMYK'
         do j0 = 1, jvpixf
         do i0 = 1, ihpixf
           idummy1 = ichar(rgb0(1,i0,j0))
           aa = 255.0D+00 - dfloat(idummy1)
           ic = int(aa * 0.99999D+00)
           idummy1 = ichar(rgb0(2,i0,j0))
           aa = 255.0D+00 - dfloat(idummy1)
           im = int(aa * 0.99999D+00)
           idummy1 = ichar(rgb0(3,i0,j0))
           aa = 255.0D+00 - dfloat(idummy1)
           iy = int(aa * 0.99999D+00)

*           ik = min(255,ic,im,iy)
*           if (ik .EQ. 255) then
*             ic = 255
*             im = 255
*             ik = 255
*           else
*             aa = dfloat(ic - ik) / dfloat(255 - ik)
*             aa = aa * 255.0D+00
*             ic = int(aa)
*             aa = dfloat(im - ik) / dfloat(255 - ik)
*             aa = aa * 255.0D+00
*             im = int(aa)
*             aa = dfloat(iy - ik) / dfloat(255 - ik)
*             aa = aa * 255.0D+00
*             iy = int(aa)
*           endif

           idummy1 = ichar(rgb0(1,i0,j0))
           rdummy1 = dfloat(idummy1) * 0.299
           idummy1 = ichar(rgb0(2,i0,j0))
           rdummy1 = dfloat(idummy1) * 0.587 + rdummy1
           idummy1 = ichar(rgb0(3,i0,j0))
           rdummy1 = dfloat(idummy1) * 0.114 + rdummy1 ! brightness
*           ik = 255 - int(rdummy1)
           ik = int((255.0D+00 - rdummy1) * 0.5D+00)

           cmyk(1,i0,j0) = char(ic)
           cmyk(2,i0,j0) = char(im)
           cmyk(3,i0,j0) = char(iy)
           cmyk(4,i0,j0) = char(ik)
         enddo
         enddo

         write(2,'(A)') '/readstring {'
         write(2,'(A)') ' currentfile exch readhexstring pop'
         write(2,'(A)') '} bind def'
         idummy1 = ihpixf
         write(2,'(''/cpicstr'',i5,'' string def'')') idummy1
         write(2,'(''/mpicstr'',i5,'' string def'')') idummy1
         write(2,'(''/ypicstr'',i5,'' string def'')') idummy1
         write(2,'(''/kpicstr'',i5,'' string def'')') idummy1
         write(2,'(A)') '%'
         write(2,'(''gsave '',f6.1,1x,f6.1,'' translate'')') xps0,yps0
         write(2,'(2f8.1,'' scale'')') fxps, fyps
         idummy1 = ihpixf
         idummy2 = jvpixf
         write(2,'(2i6,'' 8 '')') idummy1, idummy2
         write(2,'(''[ '',6i6,'' ]'')')
     &     idummy1, 0, 0, -idummy2, 0, idummy2
         write(2,'(A)') '{ cpicstr readstring }'
         write(2,'(A)') '{ mpicstr readstring }'
         write(2,'(A)') '{ ypicstr readstring }'
         write(2,'(A)') '{ kpicstr readstring }'
         write(2,'(A)') 'true 4 colorimage'

         ncol = 0
         do 480 j0 = jvpixf, 1, -1
         do 480 k0 = 1, 4
         do 480 i0 = 1, ihpixf
           idummy1 = ichar(cmyk(k0,i0,j0))
           hx = int2hex2(idummy1)
           ncol = ncol + 1
           if (ncol .EQ. hexcol) then
             write(2,'(A2)') hx
             ncol = 0
           else
             write(2,'(A2,$)') hx
           endif
 480     continue

       else
         write(*,*) ' Iout is wrong at jfjf'
         stop
       endif

       write(2,'(A)') '%'
       write(2,'(A)') 'grestore'

       write(2,'(A)') '/Times-Roman findfont 14 scalefont setfont'
       if (iout .EQ. 6) then
         write(2,'(A)') '1 setgray'
       else
         write(2,'(A)') '0 setgray'
       endif

       do i0 = 1, numpstr
         idummy1 = ihpixf
         idummy2 = jvpixf
         if ((psix(i0) .NE. -7777) .AND. (psjy(i0) .NE. -7777))
     &     write(2,'(2f8.2,'' moveto ('',A33,'') show'')')
     &       xps0 + psix(i0) / dfloat(idummy1) * fxps,
     &       yps0 + psjy(i0) / dfloat(idummy2) * fyps,
     &       psstr(i0)
       enddo

       write(2,'(A)') 'showpage'
       write(2,'(A)') '%%Trailer'
       close(2)

       return
       end

***=======================================================================
*
       subroutine bmpout(iout,rgb0,rot2)
       implicit none
* interface
       integer iout
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpixf, jvpixf, iframe, jframe
       parameter(iframe = 1, jframe = 1)
       parameter(ihpixf = ihpix0 * iframe, jvpixf = jvpix0 * jframe)
       character*1 rgb0(3,ihpixf,jvpixf) ! RGB data for output (2)
       integer rot2
* local
       logical filexist
       character*12 fnameout
       integer i0, j0, k0
       integer idummy, idummy1, idummy2
       character*14 frmtstr
       character*54 headmsw
       character*4  byt4
       character*2  byt2

* BMP (24bit depth)... this part works only when width is multiple of 4.
       idummy = mod(ihpixf, 4)
       if (idummy .NE. 0) then
         write(*,*) 'width must be multiple of 4'
         stop
       endif

* determine file name
       if (iout .GE. 1) then
         idummy2 = 0
         filexist = .true.
         do while (filexist)
           idummy2 = idummy2 + 1
           if (idummy2 .GE. 1000) then
             write(*,*) 'Not available number at 999'
             stop
           endif
           if (rot2 .LT. 1000) then
             idummy1 = 0
           else
             idummy1 = mod(rot2,10000)
           endif
           idummy1 = idummy2 ! * 10000 + idummy1
           write(fnameout,'(''sw'',i4.4,''.bmp'')') idummy1
           inquire(file = fnameout, exist = filexist)
         enddo
       endif

       open(unit=2,file=fnameout,status='unknown')
       write(*,*) 'Now generating BMP(24bit) file : ', fnameout
* header (file header ; 1--14 byte)
       headmsw( 1: 2) = 'BM' !             formula
       idummy = 54 + ihpixf * jvpixf * 3 ! file size : header and data
       call num2bit4(idummy,byt4)
       headmsw( 3: 6) = byt4(1:4)
       idummy = 0 !                        must be 0, reserved for future use
       call num2bit2(idummy,byt2)
       headmsw( 7: 8) = byt2(1:2)
       idummy = 0 !                        must be 0, reserved for future use
       call num2bit2(idummy,byt2)
       headmsw( 9:10) = byt2(1:2)
       idummy = 54 !                       must be 54 : total length of header
       call num2bit4(idummy,byt4)
       headmsw(11:14) = byt4(1:4)
* header (bit-map header ; 13--54 byte)
       idummy = 40  !                      must be 40 : length of bit-map header
       call num2bit4(idummy,byt4)
       headmsw(15:18) = byt4(1:4)
       idummy = ihpixf !                   pixel width
       call num2bit4(idummy,byt4)
       headmsw(19:22) = byt4(1:4)
       idummy = jvpixf !                   pixel height
       call num2bit4(idummy,byt4)
       headmsw(23:26) = byt4(1:4)
       idummy = 1 !                        must be 1 : pixel plane ?
       call num2bit2(idummy,byt2)
       headmsw(27:28) = byt2(1:2)
       idummy = 24 !                       must be 24 here : color depth in bit.
       call num2bit2(idummy,byt2)
       headmsw(29:30) = byt2(1:2)
       idummy = 0 !                        must be 0 here : compression method index
       call num2bit4(idummy,byt4)
       headmsw(31:34) = byt4(1:4)
       idummy = 0 !                        must be 0 here : file size if compressed
       call num2bit4(idummy,byt4)
       headmsw(35:38) = byt4(1:4)
       idummy = 0 !                        arbit. : pixel per meter, horizontal
       call num2bit4(idummy,byt4)
       headmsw(39:42) = byt4(1:4)
       idummy = 0 !                        arbit. : pixel per meter, vertical
       call num2bit4(idummy,byt4)
       headmsw(43:46) = byt4(1:4)
       idummy = 0 !                        may be zero here : num. of color used
       call num2bit4(idummy,byt4)
       headmsw(47:50) = byt4(1:4)
       idummy = 0 !                        arbit. : num. of important color
       call num2bit4(idummy,byt4)
       headmsw(51:54) = byt4(1:4)

* writting header part
       write(2,'(a54,$)') headmsw(1:54)
* image data
       idummy = ihpixf * jvpixf * 3
       write(frmtstr,'(''('',i8.8,''A,$)'')') idummy
       write(2,fmt=frmtstr)
     &     (((rgb0(k0,i0,j0),k0=3,1,-1),i0=1,ihpixf),j0=1,jvpixf)
       close(2)

       return
       end

* --------------------------------------
* convert number to 8-bit characters
* --------------------------------------

       subroutine num2bit4(inum,byt4) ! MIND lower bit comes first, Intel x86 ......
       implicit none
       integer inum
       character*4 byt4
       integer idummy1, idummy2
       idummy1 = inum
       idummy2 = idummy1 / 256**3
       byt4(4:4) = char(idummy2)
       idummy1 =-idummy2 * 256**3 +idummy1
       idummy2 = idummy1 / 256**2
       byt4(3:3) = char(idummy2)
       idummy1 =-idummy2 * 256**2 +idummy1
       idummy2 = idummy1 / 256
       byt4(2:2) = char(idummy2)
       idummy1 =-idummy2 * 256    +idummy1
       byt4(1:1) = char(idummy1)
       return
       end subroutine

* ------

       subroutine num2bit2(inum,byt2) ! MIND lower bit comes first, Intel x86 ......
       implicit none
       integer inum
       character*2 byt2
       integer idummy1, idummy2
       idummy1 = inum
       idummy2 = idummy1 / 256
       byt2(2:2) = char(idummy2)
       idummy1 =-idummy2 * 256 + idummy1
       byt2(1:1) = char(idummy1)
       return
       end subroutine


***=======================================================================
*
       subroutine rgbout(iout,rgb0,rot2)
       implicit none
* interface
       integer iout
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpixf, jvpixf, iframe, jframe
       parameter(iframe = 1, jframe = 1)
       parameter(ihpixf = ihpix0 * iframe, jvpixf = jvpix0 * jframe)
       character*1 rgb0(3,ihpixf,jvpixf) ! RGB data for output (2)
       integer rot2
* local
       character*48 fnameout, fnamefal
       logical filexist
       integer i0, j0, k0
       character*2 hx
       integer idummy1, idummy2, ncol
       integer ip
       integer hexcol, ppmcol
       parameter(hexcol = 50, ppmcol = 60)
       character*14 frmtstr
* dummy
*       real*8   rdummy
* function
       character*2 int2hex2

* determine file name
       if (iout .GE. 1) then
         idummy2 = 0
         filexist = .true.
         do while (filexist)
           idummy2 = idummy2 + 1
           if (idummy2 .GE. 1000) then
             write(*,*) 'Not available number at 999'
             stop
           endif
           if (rot2 .LT. 1000) then
             idummy1 = 0
           else
             idummy1 = mod(rot2,10000)
           endif
           idummy1 = idummy2 ! * 10000 + idummy1
           if (iout .EQ. 1) then
             write(fnameout,'(''sw'',i4.4,''.rgb'')') idummy1
             write(fnamefal,'(''sw'',i4.4,''.fal'')') idummy1
           endif
           if (iout .EQ. 2)
     &       write(fnameout,'(''sw'',i4.4,''.ppm'')') idummy1
           if (iout .EQ. 3)
     &       write(fnameout,'(''sw'',i4.4,''.ppm'')') idummy1
           if (iout .EQ. 4)
     &       write(fnameout,'(''sw'',i4.4,''.hex'')') idummy1
           if (iout .EQ. 5)
     &       write(fnameout,'(''sw'',i4.4,''.pgm'')') idummy1
           inquire(file = fnameout, exist = filexist)
         enddo
       endif

* RGB
       if (iout .EQ. 1) then
         open(unit=2,file=fnamefal,status='unknown')
         write(2,'(2i5)') ihpixf, jvpixf
         close(2)
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now open RGB file : ', fnameout
         idummy1 = ihpixf * jvpixf * 3
         write(frmtstr,'(''($'',i8.8,''A,$)'')') idummy1
         write(2,fmt=frmtstr)
     &     (((rgb0(k0,i0,j0),k0=1,3),i0=1,ihpixf),j0=1,jvpixf)
         close(2)
       endif
       if (iout .EQ. 0) then
         idummy1 = ihpixf * jvpixf * 3
         write(frmtstr,'(''($'',i8.8,''A,$)'')') idummy1
         write(*,fmt=frmtstr)
     &     (((rgb0(k0,i0,j0),k0=1,3),i0=1,ihpixf),j0=1,jvpixf)
       endif

** PPM P3
*       if (iout .EQ. 2) then
*         open(unit=2,file=fnameout,status='unknown')
*         write(*,*) 'Now open PPM (P3) file : ', fnameout
*         write(2,'(''P3'', 2(1x,i4),'' 255'')') ihpixf, jvpixf
*         ncol = 0
*         do 10 j0 = jvpixf, 1, -1
*         do 10 i0 = 1, ihpixf, 1
*         do 10 k0 = 1, 3
*           idummy1 = ichar(rgb0(k0,i0,j0))
*           if (idummy1 .EQ. 0) then
*             idummy2 = 1
*           else
*             rdummy  = dfloat(idummy1) + 0.00001
*             rdummy  = dlog10(rdummy) + 1
*             idummy2 = int(rdummy)
*           endif
*           ncol = ncol + idummy2 + 1
*           if (ncol .LE. ppmcol) then
*             if (idummy2 .EQ. 1) frmtstr = '(1x,i1,$)'
*             if (idummy2 .EQ. 2) frmtstr = '(1x,i2,$)'
*             if (idummy2 .EQ. 3) frmtstr = '(1x,i3,$)'
*           else
*             if (idummy2 .EQ. 1) frmtstr = '(1x,i1)'
*             if (idummy2 .EQ. 2) frmtstr = '(1x,i2)'
*             if (idummy2 .EQ. 3) frmtstr = '(1x,i3)'
*             ncol = 0
*           endif
*           write(2,fmt=frmtstr) idummy1
* 10      continue
*         close(2)
*       endif

* PPM P3
       if (iout .EQ. 2) then
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now open PPM (P3) file : ', fnameout
         write(2,'(''P3'', 2(1x,i4),'' 255'')') ihpixf, jvpixf
         ncol = 0
         do 10 j0 = jvpixf, 1, -1
         do 10 i0 = 1, ihpixf, 1
         do 10 k0 = 1, 3
           idummy1 = ichar(rgb0(k0,i0,j0))
           ncol = ncol + 4
           if (ncol .LT. ppmcol) then
             write(2,fmt='(1x,i3,$)') idummy1
           else
             write(2,fmt='(1x,i3)') idummy1
             ncol = 0
           endif
 10      continue
         write(2,'(A)') ' '
         close(2)
       endif


* PPM P6 of "cinema"-mode..
*       if (iout .EQ. 3) then
*         open(unit=2,file=fnameout,status='unknown')
*         write(*,*) 'Now open PPM (P6) file : ', fnameout
*** header part for running on MS-DOS
**         write(2,'(''P6'', 2(1x,i4),'' 255 '',$)') ihpixf, jvpixf
*** header part for IDL readable version
*         write(2,'(A)') 'P6'
*         write(2,'(A)') '# Creator : Fortran code by Keiji Hayashi'
*         write(2,'(2(1x,i4),'' 255 '',$)') ihpixf, jvpixf * 3 / 4
*         idummy1 = ihpixf * jvpixf * 3 * 3 / 4
*         write(frmtstr,'(''('',i8.8,''A,$)'')') idummy1
*         write(2,fmt=frmtstr)
*     +     (((rgb0(k0,i0,j0),k0=1,3),i0=1,ihpixf),
*     +                               j0=jvpixf*7/8,jvpixf/8+1,-1)
*         close(2)
*       endif
* PPM P6
       if (iout .EQ. 3) then
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now open PPM (P6) file : ', fnameout
** header part for running on MS-DOS
*         write(2,'(''P6'', 2(1x,i4),'' 255 '',$)') ihpixf, jvpixf
** header part for IDL readable version
         write(2,'(A)') 'P6'
*         write(2,'(A)') '# Creator : Fortran code by Keiji Hayashi'
*         write(2,'(2(1x,i4),'' 255 '',$)') ihpixf, jvpixf
         write(2,'(2(1x,i4))') ihpixf, jvpixf
         write(2,'(i4)') 255
         idummy1 = ihpixf * jvpixf * 3
         write(frmtstr,'(''('',i8.8,''A,$)'')') idummy1
         write(2,fmt=frmtstr) ! write(2,fmt='(786432A)')
     &     (((rgb0(k0,i0,j0),k0=1,3),i0=1,ihpixf),j0=jvpixf,1,-1)
         close(2)
       endif

* HEX
       if (iout .EQ. 4) then
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now open HEX file : ', fnameout
         write(2,'(3i8)') ihpixf, jvpixf, hexcol       ! header
         ncol = 0
         do 400 j0 = 1, jvpixf
         do 400 i0 = 1, ihpixf
         do 400 k0 = 1, 3
           idummy1 = ichar(rgb0(k0,i0,j0))
           hx = int2hex2(idummy1)
           ncol = ncol + 1
           if (ncol .EQ. hexcol) then
             write(2,'(A2)') hx
             ncol = 0
           else
             write(2,'(A2,$)') hx
           endif
 400     continue
         do ip = ncol + 1, hexcol - 1  ! Padding
           write(2,'(A2,$)') '--'
         enddo
         write(2,'(A2)') '--'
         close(2)
       endif

* PGM P5 : gray and negative
       if (iout .EQ. 5) then
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now open PGM (P5) file : ', fnameout
         write(2,'(''P5'', 2(1x,i4),'' 255 '',$)') ihpixf, jvpixf
         idummy1 = ihpixf * jvpixf
         write(frmtstr,'(''('',i8.8,''A,$)'')') idummy1
         do 500 j0 = jvpixf, 1, -1
         do 500 i0 = 1, ihpixf
           idummy1 = 0
           do k0 = 1, 3
             idummy1 = idummy1 + ichar(rgb0(k0,i0,j0))
           enddo
           idummy1 = 255 - idummy1 / 3
           rgb0(1,i0,j0) = char(idummy1)
 500     continue
         write(2,fmt=frmtstr)
     &     ((rgb0(1,i0,j0),i0=1,ihpixf),j0=jvpixf,1,-1)
         close(2)
       endif

* initialize and clean up
       do 20 j0 = 1, jvpixf
       do 20 i0 = 1, ihpixf
       do 20 k0 = 1, 3
         rgb0(k0,i0,j0) = char(0)
 20    continue

       return
       end



*
*** -------------
*
       character*2 function int2hex2(i)
       implicit none
* interface
       integer i
* local
       character*2 sdummy
       integer  idummy(2)
       integer  m

       idummy(1) =     i/16
       idummy(2) = mod(i,16)

       do m = 1, 2
         if (idummy(m) .LT. 10) then
           idummy(m) = idummy(m) + 48
         else
           idummy(m) = idummy(m) + 55
         endif
         sdummy(m:m) = char(idummy(m))
       enddo

       int2hex2 = sdummy

       return
       end



*
***=======================================================================
*
       subroutine surf2rgb(isunsurf,v0)
       implicit  none
* interface
       integer isunsurf
       real*8  v0
* local
       integer j, k, iii, m
       real*8  aa
       integer ic, jc          ! pixel potition of center of project plane
       real*8  bbb, ccc
       real*8  xp1, yp1
       integer ix1, jy1
       integer red, gre, blu
       integer ijk(3)
       integer j1, k1, j2, k2
       real*8  rtp(3), xyz(3), obj2lght(3)
       real*8  xyzeye(3), nprj(3)
       real*8  surfbr, dstnt2
       real*8  dth1, dth2, dph1, dph2
       real*8  maxmag, minmag ! actually. not only mag....
* dummy
       integer idummy1, idummy2, idummy3
       real*8  rdummy1
       logical ldummy
* constant
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* functions
       real*8  naiseki3
       logical sameside
* common
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       common  /images/ dstnt, rgb2
       common  /coord/  rr, theta, phi
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

* pixel position on the projection plane
       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,aa)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo

       maxmag = -1.0D+10
       minmag =  1.0D+10
       if (isunsurf .EQ. 1) then ! Br
         do 112 k1 = 1, kk
         do 112 j1 = 0, jj
           surfbr = br1(0,j1,k1)
           surfbr = dabs(surfbr)
           if (surfbr .GT. maxmag) maxmag = surfbr
 112     continue
       else if (isunsurf .EQ. 6) then ! P
         do 116 k1 = 1, kk
         do 116 j1 = 0, jj
           surfbr = pg1(0,j1,k1)
           surfbr = dabs(surfbr)
           if (surfbr .GT. maxmag) maxmag = surfbr
           if (surfbr .LT. minmag) minmag = surfbr
 116     continue
       else if (isunsurf .EQ. 5) then ! N
         do 115 k1 = 1, kk
         do 115 j1 = 0, jj
           surfbr = ro1(0,j1,k1)
           surfbr = dabs(surfbr)
           if (surfbr .GT. maxmag) maxmag = surfbr
           if (surfbr .LT. minmag) minmag = surfbr
 115     continue
       else if (isunsurf .EQ. 4) then ! T
         do 114 k1 = 1, kk
         do 114 j1 = 0, jj
           surfbr = pg1(0,j1,k1) / ro1(0,j1,k1)
           surfbr = dabs(surfbr)
           if (surfbr .GT. maxmag) maxmag = surfbr
           if (surfbr .LT. minmag) minmag = surfbr
 114     continue
       else
         do 113 k1 = 1, kk
         do 113 j1 = 0, jj
           surfbr = ur1(0,j1,k1) ! Vr
           surfbr = dabs(surfbr)
           if (surfbr .GT. maxmag) maxmag = surfbr
 113     continue
       endif

* drawing the solar surface
       aa = rbottom / ul
       iii = int(aa * 10.0)
       if (iii .LE. 0) iii = 1
       do 111 k = 0, (iii * 2)
       do 111 j = 0, iii
         rtp(1) = rr(0) * 1.002 ! rbottom
         rtp(2) = pi * dfloat(j) / dfloat(iii)
         rtp(3) = pi * dfloat(k) / dfloat(iii)
         call rtp2xyz(rtp,xyz)

         call prjpnt1(xyz,xp1,yp1,dstnt2)
         call pxy2ij(xp1,yp1,ul,ix1,jy1,ic,jc)

         ldummy = sameside(xyzcp,nprj,xyz,xyzeye)

         if ((ix1 .GT. 0) .AND. (ix1 .LE. ihpix2) .AND.
     &       (jy1 .GT. 0) .AND. (jy1 .LE. jvpix2) .AND.
     &       (dstnt2 .LT. dstnt(ix1,jy1)) .AND. (.NOT. ldummy)) then

           dstnt(ix1,jy1) = dstnt2

* estimate the brightness of the solar surface : colar table
           if (isunsurf .EQ. 3) then ! Vr
             call rtp2ijk(rtp,ijk)
             j1 = ijk(2)
             j2 = ijk(2) + 1
             k1 = ijk(3)
             k2 = ijk(3) + 1
             dth1 =(rtp(2) - theta(j1))/(theta(j2)-theta(j1))
             dth2 =(theta(j2) - rtp(2))/(theta(j2)-theta(j1))
             dph1 =(rtp(3) - phi(k1))/(phi(k2)-phi(k1))
             dph2 =(phi(k2) - rtp(3))/(phi(k2)-phi(k1))
             surfbr = (ur1(0,j1,k1)*dth2+ur1(0,j2,k1)*dth1)*dph2
     &               +(ur1(0,j1,k2)*dth2+ur1(0,j2,k2)*dth1)*dph1
             surfbr = surfbr * v0 / 1.0D+05
             rdummy1 = 1.0
             call ipscolor(red,gre,blu,surfbr,rdummy1)
           else if (isunsurf .EQ. 6) then ! P
             call rtp2ijk(rtp,ijk)
             j1 = ijk(2)
             j2 = ijk(2) + 1
             k1 = ijk(3)
             k2 = ijk(3) + 1
             dth1 =(rtp(2) - theta(j1))/(theta(j2)-theta(j1))
             dth2 =(theta(j2) - rtp(2))/(theta(j2)-theta(j1))
             dph1 =(rtp(3) - phi(k1))/(phi(k2)-phi(k1))
             dph2 =(phi(k2) - rtp(3))/(phi(k2)-phi(k1))
             surfbr = (pg1(0,j1,k1)*dth2
     &                +pg1(0,j2,k1)*dth1)*dph2
     &               +(pg1(0,j1,k2)*dth2
     &                +pg1(0,j2,k2)*dth1)*dph1
             surfbr = (surfbr - minmag) / (maxmag - minmag + 0.000001)
             surfbr = surfbr**0.5D+00 * 255.0D+00
             if (surfbr .LT. 191.0D+00) then
               surfbr = surfbr * 1.33
               red = int(surfbr) ! purple
               gre = 0
               blu = int(surfbr)
             else
               surfbr = (surfbr - 191.0D+00) * 4.0
               red = 255
               gre = int(surfbr)
               blu = 255
             endif
           else if (isunsurf .EQ. 5) then ! N
             call rtp2ijk(rtp,ijk)
             j1 = ijk(2)
             j2 = ijk(2) + 1
             k1 = ijk(3)
             k2 = ijk(3) + 1
             dth1 =(rtp(2) - theta(j1))/(theta(j2)-theta(j1))
             dth2 =(theta(j2) - rtp(2))/(theta(j2)-theta(j1))
             dph1 =(rtp(3) - phi(k1))/(phi(k2)-phi(k1))
             dph2 =(phi(k2) - rtp(3))/(phi(k2)-phi(k1))
             surfbr = (ro1(0,j1,k1)*dth2
     &                +ro1(0,j2,k1)*dth1)*dph2
     &               +(ro1(0,j1,k2)*dth2
     &                +ro1(0,j2,k2)*dth1)*dph1
             surfbr = (surfbr - minmag) / (maxmag - minmag + 0.000001)
             surfbr = surfbr**0.5D+00 * 255.0D+00
             if (surfbr .LT. 127.0D+00) then
               surfbr = surfbr * 2.0
               red = 0
               gre = 0
               blu = int(surfbr) ! blue
             else
               surfbr = (surfbr - 127.0) * 2.0
               red = int(surfbr)
               gre = int(surfbr)
               blu = 255
             endif
           else if (isunsurf .EQ. 4) then ! T
             call rtp2ijk(rtp,ijk)
             j1 = ijk(2)
             j2 = ijk(2) + 1
             k1 = ijk(3)
             k2 = ijk(3) + 1
             dth1 =(rtp(2) - theta(j1))/(theta(j2)-theta(j1))
             dth2 =(theta(j2) - rtp(2))/(theta(j2)-theta(j1))
             dph1 =(rtp(3) - phi(k1))/(phi(k2)-phi(k1))
             dph2 =(phi(k2) - rtp(3))/(phi(k2)-phi(k1))
             surfbr = (pg1(0,j1,k1)/ro1(0,j1,k1)*dth2
     &                +pg1(0,j2,k1)/ro1(0,j2,k1)*dth1)*dph2
     &               +(pg1(0,j1,k2)/ro1(0,j1,k2)*dth2
     &                +pg1(0,j2,k2)/ro1(0,j2,k2)*dth1)*dph1
             surfbr = (surfbr - minmag) / (maxmag - minmag + 0.000001)
             surfbr = surfbr**0.5D+00 * 255.0D+00
             if (surfbr .LT. 127.0D+00) then
               surfbr = surfbr * 2.0
               red = int(surfbr) ! red
               gre = 0
               blu = 0
             else
               surfbr = (surfbr - 127.0) * 2.0
               red = 255 !  red !!
               gre = int(surfbr)
               blu = int(surfbr)
             endif
           else if (isunsurf .EQ. 1) then ! Br
             call rtp2ijk(rtp,ijk)
             j1 = ijk(2)
             j2 = ijk(2) + 1
             k1 = ijk(3)
             k2 = ijk(3) + 1
             dth1 =(rtp(2) - theta(j1))/(theta(j2)-theta(j1))
             dth2 =(theta(j2) - rtp(2))/(theta(j2)-theta(j1))
             dph1 =(rtp(3) - phi(k1))/(phi(k2)-phi(k1))
             dph2 =(phi(k2) - rtp(3))/(phi(k2)-phi(k1))
             surfbr = (br1(0,j1,k1)*dth2+br1(0,j2,k1)*dth1)*dph2
     &               +(br1(0,j1,k2)*dth2+br1(0,j2,k2)*dth1)*dph1
             if (surfbr .GE. 0.0D+00) then
               surfbr = (surfbr /  maxmag)**0.3D+00 * 255.0D+00
               red = 255 - int(surfbr)
               gre = 255 - int(surfbr)
               blu = 255
             else
               surfbr = dabs(surfbr)
               surfbr = (surfbr /  maxmag)**0.3D+00 * 255.0D+00
               red = 255
               gre = 255 - int(surfbr)
               blu = 255 - int(surfbr)
             endif
           else ! yellow
             do m = 1, 3
               obj2lght(m) = xyzlght(m) - xyz(m)
             enddo
             call normvec3(obj2lght,rdummy1)
             call normvec3(xyz,rdummy1)
             bbb = naiseki3(obj2lght,xyz)
             ccc = 0.5D+00 * (bbb + 1.0D+00) * 255.0D+00
             red = 255
             gre = 255
             blu = int(ccc)
* white ball !!
*             red = int(ccc)
*             gre = int(ccc)
*             blu = int(ccc)
           endif
* put color value to rgb-data-array
           idummy1 = red
           idummy2 = gre
           idummy3 = blu
           if (idummy1 .GT. 255) idummy1 = 255
           if (idummy2 .GT. 255) idummy2 = 255
           if (idummy3 .GT. 255) idummy3 = 255
           if (idummy1 .LT.   0) idummy1 =   0
           if (idummy2 .LT.   0) idummy2 =   0
           if (idummy3 .LT.   0) idummy3 =   0
           rgb2(1,ix1,jy1) = char(idummy1)
           rgb2(2,ix1,jy1) = char(idummy2)
           rgb2(3,ix1,jy1) = char(idummy3)
         endif
 111   continue

       return
       end



*
***=======================================================================
*
       subroutine intg2rgb(isunsurf,ivar,rotref,nref,ienhintg,length0)
       implicit  none
* interface
       integer isunsurf,ivar,rotref,nref,ienhintg
       real*8  length0
* local
       integer istart, iend, jstart, jend
       integer ic, jc          ! pixel potition of center of project plane
       real*8  fff
       real*8  aa, bb, cc, tt, tt0, rrnear
       integer clrbtm, clrmdl  ! minimum and miduim color value of current sheet
       real*8  zza, zzb, dss
       real*8  sumsum, dpi, sumdss, rada, radb, fstep
       real*8  dstnta, dstntb
       integer sumtime
       real*8  rbottom2
       integer i, j, m, i0, j0
       integer red, gre, blu
       real*8  xyza(3), xyz0(3), rtpa(3), rtpb(3), xyzb(3)
       real*8  lofsxyz(3), xyzeye(3), nprj(3)
       real*8  xpt, ypt
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       parameter(clrbtm = 100, clrmdl = 175) ! color table control index
       real*8  fnewkirk(0:15)
       real*8  rupper2
* dummy
       integer idummy1, idummy2, idummy3
       real*8  xyzdummy(3)
       real*8  rdummy1, rdummy2, rdummy3
       real*8  rrnear2
* functions
       real*8  naiseki3, absvec3, parkvr
       real*8  sclatrtp
* temporal
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       real*8   intgval(ihpix2, jvpix2)
       real*8   intgval2(ihpix2, jvpix2)
       real*8   maxig, minig, aveig, avei2, devig
       real*8   maxig2, minig2
       integer  numig
       real*8   intghst(ihpix2), avecirc(ihpix2), ave2circ(ihpix2)
       real*8   maxcirc(ihpix2), mincirc(ihpix2)
       real*8   devhst(ihpix2)
       integer  numcirc(ihpix2), numhist(0:255), sumhist(0:255)
       real*8   parkden ! density of Parker, at particular hight.
* common
       integer iout, irlog
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  ro0(0:ii+1,0:jj,0:kk+1)
       real*8  pg0(0:ii+1,0:jj,0:kk+1)
       real*8  ur0(0:ii+1,0:jj,0:kk+1)
       real*8  ut0(0:ii+1,0:jj,0:kk+1)
       real*8  up0(0:ii+1,0:jj,0:kk+1)
       real*8  br0(0:ii+1,0:jj,0:kk+1)
       real*8  bt0(0:ii+1,0:jj,0:kk+1)
       real*8  bp0(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       real*8  ro3(0:ii+1,0:jj,0:kk+1)
       real*8  pg3(0:ii+1,0:jj,0:kk+1)
       real*8  ur3(0:ii+1,0:jj,0:kk+1)
       real*8  ut3(0:ii+1,0:jj,0:kk+1)
       real*8  up3(0:ii+1,0:jj,0:kk+1)
       real*8  br3(0:ii+1,0:jj,0:kk+1)
       real*8  bt3(0:ii+1,0:jj,0:kk+1)
       real*8  bp3(0:ii+1,0:jj,0:kk+1)
*
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom, gamma
       common  /cgcfg/  iout, irlog
       common  /coord/  rr, theta, phi
       common  /var01/  ro0, pg0
       common  /var02/  ur0, ut0, up0
       common  /var04/  br0, bt0, bp0
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /var31/  ro3, pg3
       common  /var32/  ur3, ut3, up3
       common  /var34/  br3, bt3, bp3
       common  /images/ dstnt, rgb2
       common  /sclscl/ sclr
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom
       common  /gasgas/ gamma

* newkirk filter...approx.
       fnewkirk(0)  = 3.10 ! 2.90
       fnewkirk(1)  = 2.90
       fnewkirk(2)  = 1.85
       fnewkirk(3)  = 1.27
       fnewkirk(4)  = 0.95
       fnewkirk(5)  = 0.77
       fnewkirk(6)  = 0.62
       fnewkirk(7)  = 0.51
       fnewkirk(8)  = 0.41
       fnewkirk(9)  = 0.33
       fnewkirk(10) = 0.25
       fnewkirk(11) = 0.18
       fnewkirk(12) = 0.12
       fnewkirk(13) = 0.08
       fnewkirk(14) = 0.02
       fnewkirk(15) = 0.02
* initialize
       do j = 1, jvpix2
         do i = 1, ihpix2
           intgval(i,j)  = 0.0D+00
           intgval2(i,j) = 0.0D+00
         enddo
       enddo

       if ((isunsurf .EQ. 1) .OR.
     &     ((isunsurf .GE. 3) .AND.
     &      (isunsurf .LE. 7))) then
         rbottom2 = rbottom ! ocult
       else
         rbottom2 = 0.0D+00
       endif
*       rbottom2 = 2.0D+00 / length0  ! ocult exactly 2 Rs.... like C2


       aa = length0

* set variable to be drawn
       call setsclr(ivar,rotref,nref)

* pixel position on the projection plane
       istart = 1 ! (ihpix2 - jvpix2) / 2 + 1
       iend   = ihpix2 ! istart + jvpix2
       jstart = 1
       jend   = jvpix2
       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,aa)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo

* the tolerance at the singlar : rotational axis (theta = 0, pi)
       dpi  = pi / dfloat(jj)

       maxig = -1.0D+40
       minig =  1.0D+40

** Get the RGB value of each pixel
       do 777 i = istart, iend
       do 777 j = jstart, jend

         intgval(i,j) = 0.0D+00

* Get the point on the projection plane at which line of sight starts
         aa = ul * dfloat(i-ic)
         bb = ul * dfloat(j-jc)
         do m = 1, 3
           xyz0(m) = xyzcp(m) + aa * xyzxp(m) + bb * xyzyp(m)
         enddo
* unit vector of line of sight
         do m = 1, 3
           lofsxyz(m)= xyz0(m) - xyzeye(m)
         enddo
         call normvec3(lofsxyz, tt0)

* Calculation will NOT be performed when
*   line of sight runs always beyond the sphere of radius of R_upper
         aa = - naiseki3(lofsxyz,xyzeye)
         do m = 1, 3
           xyzdummy(m) = xyzeye(m) + aa * lofsxyz(m)
         enddo
         rrnear = absvec3(xyzdummy)

         if ((rrnear .LT. rupper) .AND. (rrnear .GT. rbottom2)) then

* Rupper for los integration
           if (rrnear .GT. rbottom) then ! corona graph. line
             rupper2 = rupper * 1.5D+00
             if (rupper2 .GT. rr(ii)) rupper2 = rr(ii)
           else !                          EIT like ....
             rupper2 = rbottom * 1.1D+00
             if (rupper2 .GT. rupper) rupper2 = rupper
           endif

           if (ienhintg .EQ. 6) then
             rrnear2 = rrnear / rbottom
             if (rrnear2 .LT. 1.001D+00) rrnear2 = 1.001 * rbottom
             parkden = 1.0D+00 / parkvr(rrnear,gamma) / rrnear**2
           else
             parkden = 1.0D+00
           endif
* Find start point of integration path.
* At this point, line of sight cross the nearest boundray sphere.
           bb = naiseki3(xyzeye,lofsxyz)
           cc = bb**2 - absvec3(xyzeye) * absvec3(xyzeye) + rupper2**2
           tt = - bb - dsqrt(cc)
           if (tt .LT. tt0) tt = tt0
           tt = tt * 1.001
           do m = 1, 3
             xyzb(m) = xyzeye(m) + tt * lofsxyz(m)
           enddo

           call xyz2rtp(xyzb,rtpb)
           call prjpnt1(xyzb,xpt,ypt,dstntb)
           radb = rtpb(1)
           zzb  =  sclatrtp(rtpb)

           sumtime = 0
           sumsum = 0.0D+00
           sumdss = 0.0D+00
           do while ((radb .LT. rupper2) .AND.
     &               (radb .GT. rbottom * 1.0001D+00) .AND.
     &                dstntb .LT. dstnt(i,j))
* sum up
             sumtime = sumtime + 1
             dss  = radb * pi / dfloat(jj) / 1.0D+01
* set advancing step
             do m = 1, 3
               xyza(m) = xyzb(m) + dss * lofsxyz(m)  ! newstep
             enddo
             call xyz2rtp(xyza,rtpa)
             rada = rtpa(1)
             call prjpnt1(xyza,xpt,ypt,dstnta)
             fstep = 1.0D+00
             if (dstnta .LT. dstnt(i,j)) then
               rdummy1 = (dstnt(i,j) - dstntb) / (dstnta - dstntb)
             else
               rdummy1 = 1.0D+00
             endif
             if (rada .GE. rupper2) then
               rdummy2 = (rupper2 - radb) / (rada - radb)
             else
               rdummy2 = 1.0D+00
             endif
             if (rada .LE. rbottom) then
               rdummy3 = (rbottom - radb) / (rada - radb)
             else
               rdummy3 = 1.0D+00
             endif
             if (rdummy1 .LT. fstep) fstep = rdummy1
             if (rdummy2 .LT. fstep) fstep = rdummy2
             if (rdummy3 .LT. fstep) fstep = rdummy3
*
             dss = dss * fstep
             do m = 1, 3
               xyza(m) = xyzb(m) + dss * lofsxyz(m)
             enddo
             call xyz2rtp(xyza,rtpa)
             rada = rtpa(1)
             call prjpnt1(xyza,xpt,ypt,dstnta)
*
             sumdss = sumdss + dss
             zza =  sclatrtp(rtpa)
             sumsum = sumsum + (zza + zzb) * dss / 2.0D+00
* copy
             do m = 1, 3
               xyzb(m) = xyza(m)
             enddo
             radb   = rada
             dstntb = dstnta
             zzb    = zza
           enddo !             End of integ.
*
           if (ivar .EQ. 8) sumsum = sumsum * rrnear**1.5D+00 ! IPS
           if (ienhintg .EQ. 6) sumsum = sumsum / parkden
           if (sumdss .GT. 1.0D-20) then
             intgval(i,j) = sumsum
           else
             intgval(i,j) = 0.0D+00
           endif
         else
           intgval(i,j) = 0.0D+00
         endif

         if (intgval(i,j) .GT. maxig) maxig = intgval(i,j)
         if (intgval(i,j) .LT. minig) minig = intgval(i,j)

 777   continue

       write(11,*) ' Max Min orginal Image = ', maxig, minig

* enhancing
       if (ienhintg .EQ. 1) then ! enhance at each circle bin
         do i = 1, ihpix2
           intghst(i)  = 0.0D+00
           avecirc(i)  = 0.0D+00
           ave2circ(i) = 0.0D+00
           devhst(i)   = 0.0D+00
           maxcirc(i)  =-1.0D+30
           mincirc(i)  = 1.0D+30
           numcirc(i)  = 0
         enddo
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           aa = dfloat(i0-ic)
           bb = dfloat(j0-jc)
           cc = dsqrt(aa**2 + bb**2)
           i = int(cc) + 1
           numcirc(i)  = numcirc(i) + 1
           avecirc(i)  = avecirc(i) + intgval(i0,j0)
           ave2circ(i) = ave2circ(i)+ intgval(i0,j0)**2
           i = int(cc)
           numcirc(i)  = numcirc(i) + 2
           avecirc(i)  = avecirc(i) + intgval(i0,j0)    * 2.0
           ave2circ(i) = ave2circ(i)+ intgval(i0,j0)**2 * 2.0
           i = int(cc) - 1
           numcirc(i)  = numcirc(i) + 1
           avecirc(i)  = avecirc(i) + intgval(i0,j0)
           ave2circ(i) = ave2circ(i)+ intgval(i0,j0)**2
           if (intgval(i0,j0) .GT. maxcirc(i))
     &                             maxcirc(i) = intgval(i0,j0)
           if (intgval(i0,j0) .LT. mincirc(i))
     &                             mincirc(i) = intgval(i0,j0)
         enddo
         enddo
         do i = 1, ihpix2
           if (numcirc(i) .GT. 0) then
             avecirc(i)  = avecirc(i)  / dfloat(numcirc(i))
             ave2circ(i) = ave2circ(i) / dfloat(numcirc(i))
             cc = avecirc(i)**2 - ave2circ(i)
             cc = dabs(cc)
             cc = dsqrt(cc)
             devhst(i) = cc
           else
             avecirc(i)  = 0.0D+00
             ave2circ(i) = 0.0D+00
             devhst(i)   = 1.0D+00
             maxcirc(i)  = 1.0D+00
             mincirc(i)  = 0.0D+00
           endif
         enddo
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           aa = dfloat(i0-ic)
           bb = dfloat(j0-jc)
           cc = dsqrt(aa**2 + bb**2)
           idummy1 = int(cc)
*           cc = 0.5D+00
*     +        +(intgval(i0,j0) - avecirc(idummy1)) ! ave, std
*     +        / devhst(idummy1) / 4.0D+00
           cc = (intgval(i0,j0)    - mincirc(idummy1))
     &        /  (maxcirc(idummy1) - mincirc(idummy1)) ! max, min
           if (cc .GT. 1.0D+00) cc = 1.0D+00
           if (cc .LT. 0.0D+00) cc = 0.0D+00
           intgval2(i0,j0) = cc
         enddo
         enddo
       else if (ienhintg .EQ. 2) then ! Dev. entire
         aveig = 0.0D+00
         avei2 = 0.0D+00
         numig = 0
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           cc = intgval(i0,j0)
           if (dabs(cc) .GT. 1.0D-10) then
             numig = numig + 1
             aveig = aveig + cc
             avei2 = avei2 + cc**2
           endif
         enddo
         enddo
         if (numig .GT. 0) then
           aveig = aveig / dfloat(numig)
           avei2 = avei2 / dfloat(numig)
           devig = dabs(aveig**2 - avei2)
           devig = dsqrt(devig)
         else
           aveig = 1.0D+00
           avei2 = 1.0D+00
           devig = 1.0D+00
         endif
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           cc = 0.5D+00
     &        + (intgval(i0,j0) - aveig) / devig / 4.0D+00
           if (cc .GT. 1.0D+00) cc = 1.0D+00
           if (cc .LT. 0.0D+00) cc = 0.0D+00
           intgval2(i0,j0) = cc
         enddo
         enddo
*       else if (ienhintg .EQ. 3) then ! Max/Min of Entire image
       else if ((ienhintg .EQ. 3) .OR. (ienhintg .EQ. 6)) then ! Max/Min of Entire image
         if (dabs(maxig - minig) .GT. 1.0D-10) then
           do j0 = 1, jvpix2
           do i0 = 1, ihpix2
             cc = (intgval(i0,j0) - minig) / (maxig - minig)
             if (cc .GT. 1.0D+00) cc = 1.0D+00
             if (cc .LT. 0.0D+00) cc = 0.0D+00
             intgval2(i0,j0) = cc
           enddo
           enddo
         else
           do j0 = 1, jvpix2
           do i0 = 1, ihpix2
             intgval2(i0,j0) = 0.0D+00
           enddo
           enddo
         endif
       else if (ienhintg .EQ. 4) then ! Dev-0 : force ave=0
         avei2 = 0.0D+00
         numig = 0
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           cc = intgval(i0,j0)
           if (dabs(cc) .GT. 1.0D-10) then
             numig = numig + 1
             avei2 = avei2 + cc**2
           endif
         enddo
         enddo
         if (numig .GT. 0) then
           avei2 = avei2 / dfloat(numig)
           devig = dsqrt(avei2)
         else
           devig = 1.0D+00
         endif
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           cc = intgval(i0,j0) / devig / 3.0D+00
           if (cc .GT. 1.0D+00) cc = 1.0D+00
           if (cc .LT. 0.0D+00) cc = 0.0D+00
           intgval2(i0,j0) = cc
         enddo
         enddo
       else if (ienhintg .EQ. 5) then ! Histgram Hosei
         do i = 0, 255
           numhist(i) = 0
         enddo
         numig = 0
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           cc = intgval(i0,j0)
           if (dabs(cc) .GT. 1.0D-10) then
             numig = numig + 1
             cc = (cc - minig) / (maxig - minig) * 255.0
             idummy1 = int(cc)
             numhist(idummy1) = numhist(idummy1) + 1
           endif
         enddo
         enddo
         do i = 0, 255
           sumhist(i) = 0
           do j = 0, i
             sumhist(i) = sumhist(i) + numhist(j)
           enddo
         enddo
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           aa = dfloat(i0-ic)
           bb = dfloat(j0-jc)
           cc = dsqrt(aa**2 + bb**2)
           cc = intgval(i0,j0)
           if (dabs(cc) .GT. 1.0D-10) then
             cc = (cc - minig) / (maxig - minig) * 255.0
             idummy1 = int(cc)
             intgval2(i0,j0) = dfloat(sumhist(idummy1))
     &                       / dfloat(sumhist(255)) ! dfloat(numig)
           else
             intgval2(i0,j0) = 0.0D+00
           endif
         enddo
         enddo
       else if (ienhintg .EQ. 7) then ! Newkirk filter + Max/Min of Entire image
         maxig2 = -1.0D+30
         minig2 =  1.0D+30
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           aa = (dfloat(i0-ic))**2 + (dfloat(j0-jc))**2
           aa = dsqrt(aa) * ul / rbottom ! in Rs...only for coronal simulation.
*           if ((mod(j0,100) .EQ. 0) .AND. (mod(i0,100) .EQ. 0))
*     +                    write(*,*) i0,j0,aa
           aa = (aa - 1.0) / 0.2 + 1.0
           cc = dmod(aa,1.0D+00)
           idummy1 = int(aa)
           if (idummy1 .LT.  0) then
             bb = fnewkirk(0)
           else if (idummy1 .GT. 14) then
             bb = fnewkirk(15)
           else
             bb = fnewkirk(idummy1) * (1.0D+00 - cc)
     &          + fnewkirk(idummy1 + 1) * cc
           endif
*           if ((mod(j0,100) .EQ. 0) .AND. (mod(i0,100) .EQ. 0))
*     +                    write(*,*) aa,cc,bb
           bb = (10.0D+00)**bb
           cc =  intgval(i0,j0) / bb
           intgval2(i0,j0) = cc
           if (cc .GT. maxig2) maxig2 = cc
           if (cc .LT. minig2) minig2 = cc
         enddo
         enddo
         if (dabs(maxig2 - minig2) .GT. 1.0D-10) then
           do j0 = 1, jvpix2
           do i0 = 1, ihpix2
*             cc = (intgval2(i0,j0) - minig2) / (maxig2 - minig2) ! relative brightness
             cc = intgval2(i0,j0) / maxig2 ! abs. brightness
             if (cc .GT. 1.0D+00) cc = 1.0D+00 ! for safety... not needed
             if (cc .LT. 0.0D+00) cc = 0.0D+00 ! for safety
             intgval2(i0,j0) = cc
           enddo
           enddo
         else
           do j0 = 1, jvpix2
           do i0 = 1, ihpix2
             intgval2(i0,j0) = 0.0D+00
           enddo
           enddo
         endif
       else ! no - emphasizing...
         do j0 = 1, jvpix2
         do i0 = 1, ihpix2
           intgval2(i0,j0) = intgval(i0,j0)
         enddo
         enddo
       endif

       do 888 i = istart, iend
       do 888 j = jstart, jend
* Get the point on the projection plane at which line of sight starts
         aa = ul * dfloat(i-ic)
         bb = ul * dfloat(j-jc)
         do m = 1, 3
           xyz0(m) = xyzcp(m) + aa * xyzxp(m) + bb * xyzyp(m)
         enddo
         do m = 1, 3
           lofsxyz(m)= xyz0(m) - xyzeye(m)
         enddo
         call normvec3(lofsxyz, tt0)
         aa = - naiseki3(lofsxyz,xyzeye)
         do m = 1, 3
           xyzdummy(m) = xyzeye(m) + aa * lofsxyz(m)
         enddo
         rrnear = absvec3(xyzdummy)
         if ((rrnear .LT. rupper) .AND. (rrnear .GT. rbottom2)) then

* Rupper for los integration
          if (rrnear .GT. rbottom) then ! corona graph. line

           fff = intgval2(i,j)
           if (fff .GT. 1.0D+00) fff = 1.0D+00
           if (fff .LT. 1.0D-10) fff = 1.0D-10
           fff = fff**0.5
           if (fff .LT. 0.5D+00) then
             red = 0
             gre = 0
             blu = int(fff * 255.0)
           else
             fff = (fff - 0.5D+00) * 2.0D+00
             red = int(fff * 255.0)
             gre = int(fff * 255.0)
             blu = int(fff * 127.0) + 128
           endif

          else !                          EIT like ....

           fff = intgval2(i,j) * 30.0 ! adjust ....
           if (fff .GT. 1.0D+00) fff = 1.0D+00
           if (fff .LT. 1.0D-10) fff = 1.0D-10
           fff = fff**0.5
*           fff = fff**2
           if (fff .LT. 0.5D+00) then
             red = int(fff * 255.0)
             gre = int(fff * 255.0)
             blu = 0
           else
             fff = (fff - 0.5D+00) * 2.0D+00
             red = int(fff * 127.0) + 128
             gre = int(fff * 127.0) + 128
             blu = int(fff * 255.0)
           endif

          endif

          idummy1 =  0 ! ichar(rgb2(1,i,j))
          aa = dfloat(idummy1**2 + red**2)
          aa = dsqrt(aa)
          idummy1 = int(aa)
          idummy2 = 0 ! ichar(rgb2(2,i,j))
          aa = dfloat(idummy2**2 + gre**2)
          aa = dsqrt(aa)
          idummy2 = int(aa)
          idummy3 = 0 ! ichar(rgb2(3,i,j))
          aa = dfloat(idummy3**2 + blu**2)
          aa = dsqrt(aa)
          idummy3 = int(aa)
          if (idummy1 .GT. 255) idummy1 = 255
          if (idummy2 .GT. 255) idummy2 = 255
          if (idummy3 .GT. 255) idummy3 = 255
          if (idummy1 .LT.   0) idummy1 =   0
          if (idummy2 .LT.   0) idummy2 =   0
          if (idummy3 .LT.   0) idummy3 =   0
          rgb2(1,i,j) = char(idummy1)
          rgb2(2,i,j) = char(idummy2)
          rgb2(3,i,j) = char(idummy3)
 
*         else
*           red = 0
*           gre = 0
*           blu = 0
         endif

  888  continue

       return
       end


*
***=======================================================================
*
* ivar = 1 : HCS (Br = 0)
*      = 2 : div|v| / |v| or div(rho) / rho
*
*
***=======================================================================
*
* ivar = 1 : HCS (Br = 0)
*      = 2 : div|v| / |v| or div(rho) / rho
*
       subroutine sclr2rgb(ivar)
       implicit  none
* interface
       integer ivar
* local
       real*8  dstnt2
       integer istart, iend, jstart, jend
       integer ic, jc          ! pixel potition of center of project plane
       real*8  aaa, bbb, ccc, fff, ddd
       real*8  aa, bb, cc, tt, tt0
       integer clrbtm, clrmdl  ! minimum and miduim color value of current sheet
       real*8  zza, zzb
       integer ifind
       integer i, j, m
       integer red, gre, blu
       real*8  xyza(3), xyz0(3), xyzeye(3), nprj(3)
       real*8  nrm(3)
       real*8  lofsxyz(3), mirxyz(3), obj2lght(3)
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       parameter(clrbtm = 100, clrmdl = 175) ! color table control index
* dummy
       integer idummy1, idummy2, idummy3
       real*8  rdummy1, rdummy2
       real*8  xyzdummy(3)
       logical ldummy
* functions
       real*8  naiseki3, absvec3
       logical sameside
* common
       integer iout, irlog
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom, gamma
       real*8  ro0(0:ii+1,0:jj,0:kk+1)
       real*8  pg0(0:ii+1,0:jj,0:kk+1)
       real*8  ur0(0:ii+1,0:jj,0:kk+1)
       real*8  ut0(0:ii+1,0:jj,0:kk+1)
       real*8  up0(0:ii+1,0:jj,0:kk+1)
       real*8  br0(0:ii+1,0:jj,0:kk+1)
       real*8  bt0(0:ii+1,0:jj,0:kk+1)
       real*8  bp0(0:ii+1,0:jj,0:kk+1)
       common  /var01/  ro0, pg0
       common  /var02/  ur0, ut0, up0
       common  /var04/  br0, bt0, bp0
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /cgcfg/  iout, irlog
       common  /coord/  rr, theta, phi
       common  /sclscl/ sclr
       common  /images/ dstnt, rgb2
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom
       common  /gasgas/ gamma

* pixel position on the projection plane
       istart = 1 ! (ihpix2 - jvpix2) / 2 + 1
       iend   = ihpix2 ! istart + jvpix2
       jstart = 1
       jend   = jvpix2
       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,aa)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo

       call setsclr(ivar,idummy1,idummy2) ! 2nd & 3rd var means nothig here
*       if ((ivar .EQ. 1) .OR. (ivar .EQ. 2)) call smthsclr(sclr)

** Get the RGB value of each pixel
       do 777 i = istart, iend
       do 777 j = jstart, jend

* Get the point on the projection plane at which line of sight starts
         aa = ul * dfloat(i-ic)
         bb = ul * dfloat(j-jc)
         do m = 1, 3
           xyz0(m) = xyzcp(m) + aa * xyzxp(m) + bb * xyzyp(m)
         enddo
* unit vector of line of sight
         do m = 1, 3
           lofsxyz(m)= xyz0(m) - xyzeye(m)
         enddo
         call normvec3(lofsxyz, tt0)

* Calculation will NOT be performed when
*   line of sight runs always beyond the sphere of radius of R_upper
         aa = - naiseki3(lofsxyz,xyzeye)
         do m = 1, 3
           xyzdummy(m) = xyzeye(m) + aa * lofsxyz(m)
         enddo
         aa = absvec3(xyzdummy)

         if (aa .LT. rupper) then

* Find starting point of path where line of sight cross the outer sphere.
           bb = naiseki3(xyzeye,lofsxyz)
           cc = bb**2 - absvec3(xyzeye) * absvec3(xyzeye) + rupper**2
           tt = - bb - dsqrt(cc)
           if (tt .LT. tt0) tt = tt0
           tt = tt * 1.001
           do m = 1, 3
             xyza(m) = xyzeye(m) + tt * lofsxyz(m)
           enddo

           call findcros(xyza,lofsxyz,zzb,zza,ifind)

* if tracing point out of range
           if (ifind .EQ. 0)
     &       write(11,'('' Ifind = 0 at : '',2i4)') i, j

* when found the point .......................

* check the object locate the other side than point of eye
*                     with respect to the projection plane

           ldummy = .true.
           if  (ifind .GE. 1)
     &       ldummy = sameside(xyzcp,nprj,xyza,xyzeye)

           if ((ifind .GE. 1) .AND. (.NOT. ldummy)) then

             call prjpnt1(xyza,rdummy1,rdummy2,dstnt2)
*             call gtdivsf2(xyza,nrm)
             call gtdivsuf(xyza,nrm)

* light direction
             do m = 1, 3
               obj2lght(m) = xyzlght(m) - xyza(m)
             enddo
             call normvec3(obj2lght,rdummy1)

* correct the normal if needed.
             rdummy2 = - naiseki3(nrm,lofsxyz)
             if (rdummy2 .LT. 0.0D+00) then
               do m = 1, 3
                 nrm(m) = - nrm(m)
               enddo
             endif

* color table
             aaa = naiseki3(obj2lght, nrm)
             if (aaa .GT. 0.0D+00) then
               mirxyz(1) = 2.0D+00 * aaa * nrm(1) - obj2lght(1) ! |mirxx| = 1
               mirxyz(2) = 2.0D+00 * aaa * nrm(2) - obj2lght(2)
               mirxyz(3) = 2.0D+00 * aaa * nrm(3) - obj2lght(3)
               ccc = - naiseki3(lofsxyz,mirxyz)
               if (ccc .LT. 0.0D+00) ccc = 0.0D+00
             else
               ccc = 0.0D+00
               aaa = 0.0D+00
*               aaa = dabs(aaa) / 2.0D+00
             endif

             bbb = - naiseki3(lofsxyz,nrm)
             bbb = dabs(bbb)

             aaa = 3.0D+00 * aaa**2 / (2.0D+00 + aaa**2)
             bbb = 3.0D+00 * bbb**2 / (2.0D+00 + bbb**2)
             ccc = 3.0D+00 * ccc**2 / (2.0D+00 + ccc**2)

             fff = (1.0D+00
     &            + 2.5D+00 * aaa
     &            + 2.5D+00 * ccc
     &            + 4.0D+00 * bbb) / 10.0D+00

*             fff = 0.15D+00
*     +           + 0.85D+00
*     +           *(2.00D+00 * bbb**2 / (1.0D+00 + bbb**2))**2
*     +           *(aaa**2 + ccc**2 + 0.2D+00) / 2.0D+00

*             fff =((1.0 - bbb**2) / (1.0 + bbb**2))**5

             if (fff .GT. 1.0D+00) fff = 1.0D+00
             if (fff .LT. 0.0D+00) fff = 0.0D+00

* add color to RGB data set
* The IVAR=1 (HCS) case should done before IVAR=2 case and TUBE
* The IVAR=2 case should be done at the last process
             if (ivar .EQ. 1) then ! non-transparent
               bbb = zzb - zza
               if (bbb .GE. 0.0D+00) then
                 if (fff .LT. 0.5D+00) then
                   aaa = 256.0D+00 * 2.0D+00 - clrbtm * 2.0D+00
                   red = 0
                   gre = 0
                   blu = int(fff * aaa) + clrbtm
                 endif
                 if (fff .GE. 0.5D+00) then
                   red = int(fff * 255.0D+00) - 127
                   gre = int(fff * 255.0D+00) - 127
                   blu = 255
                 endif
               else
                 if (fff .LT. 0.5D+00) then
                   red = int(fff * 480.0D+00) + 15
                   gre = 0
                   blu = 0
                 endif
                 if (fff .GE. 0.5D+00) then
                   red = 255
                   gre = int(fff * 255.0D+00) - 127
                   blu = int(fff * 255.0D+00) - 127
                 endif
               endif
               idummy1 = red
               idummy2 = gre
               idummy3 = blu
               if (idummy1 .GT. 255) idummy1 = 255
               if (idummy2 .GT. 255) idummy2 = 255
               if (idummy3 .GT. 255) idummy3 = 255
               if (idummy1 .LT.   0) idummy1 =   0
               if (idummy2 .LT.   0) idummy2 =   0
               if (idummy3 .LT.   0) idummy3 =   0
               rgb2(1,i,j) = char(idummy1)
               rgb2(2,i,j) = char(idummy2)
               rgb2(3,i,j) = char(idummy3)
               dstnt(i,j) = dstnt2
             endif

             if ((ivar .NE. 1) .AND. (dstnt2 .LT. dstnt(i,j))) then
               ddd = 7.0D-01 ! semi-transparent
               if (fff .LT. 0.5D+00) then
                 aaa = 256.0D+00 * 2.0D+00 - clrbtm * 2.0D+00
                 red = int(fff * aaa) + clrbtm
                 gre = 0
                 blu = int(fff * aaa) + clrbtm
               else
                 red = 255
                 gre = int(fff * 255.0D+00) - 127
                 blu = 255
               endif
               red = int(fff * 150) + 105 ! test
               gre = int(fff * 150) + 105 ! test
               blu = int(fff * 150) + 105 ! test
               idummy1 = ichar(rgb2(1,i,j))
               aa = ddd*dfloat(red) + (1.0D+00-ddd)*dfloat(idummy1)
               idummy1 = int(aa)
               idummy2 = ichar(rgb2(2,i,j))
               aa = ddd*dfloat(gre) + (1.0D+00-ddd)*dfloat(idummy2)
               idummy2 = int(aa)
               idummy3 = ichar(rgb2(3,i,j))
               aa = ddd*dfloat(blu) + (1.0D+00-ddd)*dfloat(idummy3)
               idummy3 = int(aa)
               if (idummy1 .GT. 255) idummy1 = 255
               if (idummy2 .GT. 255) idummy2 = 255
               if (idummy3 .GT. 255) idummy3 = 255
               if (idummy1 .LT.   0) idummy1 =   0
               if (idummy2 .LT.   0) idummy2 =   0
               if (idummy3 .LT.   0) idummy3 =   0
               rgb2(1,i,j) = char(idummy1)
               rgb2(2,i,j) = char(idummy2)
               rgb2(3,i,j) = char(idummy3)
             endif

           endif ! End of the process of determining color

         endif

 777   continue

       return
       end

* -----------
* scalar surface with shadow..
*
       subroutine scl2rgb2(ivar)
       implicit  none
* interface
       integer ivar
* local
       real*8  dstnt2
       integer istart, iend, jstart, jend
       integer ic, jc          ! pixel potition of center of project plane
       real*8  aaa, bbb, ccc, fff ! , ddd
       real*8  aa, bb, cc, tt, tt0
       integer clrbtm, clrmdl  ! minimum and miduim color value of current sheet
       real*8  zza, zzb
       integer ifind
       integer i, j, m
       integer red, gre, blu
       real*8  xyza(3), xyz0(3), xyzeye(3), nprj(3)
       real*8  nrm(3)
       real*8  lofsxyz(3), mirxyz(3), obj2lght(3)
       real*8  drct2(3)
       real*8  xyza2(3)
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       parameter(clrbtm = 100, clrmdl = 175) ! color table control index
* dummy
       integer idummy1, idummy2, idummy3
       real*8  rdummy1, rdummy2
       real*8  xyzdummy(3)
       logical ldummy, lshadow2
* functions
       real*8  naiseki3, absvec3
       logical sameside
       logical lshadow
* common
       integer iout, irlog
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom, gamma
       real*8  ro0(0:ii+1,0:jj,0:kk+1)
       real*8  pg0(0:ii+1,0:jj,0:kk+1)
       real*8  ur0(0:ii+1,0:jj,0:kk+1)
       real*8  ut0(0:ii+1,0:jj,0:kk+1)
       real*8  up0(0:ii+1,0:jj,0:kk+1)
       real*8  br0(0:ii+1,0:jj,0:kk+1)
       real*8  bt0(0:ii+1,0:jj,0:kk+1)
       real*8  bp0(0:ii+1,0:jj,0:kk+1)
       common  /var01/  ro0, pg0
       common  /var02/  ur0, ut0, up0
       common  /var04/  br0, bt0, bp0
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /cgcfg/  iout, irlog
       common  /coord/  rr, theta, phi
       common  /sclscl/ sclr
       common  /images/ dstnt, rgb2
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom
       common  /gasgas/ gamma

* pixel position on the projection plane
       istart = 1 ! (ihpix2 - jvpix2) / 2 + 1
       iend   = ihpix2 ! istart + jvpix2
       jstart = 1
       jend   = jvpix2
       ic     = ihpix2 / 2
       jc     = jvpix2 / 2

* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,aa)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo

       call setsclr(ivar,idummy1,idummy2) ! 2nd & 3rd var means nothig here
*       if ((ivar .EQ. 1) .OR. (ivar .EQ. 2)) call smthsclr(sclr)

** Get the RGB value of each pixel
       do 777 i = istart, iend
       do 777 j = jstart, jend

* Get the point on the projection plane at which line of sight starts
         aa = ul * dfloat(i-ic)
         bb = ul * dfloat(j-jc)
         do m = 1, 3
           xyz0(m) = xyzcp(m) + aa * xyzxp(m) + bb * xyzyp(m)
         enddo
* unit vector of line of sight
         do m = 1, 3
           lofsxyz(m)= xyz0(m) - xyzeye(m)
         enddo
         call normvec3(lofsxyz, tt0)

* Calculation will NOT be performed when
*   line of sight runs always beyond the sphere of radius of R_upper
         aa = - naiseki3(lofsxyz,xyzeye)
         do m = 1, 3
           xyzdummy(m) = xyzeye(m) + aa * lofsxyz(m)
         enddo
         aa = absvec3(xyzdummy)

         if (aa .LT. rupper) then

* Find starting point of path where line of sight cross the outer sphere.
           bb = naiseki3(xyzeye,lofsxyz)
           cc = bb**2 - absvec3(xyzeye) * absvec3(xyzeye) + rupper**2
           tt = - bb - dsqrt(cc)
           if (tt .LT. tt0) tt = tt0
           tt = tt * 1.001
           do m = 1, 3
             xyza(m) = xyzeye(m) + tt * lofsxyz(m)
           enddo

           call findcros(xyza,lofsxyz,zzb,zza,ifind)

* if tracing point out of range
           if (ifind .EQ. 0)
     &       write(11,'('' Ifind = 0 at : '',2i4)') i, j

* when found the point .......................

* check the object locate the other side than point of eye
*                     with respect to the projection plane

           ldummy = .true.
           if  (ifind .GE. 1)
     &       ldummy = sameside(xyzcp,nprj,xyza,xyzeye)

           if ((ifind .GE. 1) .AND. (.NOT. ldummy)) then

             call prjpnt1(xyza,rdummy1,rdummy2,dstnt2)
*             call gtdivsf2(xyza,nrm)
             call gtdivsuf(xyza,nrm)

* light direction
             do m = 1, 3
               obj2lght(m) = xyzlght(m) - xyza(m)
             enddo
             call normvec3(obj2lght,rdummy1)

             drct2(1) = obj2lght(1)
             drct2(2) = obj2lght(2)
             drct2(3) = obj2lght(3)

             aaa = naiseki3(obj2lght,nrm)
             if (aaa .LT. -0.01D+00) then
               xyza2(1) = xyza(1)
               xyza2(2) = xyza(2)
               xyza2(3) = xyza(3)
               lshadow2 = lshadow(xyza2,drct2)
             else
               lshadow2 = .true.
             endif


* correct the normal if needed.
*             rdummy2 = - naiseki3(nrm,lofsxyz)
*             if (rdummy2 .LT. 0.0D+00) then
*               do m = 1, 3
*                 nrm(m) = - nrm(m)
*               enddo
*             endif

* color table
             if (lshadow2) then
               aaa = 0.0D+00
               ccc = 0.0D+00
             else
               mirxyz(1) = 2.0D+00 * aaa * nrm(1) - obj2lght(1) ! |mirxx| = 1
               mirxyz(2) = 2.0D+00 * aaa * nrm(2) - obj2lght(2)
               mirxyz(3) = 2.0D+00 * aaa * nrm(3) - obj2lght(3)
               call normvec3(mirxyz,rdummy1)
               ccc = - naiseki3(lofsxyz,mirxyz)
               if (ccc .LT. 0.0D+00) ccc = 0.0D+00
               ccc = dabs(ccc) + 0.000001
               ccc = 2.0 * ccc**4 / (ccc**4 + 1.0) ! Mirror
             endif

             bbb = naiseki3(lofsxyz,nrm)
             bbb = dabs(bbb)

             if (bbb .LT. 0.0D+00) bbb = 0.0D+00
             fff =(0.5D+00
     &           + 1.0D+00 * abs(bbb)
     &           + 1.0D+00 * abs(aaa)
     &           + 1.0D+00 * ccc)
     &           / 4.0D+00
             if (fff .GT. 1.0D+00) fff = 1.0D+00
             if (fff .LT. 0.0D+00) fff = 0.0D+00
             fff = abs(fff) * 255.0
*             fff = sqrt(fff) * 255.0

* The IVAR=1 (HCS) case should done before IVAR=2 case and TUBE
             bbb = zzb - zza
             if (bbb .GE. 0.0D+00) then ! blue !
               if (fff .LE. 127.0) then
                 red = 0
                 gre = 0
                 blu = int(fff * 1.5)
               else
                 red = int(2.000 * (fff - 127.0) * 0.9999)
                 gre = int(2.000 * (fff - 127.0) * 0.9999)
                 blu = int(0.666 * (fff - 127.0) * 0.9999 + 192.0)
               endif
             else ! red
               if (fff .LE. 127.0) then
                 blu = 0
                 gre = 0
                 red = int(fff * 1.5)
               else
                 blu = int(2.000 * (fff - 127.0) * 0.9999)
                 gre = int(2.000 * (fff - 127.0) * 0.9999)
                 red = int(0.666 * (fff - 127.0) * 0.9999 + 192.0)
               endif
             endif

             idummy1 = red
             idummy2 = gre
             idummy3 = blu
             if (idummy1 .GT. 255) idummy1 = 255
             if (idummy2 .GT. 255) idummy2 = 255
             if (idummy3 .GT. 255) idummy3 = 255
             if (idummy1 .LT.   0) idummy1 =   0
             if (idummy2 .LT.   0) idummy2 =   0
             if (idummy3 .LT.   0) idummy3 =   0
             rgb2(1,i,j) = char(idummy1)
             rgb2(2,i,j) = char(idummy2)
             rgb2(3,i,j) = char(idummy3)
             dstnt(i,j) = dstnt2

           endif ! End of the process of determining color

         endif

 777   continue

       return
       end

* -----------

*
* -----------------------------------------------------------
*
*       subroutine smthsclr(sclr)
*       implicit none
*       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
*       real*8  sclr(0:ii+1,0:jj,0:kk+1)
*       real*8  scl2(1:ii,1:jj-1,0:kk+1)
*       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
*       integer i, j, k
*       real*8  bb, cc
*       common  /coord/ rr, theta, phi

*       do 234 k  =  1, kk
*       do 234 j  =  1, jj - 1
*       do 234 i  =  1, ii
*         bb =(sclr(i,j-1,k)*dsin(theta(j-1))
*     +       +sclr(i,j  ,k)*dsin(theta(j  )) * 2.0D+00
*     +       +sclr(i,j+1,k)*dsin(theta(j+1)))
*     +      /(dsin(theta(j-1))
*     +       +dsin(theta(j  )) * 2.0D+00
*     +       +dsin(theta(j+1)))
*         cc =(sclr(i,j,k-1)
*     +       +sclr(i,j,  k)*2.0D+00
*     +       +sclr(i,j,k+1))/4.0D+00
*         scl2(i,j,k) = (bb + cc) / 2.0D+00
**         scl2(i,j,k) = scl2(i,j,k) + (bb + cc) / 3.0D+00
* 234   continue

*       do 240 j = 1, jj - 1
*       do 240 i = 1, ii
*         scl2(i,j,   0) = scl2(i,j,kk)
*         scl2(i,j,kk+1) = scl2(i,j, 1)
* 240   continue

*       do 250 k = 0, kk + 1
*       do 250 j = 1, jj - 1
*       do 250 i = 1, ii
*         sclr(i,j,k) = scl2(i,j,k)
* 250   continue
*       return
*       end

*
* -----------------------------------------------------------
*
       subroutine gtdivsuf(xyza,nrm)
       implicit none
       real*8  xyza(3), nrm(3)
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  rtp(3)
       real*8  ra2(0:6), th2(0:6), ph2(0:6), kin(1:6)
       integer ijk(3)
       integer i1, i2
       integer l
       real*8  ra0, th0, ph0
       real*8  dra0, dth0, dph0
       real*8  aaa
       real*8  sclatrtp
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       common  /coord/  rr, theta, phi
       common  /sclscl/ sclr

       call xyz2rtp(xyza,rtp)
       ra0 = rtp(1)
       th0 = rtp(2)
       ph0 = rtp(3)
       ra2(0) = ra0
       th2(0) = th0
       ph2(0) = ph0
* set coordinate of 6 points in the visinity of the center point
       call rtp2ijk(rtp,ijk)
       i1 = ijk(1)
       i2 = ijk(1) + 1
       dra0 = (rr(i2) - rr(i1))
       dth0 =  pi / dfloat(jj)
       dph0 = 2.0D+00 * pi / dfloat(kk)
       do l = 1, 6
         ra2(l) = ra2(0)
         th2(l) = th2(0)
         ph2(l) = ph2(0)
       enddo
       ra2(1) = ra2(0) + dra0
       ra2(2) = ra2(0) - dra0
       th2(3) = th2(0) + dth0
       th2(4) = th2(0) - dth0
       ph2(5) = ph2(0) + dph0
       ph2(6) = ph2(0) - dph0
       if (ph2(6) .LT. 0.0D+00   ) ph2(6) = ph2(6) + 2.0D+00 * pi
       if (ph2(6) .GT. 2.0D+00*pi) ph2(6) = ph2(6) - 2.0D+00 * pi
* get properties of 6 points in the visinity of the center point
       do l = 1, 6
         rtp(1) = ra2(l)
         rtp(2) = th2(l)
         rtp(3) = ph2(l)
         kin(l) = sclatrtp(rtp)
       enddo
* Get divergence
       nrm(1) = (kin(1)-kin(2))/(2.0D+00*dra0)
       nrm(2) = (kin(3)-kin(4))/(2.0D+00*dth0*ra0)
       nrm(3) = (kin(5)-kin(6))/(2.0D+00*dph0*ra0*dsin(th0))
       rtp(1) = ra2(0)
       rtp(2) = th2(0)
       rtp(3) = ph2(0)
       call vrtp2xyz(nrm,rtp)
       call normvec3(nrm,aaa)

       return
       end


*
* -----------------------------------------------------------
*
       subroutine gtdivsf2(xyza,nrm)
       implicit none
*
       real*8  xyza(3), nrm(3)
*
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
*
       real*8  xx2(0:6), yy2(0:6), zz2(0:6), kin(1:6)
       integer l
       real*8  dx0, dy0, dz0
       real*8  rtp2(3), xyz2(3) ! , rtpa(3)
       real*8  divxxx, divyyy, divzzz
       real*8  aaa
       real*8  sclatrtp
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
*
       common  /coord/  rr, theta, phi
       common  /sclscl/ sclr

* set coordinate of 6 points in the visinity of the center point
       aaa = xyza(1)**2 + xyza(2)**2 + xyza(3)**2
       aaa = dsqrt(aaa) / 1.0D+04
       dx0 = aaa
       dy0 = aaa
       dz0 = aaa
       do l = 0, 6
         xx2(l) = xyza(1)
         yy2(l) = xyza(2)
         zz2(l) = xyza(3)
       enddo
       xx2(1) = xx2(0) + dx0
       xx2(2) = xx2(0) - dx0
       yy2(3) = yy2(0) + dy0
       yy2(4) = yy2(0) - dy0
       zz2(5) = zz2(0) + dz0
       zz2(6) = zz2(0) - dz0
* get properties of 6 points in the visinity of the center point
       do l = 1, 6
         xyz2(1) = xx2(l)
         xyz2(2) = yy2(l)
         xyz2(3) = zz2(l)
         call xyz2rtp(xyz2,rtp2)
         kin(l) = sclatrtp(rtp2)
       enddo
* Get divergence
       divxxx=(kin(1)-kin(2))/(2.0D+00*dx0)
       divyyy=(kin(3)-kin(4))/(2.0D+00*dy0)
       divzzz=(kin(5)-kin(6))/(2.0D+00*dz0)
       aaa = dsqrt(divxxx**2 + divyyy**2 + divzzz**2)
       nrm(1) = divxxx / aaa
       nrm(2) = divyyy / aaa
       nrm(3) = divzzz / aaa

       return
       end


*** =================================================================
* ifind = -1         trancing point out of domain of data
* ifind =  positive  find !!
* ifind =  0         something is wrong !???
*** =================================================================
*
       subroutine findcros(xyza,drct,zzb,zza,ifind)
       implicit none
* interface
       real*8  xyza(3), drct(3), zzb, zza
       integer ifind
* local
       integer sumtime
       integer i1, i2, j1, j2, k1, k2
       integer m, istop1
       integer ijk(3)
       real*8  dpi, dpi2
       real*8  rtp(3), xyzb(3)
       real*8  sclatrtp
       real*8  ra2, th2, ph2
       real*8  dss, dss1, dss2
       real*8  eee, aaa, bbb, ccc
       real*8  crtrn, pi
       parameter(crtrn = 0.0D+00)
       parameter(pi = 3.14159265358979D+00)
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  rupper, rbottom
       common  /coord/  rr, theta, phi
       common  /sclscl/ sclr
       common  /varfm1/ rupper, rbottom

* the tolerance at the singlar : rotational axis (theta = 0, pi)
       dpi  = pi / dfloat(jj)
       dpi2 = pi / dfloat(jj) * 1.2D+00 ! modify the factor

       istop1 = -1 ! off

       sumtime = -1
       ifind = 0
       zzb = 0.0D+00 ! dummy
 111   continue
         sumtime = sumtime + 1
         call xyz2rtp(xyza,rtp)
         ra2 = rtp(1)
         th2 = rtp(2)
         ph2 = rtp(3)
         if ((ra2 .GT. rupper) .OR.
     &       (ra2 .LT. rbottom) .OR. (sumtime .GE. 5000)) then
           ifind = -1
           goto 444
         endif
         if ((th2 .LT. dpi) .OR. (th2 .GT. (pi-dpi))) goto 333
* obtain the corresponding interger 'i j k' to the point 'a'
         zza =  sclatrtp(rtp)

* check if the value cross the criterion
         eee = (zzb - crtrn)*(zza - crtrn)  ! if first step then skip
         if ((sumtime .GT. 0) .AND. (eee .LT. 0.0D+00) .AND.
     &       (th2 .GE. dpi2)  .AND. (th2 .LE. (pi - dpi2))) then
           ifind = sumtime                  ! set frag of 'find' ON
           do m = 1, 3                      ! Find then get location(x-y-z)
             xyza(m) = xyzb(m)
     &               +(xyza(m)-xyzb(m))*(crtrn-zzb)/(zza-zzb)
           enddo
           goto 444
         endif

* copy properties for the next step (if not find yet) : 'a' --> 'b'
 333     continue
         zzb = zza
* set advancing step
         call rtp2ijk(rtp,ijk)
         i1 = ijk(1)
         i2 = ijk(1) + 1
         j1 = ijk(2)
         j2 = ijk(2) + 1
         k1 = ijk(3)
         k2 = ijk(3) + 1
*         dss = rtp(1) * dsin(rtp(2)) * pi / dfloat(kk) / 3.0D+00
         dss1 = rtp(1) * pi / dfloat(kk) / 5.0D+00
         dss2 = (rr(i2) - rr(i1)) / 2.0D+00
*         dss = dmin1(dss1,dss2) !                    slow .....
*         dss = dmin1(dss,dss1,dss2) !               very slow....
         dss = rtp(1) * pi / dfloat(kk) / 5.0D+00 ! faster but ...
         do m = 1, 3
           xyzb(m) = xyza(m)
           xyza(m) = xyza(m) + dss * drct(m)  ! newstep
         enddo
         aaa = dsqrt(xyza(1)**2 + xyza(2)**2 + xyza(3)**2)
         bbb = dsqrt(xyzb(1)**2 + xyzb(2)**2 + xyzb(3)**2)
         if ((aaa .GT. rupper) .AND. (istop1 .EQ. -1)) then
           istop1 = 1
           ccc = (rupper - bbb) / (aaa - bbb) * 0.99999D+00
           do m = 1, 3
             xyza(m) = xyzb(m) + dss * drct(m) * ccc
           enddo
         endif
       goto 111    ! try again

 444   continue

       return
       end

*
* -------------------------
*

       logical function lshadow(xyza,drct)
       implicit none
* interface
       real*8  xyza(3), drct(3)
* local
       real*8  zza, zzb
       integer ifind
       integer sumtime
       integer i1, i2, j1, j2, k1, k2
       integer m, istop1
       integer ijk(3)
       real*8  dpi, dpi2
       real*8  rtp(3), xyzb(3)
       real*8  sclatrtp
       real*8  ra2, th2, ph2
       real*8  dss, dss1, dss2
       real*8  eee, aaa, bbb, ccc
       real*8  crtrn, pi
       parameter(crtrn = 0.0D+00)
       parameter(pi = 3.14159265358979D+00)
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  rupper, rbottom
       common  /coord/  rr, theta, phi
       common  /sclscl/ sclr
       common  /varfm1/ rupper, rbottom

* the tolerance at the singlar : rotational axis (theta = 0, pi)
       dpi  = pi / dfloat(jj)
       dpi2 = pi / dfloat(jj) * 1.2D+00 ! modify the factor

       istop1 = -1 ! off

       sumtime = -1
       ifind = 0
       zzb = 0.0D+00 ! dummy
 111   continue
         sumtime = sumtime + 1
         call xyz2rtp(xyza,rtp)
         ra2 = rtp(1)
         th2 = rtp(2)
         ph2 = rtp(3)
         if ((ra2 .GT. rupper) .OR.
     &       (ra2 .LT. rbottom) .OR. (sumtime .GE. 5000)) then
           ifind = -1
           goto 444
         endif
         if ((th2 .LT. dpi) .OR. (th2 .GT. (pi-dpi))) goto 333
* obtain the corresponding interger 'i j k' to the point 'a'
         zza =  sclatrtp(rtp)

* check if the value cross the criterion
         eee = (zzb - crtrn)*(zza - crtrn)  ! if first step then skip
*         if ((sumtime .GT. 0) .AND. (eee .LT. 0.0D+00) .AND.
         if ((sumtime .GT. 1) .AND. (eee .LT. 0.0D+00) .AND.
     &       (th2 .GE. dpi2)  .AND. (th2 .LE. (pi - dpi2))) then
           ifind = sumtime                  ! set frag of 'find' ON
           do m = 1, 3                      ! Find then get location(x-y-z)
             xyza(m) = xyzb(m)
     &               +(xyza(m)-xyzb(m))*(crtrn-zzb)/(zza-zzb)
           enddo
           goto 444
         endif

* copy properties for the next step (if not find yet) : 'a' --> 'b'
 333     continue
         zzb = zza
* set advancing step
         call rtp2ijk(rtp,ijk)
         i1 = ijk(1)
         i2 = ijk(1) + 1
         j1 = ijk(2)
         j2 = ijk(2) + 1
         k1 = ijk(3)
         k2 = ijk(3) + 1
*         dss = rtp(1) * dsin(rtp(2)) * pi / dfloat(kk) / 3.0D+00
         dss1 = rtp(1) * pi / dfloat(kk) / 5.0D+00
         dss2 = (rr(i2) - rr(i1)) / 2.0D+00
         dss = dmin1(dss1,dss2)
         dss = rtp(1) * pi / dfloat(kk) / 5.0D+00
         do m = 1, 3
           xyzb(m) = xyza(m)
           xyza(m) = xyza(m) + dss * drct(m)  ! newstep
         enddo
         aaa = dsqrt(xyza(1)**2 + xyza(2)**2 + xyza(3)**2)
         bbb = dsqrt(xyzb(1)**2 + xyzb(2)**2 + xyzb(3)**2)
         if ((aaa .GT. rupper) .AND. (istop1 .EQ. -1)) then
           istop1 = 1
           ccc = (rupper - bbb) / (aaa - bbb) * 0.99999D+00
           do m = 1, 3
             xyza(m) = xyzb(m) + dss * drct(m) * ccc
           enddo
         endif
       goto 111    ! try again

 444   continue

       lshadow = (ifind .GE. 1)

       return
       endfunction


**=======================================================================
*
       subroutine ocsf2rgb()
       implicit none
* local
       integer i, j, k
       real*8  pi, dpi
       parameter(pi = 3.14159265358979D+00)
       real*8  ra0, th0, ph0
       real*8  aa
       integer linemax
       parameter(linemax = 10000)
       real*8  rtpline(3,0:linemax), avec(0:linemax)
       integer nline, nfoot2
       real*8  rupper2, rbottom2
       logical lfromsun
       integer iii, m
       integer ic, jc          ! pixel potition of center of project plane
       real*8  bbb, ccc
       real*8  xp1, yp1
       integer ix1, jy1
       integer red, gre, blu
       integer ijk(3)
       integer j1, k1, j2, k2
       real*8  rtp(3), xyz(3), obj2lght(3)
       real*8  xyzeye(3), nprj(3)
       real*8  surfbr, dstnt2
       real*8  dth1, dth2, dph1, dph2
       integer ivar
* dummy
       integer idummy1, idummy2, idummy3
       real*8  rdummy1
       logical ldummy
* functions
       real*8  naiseki3
       logical sameside
* common
       integer iout, irlog
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  oc1(0:jj,0:kk+1)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       integer nfoot
       parameter(nfoot = 10000) ! approxiamte..
       logical lopen1st(nfoot)
       common  /lopline/ lopen1st
       common  /sclscl/ sclr
       common  /cgcfg/  iout, irlog
       common  /coord/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /images/ dstnt, rgb2
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

       ivar = 1
* pixel position on the projection plane
       ic     = ihpix2 / 2
       jc     = jvpix2 / 2
* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,aa)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo
       dpi = pi / dfloat(jj)

       do 300 k = 0, kk + 1
       do 300 j = 0, jj
       do 300 i = 0, ii + 1
         vcr(i,j,k) = br1(i,j,k)
         vct(i,j,k) = bt1(i,j,k)
         vcp(i,j,k) = bp1(i,j,k)
 300   continue
       lfromsun = .true.

*       rupper2 = rupper
*       rupper2 = rbottom * 2.5D+00 ! MIND this may not be exactly 2.5 Rs
       rupper2 = rbottom * 5.0D+00 ! MIND this may not be exactly 2.5 Rs
       if (rupper2 .GT. rr(ii)) rupper2 = rr(ii)
       rbottom2 = rbottom
       if (rbottom2 .LT. rr(1)) rbottom2 = rr(1)

       ra0 = rupper2 / 100.0
       if (ra0 .LT. rbottom2) ra0 = rbottom2

       write(11,'('' R_0   = '',f9.3)') ra0
       write(11,'('' R_bot = '',f9.3)') rbottom2
       write(11,'('' R_up  = '',f9.3)') rupper2

       nfoot2 = 0 ! counter
       do j = 0, jj     ! latitude
         th0 = pi * dfloat(j) / dfloat(jj)
         do k = 1, kk   ! longitude

           if ((j .EQ. 0) .OR. (j .EQ. jj)) then
             oc1(j,k) = 0.0D+00 ! open
*             oc1(j,k) = 0.5D+00 ! half
           else
             nfoot2 = nfoot2 + 1
             ph0 = 2.00D+00 * pi * dfloat(k) / dfloat(kk)

             if (ph0 .LT. 0.0D+00) then
               ph0 = ph0 + 2.0D+00 * pi
             else if (ph0 .GT. 2.0D+00 * pi) then
               ph0 = dmod(ph0, 2.0D+00 * pi)
             endif

             rtp(1) = ra0
             rtp(2) = th0
             rtp(3) = ph0
             call rtp2xyz(rtp,xyz)

             call prjpnt1(xyz,xp1,yp1,dstnt2)
             call pxy2ij(xp1,yp1,ul,ix1,jy1,ic,jc)
             ldummy = sameside(xyzcp,nprj,xyz,xyzeye)
             if (.NOT. ldummy) then
               rtpline(1,0) = ra0
               rtpline(2,0) = th0
               rtpline(3,0) = ph0
               call trcline(rtpline,nline,rupper2,rbottom2,
     &                      ivar,avec,lfromsun)

               if (rtpline(1,nline) .GT. rupper2 * 0.99D+00)  then
                 oc1(j,k) = 0.0D+00 ! open !
               else if ((rtpline(2,nline) .LT. 0.05 * pi) .OR.
     &                  (rtpline(2,nline) .GT. 0.95 * pi)) then
                 oc1(j,k) = 0.0D+00 ! open
*                 oc1(j,k) = 0.5D+00 ! half
               else
                 oc1(j,k) = 1.0D+00 ! closed
               endif
             else
               oc1(j,k) = 1.0D+00 ! closed ... far side
             endif
           endif
         enddo
         oc1(j,   0) = oc1(j,kk)
         oc1(j,kk+1) = oc1(j, 1)
       enddo

* drawing the solar surface
       aa = rbottom / ul
       iii = int(aa * 10.0)
       if (iii .LE. 0) iii = 1
       do 111 k = 0, (iii * 2)
       do 111 j = 0, iii
         rtp(1) = rr(0) * 1.002 ! rbottom
         rtp(2) = pi * dfloat(j) / dfloat(iii)
         rtp(3) = pi * dfloat(k) / dfloat(iii)
         call rtp2xyz(rtp,xyz)

         call prjpnt1(xyz,xp1,yp1,dstnt2)
         call pxy2ij(xp1,yp1,ul,ix1,jy1,ic,jc)

         ldummy = sameside(xyzcp,nprj,xyz,xyzeye)

         if ((ix1 .GT. 0) .AND. (ix1 .LE. ihpix2) .AND.
     &       (jy1 .GT. 0) .AND. (jy1 .LE. jvpix2) .AND.
     &       (dstnt2 .LT. dstnt(ix1,jy1)) .AND. (.NOT. ldummy)) then

           dstnt(ix1,jy1) = dstnt2

* estimate the brightness of the solar surface : colar table
           call rtp2ijk(rtp,ijk)
           j1 = ijk(2)
           j2 = ijk(2) + 1
           k1 = ijk(3)
           k2 = ijk(3) + 1
           dth1 =(rtp(2) - theta(j1))/(theta(j2)-theta(j1))
           dth2 =(theta(j2) - rtp(2))/(theta(j2)-theta(j1))
           dph1 =(rtp(3) - phi(k1))/(phi(k2)-phi(k1))
           dph2 =(phi(k2) - rtp(3))/(phi(k2)-phi(k1))
           surfbr = (oc1(j1,k1)*dth2+oc1(j2,k1)*dth1)*dph2
     &             +(oc1(j1,k2)*dth2+oc1(j2,k2)*dth1)*dph1
           do m = 1, 3
             obj2lght(m) = xyzlght(m) - xyz(m)
           enddo
           call normvec3(obj2lght,rdummy1)
           call normvec3(xyz,rdummy1)
           bbb = naiseki3(obj2lght,xyz)
           ccc = 0.25D+00 * (bbb + 3.0D+00) * 255.0D+00
     &          *(surfbr + 1.0D+00) * 0.5D+00
           red = int(ccc)
           gre = int(ccc)
           blu = int(ccc)
* put color value to rgb-data-array
           idummy1 = red
           idummy2 = gre
           idummy3 = blu
           if (idummy1 .GT. 255) idummy1 = 255
           if (idummy2 .GT. 255) idummy2 = 255
           if (idummy3 .GT. 255) idummy3 = 255
           if (idummy1 .LT.   0) idummy1 =   0
           if (idummy2 .LT.   0) idummy2 =   0
           if (idummy3 .LT.   0) idummy3 =   0
           rgb2(1,ix1,jy1) = char(idummy1)
           rgb2(2,ix1,jy1) = char(idummy2)
           rgb2(3,ix1,jy1) = char(idummy3)
         endif
 111   continue

       return
       end

*
**=======================================================================
*
       subroutine tube2rgb(irtube,ivar,itubecol,v0,nframe)
       implicit  none
* interface
       integer irtube, ivar, itubecol, nframe
       real*8  v0
* local
       integer i, j, k, je, ke
       real*8  pi, dpi
       parameter(pi = 3.14159265358979D+00)
       real*8  ra0, th0, ph0
       real*8  rdummy, aa
       integer linemax
       parameter(linemax = 10000)
       real*8  rtpline(3,0:linemax), avec(0:linemax)
       integer nline, nfoot2
       real*8  rupper2, rbottom2
       logical lfromsun, lequator
       logical lopen1, lopen2
       integer ivar2, nrotref, nref
* functions
*       real*8  b0deg ! not used now
* common
       integer iout, irlog
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0 * ifine, jvpix2 = jvpix0 * ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       integer nfoot
       parameter(nfoot = 10000) ! approxiamte..
       logical lopen1st(nfoot)
       common  /lopline/ lopen1st
       common  /sclscl/ sclr
       common  /cgcfg/  iout, irlog
       common  /coord/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /images/ dstnt, rgb2
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

       dpi = pi / dfloat(jj)

       if (ivar .EQ. 1) then
         do 300 k = 0, kk + 1
         do 300 j = 0, jj
         do 300 i = 0, ii + 1
           vcr(i,j,k) = br1(i,j,k)
           vct(i,j,k) = bt1(i,j,k)
           vcp(i,j,k) = bp1(i,j,k)
 300     continue
       else  if (ivar .EQ. 2) then
         do 400 k = 0, kk + 1
         do 400 j = 0, jj
         do 400 i = 0, ii + 1
           vcr(i,j,k) = ur1(i,j,k)
           vct(i,j,k) = ut1(i,j,k)
           vcp(i,j,k) = up1(i,j,k)
 400     continue
       else
         write(*,*) ' Ivar is wrong '
         stop
       endif

*       rupper2 = rbottom * 5.0D+00 ! MIND this may not be exactly 2.5 Rs
       rupper2 = rupper
       if (rupper2 .GT. rr(ii)) rupper2 = rr(ii)
       rbottom2 = rbottom
       if (rbottom2 .LT. rr(1)) rbottom2 = rr(1)

       if (itubecol .GT. 10) then
         do k = 0, kk + 1
         do j = 0, jj
         do i = 0, ii + 1
           aa = ur1(i,j,k)**2 + ut1(i,j,k)**2 + up1(i,j,k)**2
           sclr(i,j,k) = dsqrt(aa+0.000000001) ! usual
*           sclr(i,j,k) = aa ! enhance ....
         enddo
         enddo
         enddo
       else !  otherwize normal...
         do k = 0, kk + 1
         do j = 0, jj
         do i = 0, ii + 1
           aa = vcr(i,j,k)**2 + vct(i,j,k)**2 + vcp(i,j,k)**2
           sclr(i,j,k) = dsqrt(aa+0.000000001) ! usual
*           sclr(i,j,k) = aa ! enhance
         enddo
         enddo
         enddo
       endif

       if ((itubecol .EQ. 2) .OR. (itubecol .EQ. 16)) then
         ra0 = rupper2 - 0.000001
         je = 1
         lfromsun = .false.
         lequator = .true.
       else if (itubecol .EQ. 4) then
         ra0 = rupper2 / 100.0
         if (ra0 .LT. rbottom2) ra0 = rbottom2
         je = 1
         lfromsun = .true.
         itubecol = 2
         lequator = .true.
       else
         ra0 = rupper2 / 100.0
         if (ra0 .LT. rbottom2) ra0 = rbottom2
*         je = jj / 10 ! approx. frequently changed
*         if (je .LT. 10) je = 10
         je = 20
*         je = 50
         lfromsun = .true.
         lequator = .false.
         if (itubecol .EQ.  3) itubecol = 2
       endif

       write(11,'('' R_0   = '',f9.3)') ra0
       write(11,'('' R_bot = '',f9.3)') rbottom2
       write(11,'('' R_up  = '',f9.3)') rupper2

       nfoot2 = 0 ! counter

       if (itubecol .EQ. 18) then
         ivar2 = 18 ! temp.... this happens to be same as itubecol.
         nrotref = 0
         nref   = 0 ! not used in setsclr with ivar2=18
         call setsclr(ivar2,nrotref,nref)
       endif

       do j = 1, je     ! latitude
*       do j = 2, je - 1    ! Skip near Axis
*       do j = 3, je - 2    ! Skip near Axis

         if (.NOT. lequator) then
           th0 = pi * (dfloat(j) - 5.0D-01)/ dfloat(je)
           rdummy = dfloat(2 * je) * dsin(th0) + 1.0
           ke = int(rdummy)
         else
           th0 = pi / 2.0D+00
           ke = 20
         endif

         do k = 1, ke   ! longitude
*         do k = ke / 2 - 2, ke / 2 + 2  ! longitude

           nfoot2 = nfoot2 + 1

           if (lfromsun) then
             ph0 = 2.00D+00 * pi * dfloat(k) / dfloat(ke)
     &           + 0.05D+00 * dfloat(k) ! shift
           else
             ph0 = 2.00D+00 * pi * (dfloat(k) + 5.0D-01) / dfloat(ke)
           endif

           if (ph0 .LT. 0.0D+00) then
             ph0 = ph0 + 2.0D+00 * pi
           else if (ph0 .GT. 2.0D+00 * pi) then
             ph0 = dmod(ph0, 2.0D+00 * pi)
           endif

           rtpline(1,0) = ra0
           rtpline(2,0) = th0
           rtpline(3,0) = ph0

           if (itubecol .EQ. 18) then
             call trcline(rtpline,nline,rupper2,rbottom2,
     &                     ivar,avec,lfromsun) ! linear coordinate
             do i = 0, nline
               avec(i) = avec(i) * 1.0 ! in unit MK
             enddo
           else
             call trcline(rtpline,nline,rupper2,rbottom2,
     &                    ivar,avec,lfromsun) ! linear coordinate
             do i = 0, nline
               avec(i) = avec(i) * v0 / 1.0D+05
             enddo
           endif

           if (itubecol .LT. 10) then
             lopen1 = .false. ! not used...
           else
             if (rtpline(1,nline) .GT. rupper2 * 0.99D+00)  then
               lopen1 = .true.
             else
               lopen1 = .false.
             endif
           endif
           if (nframe .EQ. 1)    lopen1st(nfoot2) = lopen1 ! save
           if (itubecol .EQ. 15) lopen2 = lopen1st(nfoot2) ! restore

           if (nline .GT. 2) ! normal... draw all lines
**
*           if ((nline .GT. 2) .AND. (.NOT. (lopen2))) ! draw only line closed in the first frame
*           if ((nline .GT. 2) .AND. (lopen1)) ! draw only line open to interplanetary space
**
*           if (nframe .EQ. 1) lopen2 = lopen1
*           if ((nline .GT. 2) .AND.
*     &         ((.NOT. lopen1) .OR.
*     &          (.NOT. lopen2))) ! draw only line closed in the first or current frame
     &       call circ2rgb(rtpline,nline,irtube,ivar,
     &                     itubecol,avec,lopen1, lopen2)

         enddo
       enddo

       return
       end


*
** -------------
*
       subroutine circ2rgb(rtpline,nline,irtube,
     &                     ivar,itubecol,avec,lopen1,lopen2)
       implicit none
* interface
       integer linemax
       parameter(linemax = 10000)
       real*8  rtpline(3,0:linemax)
       integer nline, irtube, ivar, itubecol
       real*8  avec(0:linemax)
       logical lopen1, lopen2
* local
       integer n, i1, i2, m
       real*8  ra, th, ph
* tube
       real*8  ncc(3,2), nseg(3,3)
       real*8  xyzcc(3), nrmcc(3)
       real*8  crimx0(3), crimy0(3), crim(3)
       real*8  xptub, yptub, ntube(3)
       real*8  xyzgaiin(3)
       real*8  lofsxyz(3), obj2lght(3), mirxyz(3), xyzeye(3), nprj(3)
       real*8  xyzl(3,4)
       real*8  rtub
       integer itheta0, itheta
       integer icirc0,  icirc
* tube -> rgb
       real*8  xp(2), yp(2)        ! location of the projected line segnment
       real*8  xp2, yp2
       real*8  dstnt2
       integer ixp, jyp
       integer ic, jc              ! pixel potition of center of project plane
       real*8  logrc, logrupp, logrbot
* color table
       integer red0, gre0, blu0, red2, gre2, blu2, clr1, clr2, clrg
       real*8  red, gre, blu
       real*8  aaa, bbb, ccc, ddd, fff, cc
       real*8  clrbtm, clrmdl
       parameter(clrbtm = 100, clrmdl = 175) ! color table control index
       logical lwhite
       parameter(lwhite = .false.) ! if true, all lines will be white(gray) ... usually off
* dummy
       integer idummy1, idummy2, idummy3
       real*8  rdummy1, rdummy2, rdummy3, aa
       real*8  xyzdummy(3)
       logical ldummy
       integer ix2,jy2, ixp2(2), jyp2(2)
* constant
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
* functions
       logical lnotdrwn
       real*8  lograd, naiseki3, absvec3
* common
       integer iout, irlog
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       integer ihpix0, jvpix0, ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix0 = 256,          jvpix0 = 256)
       parameter(ihpix2 = ihpix0*ifine, jvpix2 = jvpix0*ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       real*8  rupper, rbottom
       common  /cgcfg/  iout, irlog
       common  /coord/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /images/ dstnt, rgb2
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /varfm1/ rupper, rbottom

       ic = ihpix2 / 2
       jc = jvpix2 / 2

* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,aa)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo

       rdummy1 = rupper
       logrupp = lograd(rdummy1,rbottom,rupper,irlog)
       logrc   = logrupp * 1.1D+00
       rdummy1 = rbottom
       logrbot = lograd(rdummy1,rbottom,rupper,irlog)

* radius of the tube and process number of drawing circumference(surface of tube)
       rtub = logrupp / dfloat(irtube)
       rdummy1 = (rtub * pi * 2.0D+00) / ul
       itheta0 = int(rdummy1 + 0.0001) * 3 + 1
*       if (itheta0 .LT. 8) itheta0 = 8

       do n = 1, nline - 2

* XYZ representation
         do i1 = 1, 4
           i2 = n + 3 - i1
           ra = rtpline(1,i2)
           ra = lograd(ra,rbottom,rupper,irlog)
           th = rtpline(2,i2)
           ph = rtpline(3,i2)
           xyzl(1,i1) = ra * dsin(th) * dcos(ph)
           xyzl(2,i1) = ra * dsin(th) * dsin(ph)
           xyzl(3,i1) = ra * dcos(th)
         enddo

* projection
         do i1 = 1, 2
           do m = 1, 3
             xyzdummy(m) = xyzl(m,i1)
           enddo
           call prjpnt1(xyzdummy,xp2,yp2,dstnt2)
           xp(i1) = xp2
           yp(i1) = yp2
           call pxy2ij(xp2,yp2,ul,ix2,jy2,ic,jc)
           ixp2(i1) = ix2
           jyp2(i1) = jy2
         enddo

         if (((ixp2(1) .GT. 0) .AND. (ixp2(1) .LE. ihpix2) .AND.
     &        (jyp2(1) .GT. 0) .AND. (jyp2(1) .LE. jvpix2)) .OR.
     &       ((ixp2(2) .GT. 0) .AND. (ixp2(2) .LE. ihpix2) .AND.
     &        (jyp2(2) .GT. 0) .AND. (jyp2(2) .LE. jvpix2))) then

* pixel number of projected line segment
         rdummy1 = dabs(xp(1) - xp(2))
         rdummy2 = dabs(yp(1) - yp(2))
         rdummy3 = dmax1(rdummy1, rdummy2)
*         rdummy3 = dmax1(rdummy1, rdummy2, rtub) ! for safety..
         rdummy3 = rdummy3 / ul
         icirc0 = int(rdummy3 + 0.00001) * 2 + 1

* set tangential vector of line segment
         do i2 = 1, 3
           rdummy1 = 0.0D+00
           do i1 = 1, 3
             rdummy1 = rdummy1 + (xyzl(i1,i2) - xyzl(i1,i2+1))**2
           enddo
           rdummy1 = dsqrt(rdummy1)
           do i1 = 1, 3
             nseg(i1,i2) = (xyzl(i1,i2) - xyzl(i1,i2+1)) / rdummy1
           enddo
         enddo

* set normal vector for two circles
*     (the center of two circles are
*      at the begining and ending points of line segment)
         do i1 = 1, 2
           rdummy1 = 0.0D+00
           do i2 = 1, 3
             rdummy1 = rdummy1+(nseg(i2,i1)+nseg(i2,i1+1))**2/4.0D+00
           enddo
           rdummy1 = dsqrt(rdummy1)
           do i2 = 1, 3
             ncc(i2,i1) = (nseg(i2,i1)+nseg(i2,i1+1))/2.0D+00/rdummy1
           enddo
         enddo

* draw circles as assembling tube
         do icirc = 0, icirc0

           rdummy1 = dfloat(icirc) / dfloat(icirc0)
* center position and normal of one circle
           do i1 = 1, 3
             xyzcc(i1) = xyzl(i1,2) *(1.0D+00 - rdummy1)
     &                 + xyzl(i1,3) * rdummy1
             nrmcc(i1) =  ncc(i1,1) *(1.0D+00 - rdummy1)
     &                 +  ncc(i1,2) * rdummy1
           enddo
           call normvec3(nrmcc, rdummy1)

* two normalized vector perpendicular both to each other and to normal vector
           xyzgaiin(1) = 1.0D+00
           xyzgaiin(2) = 0.0D+00
           xyzgaiin(3) = 0.0D+00
           call gaiseki3(nrmcc, xyzgaiin, crimx0)
           call normvec3(crimx0, rdummy1)
           if (rdummy1 .LT. 1.0D-05) then
             xyzgaiin(2) = 1.0D+00
             call gaiseki3(nrmcc, xyzgaiin, crimx0)
             call normvec3(crimx0, rdummy2)
           endif
           call gaiseki3(nrmcc, crimx0, crimy0)
           call normvec3(crimy0, rdummy1)

* draw part of circle
           do itheta = 1, itheta0
             rdummy3 = 2.0D+00 * pi * dfloat(itheta)/dfloat(itheta0)
             do i1 = 1, 3
               ntube(i1) = crimx0(i1) * dcos(rdummy3)
     &                   + crimy0(i1) * dsin(rdummy3)
               crim(i1) = xyzcc(i1) + ntube(i1) * rtub
             enddo
             rdummy1 = absvec3(crim)
             call prjpnt1(crim,xptub,yptub,dstnt2)

             ldummy = lnotdrwn(xyzcp,nprj,crim,xyzeye)
             call pxy2ij(xptub,yptub,ul,ixp,jyp,ic,jc)

             if ((rdummy1 .GE. logrbot) .AND.
     &           (rdummy1 .LE. logrupp) .AND.
     &           (ixp .GT. 0) .AND. (ixp .LE. ihpix2) .AND.
     &           (jyp .GT. 0) .AND. (jyp .LE. jvpix2) .AND.
     &           (.NOT. ldummy)) then

* overwrite or not
               if (dstnt2 .GT. dstnt(ixp,jyp)) then
                 ddd = 0.0D+00
               else
                 ddd = 1.0D+00
                 dstnt(ixp,jyp) = dstnt2
               endif

* light direction
               do m = 1, 3
                 obj2lght(m) = xyzlght(m) - crim(m)
               enddo
               call normvec3(obj2lght,rdummy1)
* line of sight
               do m = 1, 3
                 lofsxyz(m) = crim(m) - xyzeye(m)
               enddo
               call normvec3(lofsxyz,rdummy1)

* color table
               aaa = naiseki3(ntube, obj2lght)
               do i1 = 1, 3
                 mirxyz(i1) = 2.0D+00 * aaa * ntube(i1) - obj2lght(i1)
               enddo
               call normvec3(mirxyz, rdummy1)
               ccc = naiseki3(lofsxyz, mirxyz)
               if (ccc .LT. 0.0D+00) then
                 ccc = 0.0D+00
               else
                 ccc = ccc**2 ! mirror-like effect
               endif
               bbb = naiseki3(ntube, lofsxyz)
               bbb = dabs(bbb)
               aaa = dabs(aaa)
               aaa = aaa**2 ! 0.5D+00
               fff = 0.2D+00
     &             + (aaa + 0.3D+00 * ccc + 0.1D+00) / 1.4D+00
     &             * (bbb + 1.0D+00) / 2.0D+00

*               fff = 0.2D+00 + 0.5D+00 * bbb**3 ! MIND : aaa & ccc is not used ....

               if (fff .GT. 1.0D+00) fff = 1.0D+00
               if (fff .LT. 0.5D+00) then
                 clr1 = 0
                 clr2 = int(fff * 255.0 * 2.0)
               else
                 clr1 = int(fff * 255.0 * 2.0) - 255
                 clr2 = 255
               endif

               if (ivar .EQ. 1) then

                 if (itubecol .EQ. 11) then
                   bbb = naiseki3(ntube, lofsxyz)
                   if (bbb .GE. 0.0D+00) then
                     red0 = clr2
                     gre0 = clr2
                     blu0 = clr1
                   else
** deep green
*                     red0 = clr1
*                     gre0 = clr2 ! green
*                     blu0 = clr1
** light green
                     bbb = dabs(bbb)
                     red = 200.0 * (1.0 + bbb) / 2.0
                     gre = 255.0 * (1.0 + bbb) / 2.0
                     blu = 200.0 * (1.0 + bbb) / 2.0
                     red0 = int(red)
                     gre0 = int(gre)
                     blu0 = int(blu)
* white only
                     if (lwhite) then
                       bbb = dabs(bbb)
                       aa = 255.0D+00 *(bbb + 1.0D+00) * 0.5D+00
                       clrg = int(aa * 0.99999D+00)
                       red0 = clrg
                       gre0 = clrg
                       blu0 = clrg
                     endif

                   endif

                 else if (itubecol .EQ. 17) then
                   aa = 255.0D+00 *(fff + 1.0D+00) * 0.5D+00 ! another choice (2)
                   clrg = int(aa * 0.99999D+00)
                   red0 = clrg
                   gre0 = clrg
                   blu0 = clrg

                 else if (itubecol .EQ. 18) then ! Temp.
                   cc = avec(n) ! avec(nline)
                   call tmpcolor(red0,gre0,blu0,cc)
*                   red = dfloat(red0) * (1.0 + fff) / 2.0
*                   gre = dfloat(gre0) * (1.0 + fff) / 2.0
*                   blu = dfloat(blu0) * (1.0 + fff) / 2.0
                   bbb = naiseki3(ntube, lofsxyz)
                   bbb = dabs(bbb)
                   red = dfloat(red0) * (1.0 + bbb) / 2.0
                   gre = dfloat(gre0) * (1.0 + bbb) / 2.0
                   blu = dfloat(blu0) * (1.0 + bbb) / 2.0
                   red0 = int(red)
                   gre0 = int(gre)
                   blu0 = int(blu)

                 else if (itubecol .EQ. 12) then ! IPS color
                   cc = avec(n) ! avec(nline)
                   rdummy1 = 1.0D+00
                   call ipscolor(red0,gre0,blu0,cc,rdummy1)
                   red = dfloat(red0) * (1.0 + fff) / 2.0
                   gre = dfloat(gre0) * (1.0 + fff) / 2.0
                   blu = dfloat(blu0) * (1.0 + fff) / 2.0
                   red0 = int(red)
                   gre0 = int(gre)
                   blu0 = int(blu)

                 else if (itubecol .EQ. 13) then ! quasi-IPS color ..........................
                   cc = avec(n)
**
*                   rdummy1 = 1.0D+00
*                   call ipscolr2(red0,gre0,blu0,cc,rdummy1) ! quasi-IPS-color
*                   red = dfloat(red0) * (1.0 + fff) / 2.0
*                   gre = dfloat(gre0) * (1.0 + fff) / 2.0
*                   blu = dfloat(blu0) * (1.0 + fff) / 2.0
**
                   rdummy1 = 50.0D+00 ! max in km/s
                   call velcolor(red0,gre0,blu0,cc,rdummy1) ! black-red-white
                   bbb = naiseki3(ntube, lofsxyz)
                   bbb = dabs(bbb)
                   red = dfloat(red0) * (1.0 + bbb) / 2.0
                   gre = dfloat(gre0) * (1.0 + bbb) / 2.0
                   blu = dfloat(blu0) * (1.0 + bbb) / 2.0
                   red0 = int(red)
                   gre0 = int(gre)
                   blu0 = int(blu)

                 else if (itubecol .EQ. 14) then ! open-close in each frame
                   if (lopen1) then
                     red0 = clr2 ! red
                     gre0 = clr1
                     blu0 = clr1
                   else
                     red0 = clr1
                     gre0 = clr2 ! green
                     blu0 = clr1
                   endif
                   bbb = naiseki3(ntube, lofsxyz)
                   if (bbb .GE. 0.0D+00) then ! inside
                     red0 = red0 * 3 / 4
                     gre0 = gre0 * 3 / 4
                     blu0 = blu0 * 3 / 4
                   endif

* white only
                   if (lwhite) then
                     bbb = naiseki3(ntube, lofsxyz)
                     if (bbb .GE. 0.0D+00) then
                       red0 = clr2
                       gre0 = clr2
                       blu0 = clr2
                     else
                       bbb = dabs(bbb)
                       aa = 255.0D+00 *(bbb + 1.0D+00) * 0.5D+00
                       clrg = int(aa * 0.99999D+00)
                       red0 = clrg
                       gre0 = clrg
                       blu0 = clrg
                     endif
                   endif

                 else !                             open-close first frame
                   if (lopen1) then
                     if (lopen2) then ! if open
                       red0 = clr2 ! red
                       gre0 = clr1
                       blu0 = clr1
                     else !             if newly open
                       red0 = clr2 ! yellow
                       gre0 = clr2
                       blu0 = clr1
                     endif
                   else
                     if (lopen2) then ! if newly closed
                       red0 = clr2 ! purple
                       gre0 = clr1
                       blu0 = clr2
                     else !             if closed
** deep green
*                       red0 = clr1
*                       gre0 = clr2 ! green
*                       blu0 = clr1
** light green
                       bbb = dabs(bbb)
                       red = 200.0 * (1.0 + bbb) / 2.0
                       gre = 255.0 * (1.0 + bbb) / 2.0
                       blu = 200.0 * (1.0 + bbb) / 2.0
                       red0 = int(red)
                       gre0 = int(gre)
                       blu0 = int(blu)
                     endif
                   endif
                   bbb = naiseki3(ntube, lofsxyz)
                   if (bbb .GE. 0.0D+00) then ! inside
                     red0 = red0 * 3 / 4
                     gre0 = gre0 * 3 / 4
                     blu0 = blu0 * 3 / 4
                   endif
* white only
                   if (lwhite) then
                     bbb = naiseki3(ntube, lofsxyz)
                     if (bbb .GE. 0.0D+00) then
                       red0 = clr2
                       gre0 = clr2
                       blu0 = clr2
                     else
                       bbb = dabs(bbb)
                       aa = 255.0D+00 *(bbb + 1.0D+00) * 0.5D+00
                       clrg = int(aa * 0.99999D+00)
                       red0 = clrg
                       gre0 = clrg
                       blu0 = clrg
                     endif
                   endif
                 endif

               endif

               if (ivar .EQ. 2) then
                 if (itubecol .EQ. 2) then
                   cc = avec(n) ! avec(nline)
                   rdummy1 = 1.0D+00
                   call ipscolor(red0,gre0,blu0,cc,rdummy1)
                   red = dfloat(red0) * (1.0 + fff) / 2.0
                   gre = dfloat(gre0) * (1.0 + fff) / 2.0
                   blu = dfloat(blu0) * (1.0 + fff) / 2.0
                   red0 = int(red)
                   gre0 = int(gre)
                   blu0 = int(blu)
                 else
                   bbb = naiseki3(ntube, lofsxyz)
                   if (bbb .GE. 0.0D+00) then
                     red0 = clr2
                     gre0 = clr1
                     blu0 = clr2
                   else
                     red0 = clr1
                     gre0 = clr2
                     blu0 = clr2
                   endif
                 endif
               endif

               red2 = ichar(rgb2(1,ixp,jyp))
               gre2 = ichar(rgb2(2,ixp,jyp))
               blu2 = ichar(rgb2(3,ixp,jyp))
* put RGB data set
               rdummy3 = (1.0-ddd)*dfloat(red2)+ddd*dfloat(red0)
               idummy1 = int(rdummy3)
               rdummy3 = (1.0-ddd)*dfloat(gre2)+ddd*dfloat(gre0)
               idummy2 = int(rdummy3)
               rdummy3 = (1.0-ddd)*dfloat(blu2)+ddd*dfloat(blu0)
               idummy3 = int(rdummy3)
               if (idummy1 .GT. 255) idummy1 = 255
               if (idummy2 .GT. 255) idummy2 = 255
               if (idummy3 .GT. 255) idummy3 = 255
               if (idummy1 .LT.   0) idummy1 =   0
               if (idummy2 .LT.   0) idummy2 =   0
               if (idummy3 .LT.   0) idummy3 =   0
               rgb2(1,ixp,jyp) = char(idummy1)
               rgb2(2,ixp,jyp) = char(idummy2)
               rgb2(3,ixp,jyp) = char(idummy3)
             endif

           enddo ! End of do-loop for one circle
         enddo !   End of do-loop for one segment of tube

       endif

       enddo !     End of do-loop for one tube

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
       parameter(linemax = 10000)
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
       parameter(pi = 3.14159265358979D+00)
       real*8  sclatrtp
       logical lstop
       real*8  aa
* common
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       common  /sclscl/ sclr
       common  /coord/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp
*
       dpi = pi / dfloat(jj)

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
*     +    / dfloat(kk) * 0.5D+00 * dfloat(idrct)
       dr = (rr(1) - rr(0)) / 4.0
       call rungevec(rtp1, rtp2, dr, avecl, ivar)
       avecl = sclatrtp(rtp2)
       avec(0) = avecl

       if ((      lfromsun) .AND. (rtp2(1) .LT. rtp1(1))) idrct = -1
       if ((.NOT. lfromsun) .AND. (rtp2(1) .GT. rtp1(1))) idrct = -1

*       write(*,*) idrct

* estimate magnetic vector and get coordinate of line
       nline = 0
       ra = rtpline(1,0)
       th = rtpline(2,0)
       ph = rtpline(3,0)
       if (lfromsun) then
         rtp1(1) = ra + 1.0D-03
       else
         rtp1(1) = ra - 1.0D-03
       endif
       rtp1(2) = th
       rtp1(3) = ph
       lstop = .false.
       do while((.NOT. lstop) .AND.
     &          (th .GE. dpi) .AND.
     &          (th .LE. (pi - dpi)) .AND.
     &          (nline .LT. linemax))

         dr = ra * pi / dfloat(jj) / 1.0D+01 * dfloat(idrct)
         nline = nline + 1

         call rungevec(rtp1, rtp2, dr, avecl, ivar) ! 1 --> 2
*         call rungevec_xyz(rtp1, rtp2, dr, avecl, ivar)
         if (rtp2(1) .GT. rupper) then
           lstop = .true.
           aa = (rupper - rtp1(1)) / (rtp2(1) - rtp1(1))
           rtp2(1) = rtp1(1) + (rtp2(1) - rtp1(1)) * aa
           rtp2(2) = rtp1(2) + (rtp2(2) - rtp1(2)) * aa
           rtp2(3) = rtp1(3) + (rtp2(3) - rtp1(3)) * aa
         endif
         if (rtp2(1) .LT. rbottom) then
           lstop = .true.
           aa = (rbottom - rtp1(1)) / (rtp2(1) - rtp1(1))
           rtp2(1) = rtp1(1) + (rtp2(1) - rtp1(1)) * aa
           rtp2(2) = rtp1(2) + (rtp2(2) - rtp1(2)) * aa
           rtp2(3) = rtp1(3) + (rtp2(3) - rtp1(3)) * aa
         endif

         avecl = sclatrtp(rtp2)
         ra = rtp2(1)
         th = rtp2(2)
         ph = rtp2(3)
         if (ph .LT. 0.0D+00) then
           ph = ph + 2.0D+00 * pi
         else if (ph .GT. 2.0D+00 * pi) then
           ph = dmod(ph, 2.0D+00 * pi)
         endif
         rtp2(3) = ph

         do i = 1, 3
           rtpline(i,nline) = rtp2(i)
           rtp1(i)          = rtp2(i)
         enddo

         avec(nline) = avecl

       enddo        ! End of do-while-loop

 300   continue

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
       real*8  bbb
       real*8  k1(3), k2(3), k3(3), k4(3), k0(3)
       real*8  b1(3), b2(3), b3(3), b4(3), b0(3)
       real*8  rtp3(3), w(3)
       real*8  pi, ph
       integer k
       parameter(pi = 3.14159265358979D+00)
* common etc.
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1), theta(0:jj), phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common  /coord/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp

       w(1) = 1.0D+00
       w(2) = 1.0D+00 / rtp1(1)
       w(3) = 1.0D+00 / rtp1(1) / dsin(rtp1(2))

       do k = 1, 3
         rtp3(k) = rtp1(k)
       enddo
       ph = rtp3(3)
       ph = ph + 2.0D+00 * pi
       ph = dmod(ph, 2.0D+00 * pi)
       rtp3(3) = ph
       call vecatrtp(b1,rtp3)
       bbb = b1(1) ** 2 + b1(2) ** 2 + b1(3) ** 2
       bbb = dsqrt(bbb)

       if (ivar .EQ. 1) then
         avecl = bbb
       else if (ivar .EQ. 2) then
         avecl = b1(1)
       else
         avecl = 0.0D+00
       endif

       do k = 1, 3
         b1(k) = b1(k) / bbb
         k1(k) = b1(k) * dr * w(k)
         rtp3(k) = rtp1(k) + k1(k) * 0.5D+00
       enddo
       ph = rtp3(3)
       ph = ph + 2.0D+00 * pi
       ph = dmod(ph, 2.0D+00 * pi)
       rtp3(3) = ph
       call vecatrtp(b2,rtp3)
       bbb = b2(1) ** 2 + b2(2) ** 2 + b2(3) ** 2
       bbb = dsqrt(bbb)
       do k = 1, 3
         b2(k) = b2(k) / bbb
         k2(k) = b2(k) * dr * w(k)
         rtp3(k) = rtp1(k) + k2(k) * 0.5D+00
       enddo
       ph = rtp3(3)
       ph = ph + 2.0D+00 * pi
       ph = dmod(ph, 2.0D+00 * pi)
       rtp3(3) = ph
       call vecatrtp(b3,rtp3)
       bbb = b3(1) ** 2 + b3(2) ** 2 + b3(3) ** 2
       bbb = dsqrt(bbb)
       do k = 1, 3
         b3(k) = b3(k) / bbb
         k3(k) = b3(k) * dr * w(k)
         rtp3(k) = rtp1(k) + k3(k)
       enddo
       ph = rtp3(3)
       ph = ph + 2.0D+00 * pi
       ph = dmod(ph, 2.0D+00 * pi)
       rtp3(3) = ph
       call vecatrtp(b4,rtp3)
       bbb = b4(1) ** 2 + b4(2) ** 2 + b4(3) ** 2
       bbb = dsqrt(bbb)
       do k = 1, 3
         b4(k) = b4(k) / bbb
         k4(k) = b4(k) * dr * w(k)
         b0(k) =(b1(k) + 2.0D+00*b2(k) + 2.0D+00*b3(k) + b4(k))
     &        / 6.0D+00
       enddo
       bbb = b0(1) ** 2 + b0(2) ** 2 + b0(3) ** 2
       bbb = dsqrt(bbb)
       do k = 1, 3
         b0(k) = b0(k) / bbb
         k0(k) = b0(k) * dr * w(k)
         rtp2(k) = rtp1(k) + k0(k)
       enddo

       return
       end

*
* ---------------------------------------------------------------
*
       subroutine rungevec_xyz(rtp1, rtp2, dr, avecl, ivar)  ! rtp1 --> rtp2
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
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common  /coord/  rr, theta, phi
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
       th = dacos(xyza(3) / ra * 0.99999D+00)
       aa = dsqrt(xyza(1)**2 + xyza(2)**2)
       if (aa .LT. 1.0D-03) aa = 1.0D-03
       ph = dacos(xyza(1) / aa * 0.99999D+00)
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
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  vcr(0:ii+1,0:jj,0:kk+1)
       real*8  vct(0:ii+1,0:jj,0:kk+1)
       real*8  vcp(0:ii+1,0:jj,0:kk+1)
       common  /coord/  rr, theta, phi
       common  /vecvec/ vcr, vct, vcp

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

*       write(*,*) rtp(1),rtp(2),rtp(3)
*       write(*,*) ijk(1),ijk(2),ijk(3)
*       write(*,*) vec(1),vec(2),vec(3)

       return
       end


*
* --------------------------------------------------------------
* scalar at rtp
*
       real*8 function sclatrtp(rtp)
       implicit none
* interface
       real*8  rtp(3)
* local
       real*8  dra1, dra2, dth1, dth2, dph1, dph2
       integer i1, i2, j1, j2, k1, k2, ijk(3)
       real*8  ans, ra, th, ph ! , aaa
*
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
*
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       common  /coord/  rr, theta, phi
       common  /sclscl/ sclr

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

       ans =(sclr(i1,j1,k1)*dth2*dph2+sclr(i1,j2,k1)*dth1*dph2
     &     + sclr(i1,j1,k2)*dth2*dph1+sclr(i1,j2,k2)*dth1*dph1)*dra2
     &     +(sclr(i2,j1,k1)*dth2*dph2+sclr(i2,j2,k1)*dth1*dph2
     &     + sclr(i2,j1,k2)*dth2*dph1+sclr(i2,j2,k2)*dth1*dph1)*dra1

       sclatrtp = ans

       return
       end



* ===================================================================
* read files
*    READDATA
*    READINIT
* ===================================================================
*
       subroutine readdata(rt,nnn) ! XX3 be filled
       implicit none
*
       integer rt, nnn
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
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
       common  /var31/  ro3, pg3
       common  /var32/  ur3, ut3, up3
       common  /var34/  br3, bt3, bp3

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
           lmhd = .false.
           do 111 k = 0, kk + 1
           do 111 j = 0, jj
           do 111 i = 0, ii + 1
             ro3(i,j,k) = 1.0D+00 ! dummy !
             pg3(i,j,k) = 1.0D+00
             ur3(i,j,k) = 0.0D+00
             ut3(i,j,k) = 0.0D+00
             up3(i,j,k) = 0.0D+00
             br3(i,j,k) = 0.0D+00
             bt3(i,j,k) = 0.0D+00
             bp3(i,j,k) = 0.0D+00
 111       continue
         endif
       endif

       return
       end

*
** -------------------------
*
       subroutine setsclr(ivar,rotref,nref)
       implicit none
       integer ivar,rotref,nref
* local
       real*8  maxsclr, minsclr
       real*8  aa, bb, cc
       real*8  sumx, sumxx, stddev, ave
       integer i, j, k, m
       real*8  r1
       real*8  gamma
       real*8  pi
       parameter(pi = 3.14159265358979D+00)
       real*8  lofsxyz(3), magrtp(3), magxyz(3), cros(3)
       real*8  nprj(3), xyzeye(3)
       real*8  ra, th, ph
       real*8  rdummy1
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  prkrho(0:ii+1)
* functions
       real*8  parkvr, absvec3 ! , naiseki3
* common
       integer iout, irlog
       real*8  sclr(0:ii+1,0:jj,0:kk+1)
       real*8  ro0(0:ii+1,0:jj,0:kk+1)
       real*8  pg0(0:ii+1,0:jj,0:kk+1)
       real*8  ur0(0:ii+1,0:jj,0:kk+1)
       real*8  ut0(0:ii+1,0:jj,0:kk+1)
       real*8  up0(0:ii+1,0:jj,0:kk+1)
       real*8  br0(0:ii+1,0:jj,0:kk+1)
       real*8  bt0(0:ii+1,0:jj,0:kk+1)
       real*8  bp0(0:ii+1,0:jj,0:kk+1)
       real*8  ro1(0:ii+1,0:jj,0:kk+1)
       real*8  pg1(0:ii+1,0:jj,0:kk+1)
       real*8  ur1(0:ii+1,0:jj,0:kk+1)
       real*8  ut1(0:ii+1,0:jj,0:kk+1)
       real*8  up1(0:ii+1,0:jj,0:kk+1)
       real*8  br1(0:ii+1,0:jj,0:kk+1)
       real*8  bt1(0:ii+1,0:jj,0:kk+1)
       real*8  bp1(0:ii+1,0:jj,0:kk+1)
       real*8  ro3(0:ii+1,0:jj,0:kk+1)
       real*8  pg3(0:ii+1,0:jj,0:kk+1)
       real*8  ur3(0:ii+1,0:jj,0:kk+1)
       real*8  ut3(0:ii+1,0:jj,0:kk+1)
       real*8  up3(0:ii+1,0:jj,0:kk+1)
       real*8  br3(0:ii+1,0:jj,0:kk+1)
       real*8  bt3(0:ii+1,0:jj,0:kk+1)
       real*8  bp3(0:ii+1,0:jj,0:kk+1)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       common  /var01/  ro0, pg0
       common  /var02/  ur0, ut0, up0
       common  /var04/  br0, bt0, bp0
       common  /var11/  ro1, pg1
       common  /var12/  ur1, ut1, up1
       common  /var14/  br1, bt1, bp1
       common  /var31/  ro3, pg3
       common  /var32/  ur3, ut3, up3
       common  /var34/  br3, bt3, bp3
       common  /coord/  rr, theta, phi
       common  /gasgas/ gamma
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
       common  /sclscl/ sclr
       common  /cgcfg/  iout, irlog

* location of eye
       call gaiseki3(xyzxp,xyzyp,nprj)
       call normvec3(nprj,aa)
       do m = 1, 3
         xyzeye(m) = xyzcp(m) + nprj(m) * dproj
       enddo

       do i = 0, ii + 1
         r1 = rr(i)
         prkrho(i) = 1.0D+00 / parkvr(r1,gamma) / r1**2
       enddo

** set var
       if (ivar .EQ. 18) then
         write(11,'(A)') ' std(T) or Max/Min at i'
         do 918 i = 0, ii + 1
           sumx  = 0.0D+00
           sumxx = 0.0D+00
           maxsclr = -1.0D+10
           minsclr =  1.0D+10
           do k = 1, kk
           do j = 0, jj
             bb = pg1(i,j,k) / ro1(i,j,k)
             if (bb .GT. maxsclr) maxsclr = bb
             if (bb .LT. minsclr) minsclr = bb
             sumx  = sumx  + bb
             sumxx = sumxx + bb**2
           enddo
           enddo
           sumx  = sumx  / dfloat(kk) / dfloat(jj+1)
           sumxx = sumxx / dfloat(kk) / dfloat(jj+1)
           ave  = sumx
           stddev = sumxx - sumx**2 + 1.0D-10
           if (stddev .LT. 1.0D-20) stddev = 1.0D-20
           stddev = dsqrt(stddev)
           do j = 0, jj
           do k = 1, kk
             bb = pg1(i,j,k) / ro1(i,j,k)

*             bb = (bb - ave) / stddev
             bb = (bb - minsclr) / (maxsclr - minsclr + 1.0D-10)

             sclr(i,j,k) = bb
           enddo
           enddo
 918     continue


       else if (ivar .EQ. 14) then
         write(11,'(A)') ' rho'

         do 914 i = 0, ii + 1
           do j = 0, jj
           do k = 1, kk
             sclr(i,j,k) = ro1(i,j,k)
           enddo
           enddo
 914     continue

       else if (ivar .EQ. 17) then
         write(11,'(A)') ' std(rho)'

         do 917 i = 0, ii + 1
           sumx  = 0.0D+00
           sumxx = 0.0D+00
           do k = 1, kk
           do j = 0, jj
             bb = ro1(i,j,k)
             sumx  = sumx  + bb
             sumxx = sumxx + bb**2
           enddo
           enddo
           sumx  = sumx  / dfloat(kk) / dfloat(jj+1)
           sumxx = sumxx / dfloat(kk) / dfloat(jj+1)
           ave  = sumx
           stddev = sumxx - sumx**2 + 1.0D-10
           if (stddev .LT. 1.0D-10) stddev = 1.0D-10
           stddev = dsqrt(stddev)
           do j = 0, jj
           do k = 1, kk
             bb = ro1(i,j,k)
             bb = (bb - ave) / stddev
*             bb = bb - 0.1 ! Geta........ if needed.
*             if (bb .LT. 0.0D+00) then
*               bb = 0.0D+00
*             else
*               bb = bb
*             endif
             sclr(i,j,k) = bb
           enddo
           enddo
 917     continue

       else if (ivar .EQ. 13) then
         write(11,'(A)') ' std(V)'

         do 913 i = 0, ii + 1
           sumx  = 0.0D+00
           sumxx = 0.0D+00
           do k = 1, kk
           do j = 0, jj
             bb = ur1(i,j,k)**2 + ut1(i,j,k)**2 + up1(i,j,k)**2
             bb = 1.0D+00 / (bb * ro1(i,j,k) + 1.0D-10)
             sumx  = sumx  + bb
             sumxx = sumxx + bb**2
           enddo
           enddo
           sumx  = sumx  / dfloat(kk) / dfloat(jj+1)
           sumxx = sumxx / dfloat(kk) / dfloat(jj+1)
           ave  = sumx
           stddev = sumxx - sumx**2 + 1.0D-10
           if (stddev .LT. 1.0D-20) stddev = 1.0D-20
           stddev = dsqrt(stddev)
           do j = 0, jj
           do k = 1, kk
             bb = ur1(i,j,k)**2 + ut1(i,j,k)**2 + up1(i,j,k)**2
             bb = 1.0D+00 / (bb * ro1(i,j,k) + 1.0D-10)
             bb = (bb - ave) / stddev
             if (bb .LT. 0.0D+00) then
               bb = 0.0D+00
             else
               bb = bb
             endif
             sclr(i,j,k) = bb
           enddo
           enddo
 913     continue

       else if (ivar .EQ. 12) then
         write(11,'(A)') ' std(rho) +alpha'

         do 912 i = 0, ii + 1
           sumx  = 0.0D+00
           sumxx = 0.0D+00
           do k = 1, kk
           do j = 0, jj
             sumx  = sumx  + ro1(i,j,k)
             sumxx = sumxx + ro1(i,j,k)**2
           enddo
           enddo
           sumx  = sumx  / dfloat(kk) / dfloat(jj+1)
           sumxx = sumxx / dfloat(kk) / dfloat(jj+1)
           ave  = sumx
           stddev = sumxx - sumx**2 + 1.0D-10
           if (stddev .LT. 1.0D-20) stddev = 1.0D-20
           stddev = dsqrt(stddev)
           do j = 0, jj
           do k = 1, kk
             bb = (ro1(i,j,k) - ave) / stddev
**
*             bb = bb + 1.0D+00 ! <== filter, if needed.
**
             if (bb .LT. 0.0D+00) then
               bb = 0.0D+00
             else
               bb = bb ! **1.5D+00
             endif
**
             sclr(i,j,k) = bb
           enddo
           enddo
 912     continue

       else if (ivar .EQ. 11) then
         write(11,'(A)') ' V^2 etc'

         do 911 i = 0, ii + 1
           sumx  = 0.0D+00
           sumxx = 0.0D+00
           do k = 1, kk
           do j = 0, jj
             bb = ur1(i,j,k)**2 + ut1(i,j,k)**2 + up1(i,j,k)**2
             bb = dsqrt(bb + 1.0D-10) * 10.0D+00
             if (bb .GT. 1.0D+00) then
               bb = 0.0D+00
             else
               bb = 1.0D+00 / (bb**2 + 1.0D+00)
             endif
             sumx  = sumx  + bb
             sumxx = sumxx + bb**2
           enddo
           enddo
           sumx  = sumx  / dfloat(kk) / dfloat(jj+1)
           sumxx = sumxx / dfloat(kk) / dfloat(jj+1)
           ave  = sumx
           stddev = sumxx - sumx**2 + 1.0D-10
           if (stddev .LT. 1.0D-20) stddev = 1.0D-20
           stddev = dsqrt(stddev)
           do j = 0, jj
           do k = 1, kk
             bb = ur1(i,j,k)**2 + ut1(i,j,k)**2 + up1(i,j,k)**2
             bb = dsqrt(bb + 1.0D-10) * 10.0D+00
             if (bb .GT. 1.0D+00) then
               bb = 0.0D+00
             else
               bb = 1.0D+00 / (bb**2 + 1.0D+00)
             endif
             bb = (bb - ave) / stddev
             if (bb .LT. 0.0D+00) then
               bb = 0.0D+00
             else
               bb = bb ! **1.5D+00
             endif
             sclr(i,j,k) = bb
           enddo
           enddo
 911     continue

       else if (ivar .EQ. 10) then
         write(11,'(A)') ' rho r^2 '

         do 901 k = 1, kk
         do 901 j = 0, jj
         do 901 i = 0, ii + 1
           bb = ro1(i,j,k) * rr(i)**2
           sclr(i,j,k) = bb**2
 901     continue

       else if (ivar .EQ. 9) then ! rho (N/Nref)
         write(11,'(A)') ' B Los, and rho/rho_0 etc '

         call readdata(rotref,nref) ! XX3 be filled
         do k = 0, kk + 1
         do j = 0, jj
         do i = 0, ii + 1
           ro0(i,j,k) = ro3(i,j,k)
           pg0(i,j,k) = pg3(i,j,k)
           ur0(i,j,k) = ur3(i,j,k)
           ut0(i,j,k) = ut3(i,j,k)
           up0(i,j,k) = up3(i,j,k)
           br0(i,j,k) = br3(i,j,k)
           bt0(i,j,k) = bt3(i,j,k)
           bp0(i,j,k) = bp3(i,j,k)
         enddo
         enddo
         enddo
         do 900 k = 1, kk
         do 900 j = 0, jj
         do 900 i = 0, ii + 1

           ra = rr(i)
           th = theta(j)
           ph = phi(k)
           magrtp(1) = br1(i,j,k)
           magrtp(2) = bt1(i,j,k)
           magrtp(3) = bp1(i,j,k)
           magxyz(1) = magrtp(1) * dsin(th) * dcos(ph) ! UNDER
     &               + magrtp(2) * dcos(th) * dcos(ph)
     &               - magrtp(3)            * dsin(ph)
           magxyz(2) = magrtp(1) * dsin(th) * dsin(ph)
     &               + magrtp(2) * dcos(th) * dsin(ph)
     &               + magrtp(3)            * dcos(ph)
           magxyz(3) = magrtp(1) * dcos(th)
     &               - magrtp(2) * dsin(th)

           lofsxyz(1) = ra * dsin(th) * dcos(ph) - xyzeye(1)
           lofsxyz(2) = ra * dsin(th) * dsin(ph) - xyzeye(2)
           lofsxyz(3) = ra * dcos(th)            - xyzeye(3)
           call normvec3(lofsxyz,rdummy1)
           call gaiseki3(magxyz, lofsxyz, cros)
           aa = absvec3(cros) ! * rr(i)**2
           bb = ro1(i,j,k) / ro0(i,j,k) - 1.0D+00
           if (bb .GT. 1.00D-20) then
             sclr(i,j,k) =  aa * (bb + 1.0)
           else
             sclr(i,j,k) = 0.0D+00
           endif

* reset !
           rdummy1 = ro1(i,j,k) / ro0(i,j,k)
           if (rdummy1 .GT. 1.07D+00) then
             sclr(i,j,k) = rdummy1 - 1.0D+00
           else
             sclr(i,j,k) = 0.0D+00
           endif

 900     continue

       else if (ivar .EQ. 8) then ! expected IPS
         write(11,'(A)') ' IPS-vel '

         do 800 k = 1, kk
         do 800 j = 0, jj
         do 800 i = 0, ii + 1
           ra = rr(i)
           th = theta(j)
           ph = phi(k)
           magrtp(1) = ur1(i,j,k) ! Mind "mag" does not yet stand for Magnetic field.
           magrtp(2) = ut1(i,j,k)
           magrtp(3) = up1(i,j,k)
           magxyz(1) = magrtp(1) * dsin(th) * dcos(ph) ! UNDER
     &               + magrtp(2) * dcos(th) * dcos(ph)
     &               - magrtp(3)            * dsin(ph)
           magxyz(2) = magrtp(1) * dsin(th) * dsin(ph)
     &               + magrtp(2) * dcos(th) * dsin(ph)
     &               + magrtp(3)            * dcos(ph)
           magxyz(3) = magrtp(1) * dcos(th)
     &               - magrtp(2) * dsin(th)
           lofsxyz(1) = ra * dsin(th) * dcos(ph) - xyzeye(1)
           lofsxyz(2) = ra * dsin(th) * dsin(ph) - xyzeye(2)
           lofsxyz(3) = ra * dcos(th)            - xyzeye(3)
           call gaiseki3(magxyz, lofsxyz, cros)
           rdummy1 = absvec3(cros)
           call normvec3(lofsxyz,rdummy1)
*           rdummy1 = naiseki3(magxyz, lofsxyz)
*           rdummy1 = dabs(rdummy1)
           sclr(i,j,k) = rdummy1 * ro1(i,j,k) !  * rr(i)**2

 800     continue

       else if (ivar .EQ. 16) then
         write(11,'(A)') ' Pg/park etc '

         do 716 k = 1, kk
         do 716 j = 0, jj
         do 716 i = 0, ii + 1
           rdummy1 = pg1(i,j,k) / prkrho(i)**gamma - 0.85
           if (rdummy1 .GT. 0.0) then
             sclr(i,j,k) = rdummy1 ! **3
           else
             sclr(i,j,k) = 0.0D+00
           endif
 716     continue

       else if (ivar .EQ. 7) then
         write(11,'(A)') ' rho/park etc '

         do 700 k = 1, kk
         do 700 j = 0, jj
         do 700 i = 0, ii + 1
*           rdummy1 = ro1(i,j,k) / prkrho(i) - 1.00
           rdummy1 = ro1(i,j,k) / prkrho(i) - 0.85
           if (rdummy1 .GT. 0.0) then
             sclr(i,j,k) = rdummy1 ! **3
           else
             sclr(i,j,k) = 0.0D+00
           endif
 700     continue

       else if (ivar .EQ. 6) then
         write(11,'(A)') ' dV '

         do 600 k = 1, kk
         do 600 j = 1, jj - 1
           sclr(0,j,k) = 0.0D+00
           do 600 i = 0, ii + 1
             aa = (theta(j+1) - theta(j-1)) / 2.0D+00 * rr(i)
             bb = (phi(k+1) - phi(k-1))
     &          / 2.0D+00 * rr(i) * dsin(theta(j))
             cc = (ur1(i,j+1,k) - ur1(i,j,k))**2 / aa**2
     &          + (ur1(i,j,k) - ur1(i,j-1,k))**2 / aa**2
     &          + (ur1(i,j,k+1) - ur1(i,j,k))**2 / bb**2
     &          + (ur1(i,j,k) - ur1(i,j,k+1))**2 / bb**2
             sclr(i,j,k) = cc * rr(i)**2 - 20.0D+00
 600     continue
         do 610 j = 1, jj - 1
         do 610 i = 0, ii + 1
           sclr(i,j,   0) = sclr(i,j,kk)
           sclr(i,j,kk+1) = sclr(i,j, 1)
 610     continue
         do i = 0, ii + 1
           aa = 0.0D+00
           bb = 0.0D+00
           do k = 1, kk
             aa = aa + sclr(i,   1,k)
             bb = bb + sclr(i,jj-1,k)
           enddo
           aa = aa / dfloat(kk)
           bb = bb / dfloat(kk)
           do k = 1, kk
             sclr(i, 0,k) = aa
             sclr(i,jj,k) = bb
           enddo
         enddo

       else if (ivar .EQ. 5) then
         write(11,'(A)') ' B Los etc..'

         do 500 k = 1, kk
         do 500 j = 0, jj
         do 500 i = 0, ii + 1
           ra = rr(i)
           th = theta(j)
           ph = phi(k)
           magrtp(1) = br1(i,j,k)
           magrtp(2) = bt1(i,j,k)
           magrtp(3) = bp1(i,j,k)
           magxyz(1) = magrtp(1) * dsin(th) * dcos(ph) ! UNDER
     &               + magrtp(2) * dcos(th) * dcos(ph)
     &               - magrtp(3)            * dsin(ph)
           magxyz(2) = magrtp(1) * dsin(th) * dsin(ph)
     &               + magrtp(2) * dcos(th) * dsin(ph)
     &               + magrtp(3)            * dcos(ph)
           magxyz(3) = magrtp(1) * dcos(th)
     &               - magrtp(2) * dsin(th)

           lofsxyz(1) = ra * dsin(th) * dcos(ph) - xyzeye(1)
           lofsxyz(2) = ra * dsin(th) * dsin(ph) - xyzeye(2)
           lofsxyz(3) = ra * dcos(th)            - xyzeye(3)
           call normvec3(lofsxyz,rdummy1)
           call gaiseki3(magxyz, lofsxyz, cros)
           aa = absvec3(cros) ! * rr(i)**2
           bb = ro1(i,j,k) / prkrho(i)
           sclr(i,j,k) = (aa * bb)**2
 500     continue

       else if (ivar .EQ. 4) then
         write(11,'(A)') ' (drho/prk)^20 '

         do 400 k = 1, kk
         do 400 j = 0, jj
         do 400 i = 0, ii + 1
           sclr(i,j,k) =  ro1(i,j,k)
     &                  *(ro1(i,j,k) / prkrho(i))**20
 400     continue

       else if (ivar .EQ. 3) then
         write(11,'(A)') ' drho/prk '

         do 300 k = 1, kk
         do 300 j = 0, jj
         do 300 i = 0, ii + 1
           ra = rr(i)
           th = theta(j)
           ph = phi(k)
           magrtp(1) = br1(i,j,k)
           magrtp(2) = bt1(i,j,k)
           magrtp(3) = bp1(i,j,k)
           magxyz(1) = magrtp(1) * dsin(th) * dcos(ph) ! UNDER
     &               + magrtp(2) * dcos(th) * dcos(ph)
     &               - magrtp(3)            * dsin(ph)
           magxyz(2) = magrtp(1) * dsin(th) * dsin(ph)
     &               + magrtp(2) * dcos(th) * dsin(ph)
     &               + magrtp(3)            * dcos(ph)
           magxyz(3) = magrtp(1) * dcos(th)
     &               - magrtp(2) * dsin(th)

           lofsxyz(1) = ra * dsin(th) * dcos(ph) - xyzeye(1)
           lofsxyz(2) = ra * dsin(th) * dsin(ph) - xyzeye(2)
           lofsxyz(3) = ra * dcos(th)            - xyzeye(3)
           call normvec3(lofsxyz,rdummy1)
           call gaiseki3(magxyz, lofsxyz, cros)
           aa = absvec3(cros)
           aa = dabs(aa)
           aa = (aa**0.5 + 0.5)

           bb = ro1(i,j,k) / prkrho(i)
           if (bb .GT. 1.0) then
             if (ro1(i,j,k) .GE. 10.0)
     &         bb = (bb + 0.4)**3 * dlog10(ro1(i,j,k))
           else
             bb = 0.0
           endif

           sclr(i,j,k) = aa * bb / 10.0
 300     continue

       else if (ivar .EQ. 2) then
         write(11,'(A)') ' drho/rho - 0.15 '

         do 200 j = 1, jj - 1
         do 200 k = 1, kk
           sclr(0,j,k) = 0.0D+00
           do 200 i = 0, ii + 1
             sclr(i,j,k) =(dabs(ro1(i,j+1,k  )-ro1(i,j  ,k  ))
     &                    +dabs(ro1(i,j  ,k  )-ro1(i,j-1,k  ))
     &                    +dabs(ro1(i,j  ,k+1)-ro1(i,j  ,k  ))
     &                    +dabs(ro1(i,j  ,k  )-ro1(i,j  ,k-1)))
     &                   / 4.0D+00 / ro1(i,j,k)
     &                   - 0.15D+00
 200     continue
         do 210 i = 0, ii + 1
         do 210 j = 1, jj - 1
           sclr(i,j,   0) = sclr(i,j,kk)
           sclr(i,j,kk+1) = sclr(i,j, 1)
 210     continue
         do i = 0, ii + 1
           aa = 0.0D+00
           bb = 0.0D+00
           do k = 1, kk
             aa = aa + sclr(i,   1,k)
             bb = bb + sclr(i,jj-1,k)
           enddo
           aa = aa / dfloat(kk)
           bb = bb / dfloat(kk)
           do k = 1, kk
             sclr(i, 0,k) = aa
             sclr(i,jj,k) = bb
           enddo
         enddo

       else if (ivar .EQ. 1) then

         write(11,'(A)') ' Br '
         do 100 k = 1, kk
         do 100 j = 0, jj
         do 100 i = 0, ii + 1
           sclr(i,j,k) = br1(i,j,k)
 100     continue

       else
         write(*,*) ' Ivar is not correct ! '
       endif

* overlap
       do 10 j = 0, jj
       do 10 i = 0, ii + 1
         sclr(i,j,   0) = sclr(i,j,kk)
         sclr(i,j,kk+1) = sclr(i,j, 1)
 10    continue

       maxsclr = -1.0D+10
       minsclr =  1.0D+10
       do 20 k = 1, kk
       do 20 j = 0, jj
       do 20 i = 0, ii + 1
         if (sclr(i,j,k) .GT. maxsclr) maxsclr = sclr(i,j,k)
         if (sclr(i,j,k) .LT. minsclr) minsclr = sclr(i,j,k)
 20    continue

       write(11,'('' Scalar Index = '',i2,$)') ivar
       write(11,'('' : Max / Min (i>0) = '',e11.5,'' / '',e11.5)')
     &   maxsclr, minsclr
       if (iout .GT. 0) then
         write(*,'('' Scalar Index = '',i2,$)') ivar
         write(*,'('' : Max / Min (i>0) = '',e11.5,'' / '',e11.5)')
     &     maxsclr, minsclr
       endif

       return
       end


*
** -----------------------------------------------------------
*
       real*8 function parkvr(r1,gamma)
       implicit none
* interface
       real*8  r1, gamma
* local
       real*8  vp, dd, vg
* functions
       real*8  parkff, dparkff

       vp = 1.0D-03
       dd = 1.0D+02
       if (r1 .GT. 1.0D+00) vp = 1.0D+01
       do while(dd .GE. 1.0D-04)
         vg = vp - parkff(r1,vp,gamma) / dparkff(r1,vp,gamma)
         dd = dabs(1.0D+00 - vg / vp)
         vp = vg
       enddo

       parkvr = vp

       return
       end


*
** -------------
*
       real*8 function parkff(rr,vr,gamma)
       implicit none
       real*8 rr, vr, gamma, cc
       real*8 ans, aa
       cc = 0.0D+00
       aa = gamma - 1.0D+00
       aa = dabs(aa)
       if (aa .LT. 1.0D-05) then
         ans = 0.5D+00 *(vr**2 - 1.0D+00)
     &       - dlog(vr) - 2.0D+00 * dlog(rr)
     &       - 2.0D+00 *(1.0D+00 / rr - 1.0D+00) + cc
       else
         ans = 0.5D+00 *(vr**2 - 1.0D+00)
     &       +((rr**2 * vr)**(-aa) - 1.0D+00) / aa
     &       - 2.0D+00 *(1.0D+00 / rr - 1.0D+00) + cc
       endif
       parkff = ans
       return
       end


*
** ------------
*
       real*8 function dparkff(rr,vr,gamma)
       implicit none
       real*8 rr, vr, gamma
       real*8 ans, aa
       aa = gamma - 1.0D+00
       aa = dabs(aa)
       if (aa .LT. 1.0D-05) ans = vr - 1.0D+00 / vr
       if (aa .GE. 1.0D-05)
     &   ans = vr - rr**(-2.0D+00*aa) * vr**(-gamma)
       dparkff = ans
       return
       end


*
* ----------------------------------------------------
*
       subroutine readinit(gamma,t0,length0,v0,omega)
       implicit none
* interface
       real*8  gamma, t0, length0, v0, omega
* local
       real*8  a1, a2 ! , aa
       integer idummy
*       integer iorien, igrid
       integer i1, i2, i, j, k
       integer imax, jmax, kmax, cmax
       character*50 strdummy
       character*50 strdumm2
       character*50 cmnt(30)
*
       real*8  rsun
       parameter(rsun=6.96D+10)     ! Solar radius (cm)
*
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       common  /coord/  rr, theta, phi

       open(unit=1,file='init.dat',status='old')
* read coordinate
       imax = ii
       jmax = jj
       kmax = kk
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
* read comment
       cmax = 0
 300   continue
         read(1,'(A)',END=399) strdummy
         cmax = cmax + 1
         cmnt(cmax) = strdummy
*
         i1 = index(strdummy, 'gamma')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') gamma
         endif
*
         i1 = index(strdummy, 'omega')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') omega
         endif
*
         i1 = index(strdummy, 't0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') t0
         endif
*
         i1 = index(strdummy, 'r0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') length0
         endif
*
         i1 = index(strdummy, 'v0')
         if (i1 .NE. 0) then
           i2 = index(strdummy,' = ')
           strdumm2 = strdummy(i2+3:50)
           read(strdumm2,'(e11.5)') v0
         endif
*
         goto 300
 399   continue
       close(1)

       length0 = length0 / rsun
       omega = omega * t0 / (3.6D+03 * 2.4D+01 * 1.80D+02) * 3.14D+00

       write(11,'(A)') ' Reading INIT.DAT'
       write(11,'('' Gamma   = '',f9.4)') gamma
       write(11,'('' Omega   = '',f9.4)') omega
       write(11,'('' t0      = '',e9.3,'' sec'')') t0
       write(11,'('' Length0 = '',f9.4,'' R_sun'')') length0
       write(11,'('' RR : '',f9.4,'' to '',f9.4)') rr(0),rr(ii+1)

       return
       end


*
*** -----------------------------------------------------------
*
* 0 : 48   a :  97    A : 65
* 1 : 49   z : 122    Z : 90
* 9 : 57
*
* Itrape = 0      : non-transparent
*          1      :     transparent
*          others : not drawn
*
*** -----------------------------------------------------------
*
       subroutine charwri(charout,charlen,ix,jy,itrape,rgb0)
       implicit none
       character*20 charout
       integer charlen, ix, jy, itrape
       integer ihpix0, jvpix0
       parameter(ihpix0 = 256, jvpix0 = 256)
       integer ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix2 = ihpix0*ifine, jvpix2 = jvpix0*ifine)
       character*1 rgb0(3,ihpix2, jvpix2)     ! RGB data array
       integer ich, jch
       parameter(ich = 8, jch = 11)
       character*1 cdummy
       integer i, j, m, i0, j0, n, l
       integer idummy1, idummy2, idummy3
       integer num(8,11,0:9), cap(8,11,1:26)
       integer num0(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer num1(8,11)    /0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,
     &        0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0/
       integer num2(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,1,0,
     &        0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,
     &        0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1/
       integer num3(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,
     &        0,0,0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer num4(8,11)    /0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,
     &        0,0,0,1,1,1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,
     &        0,1,1,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,1,1,1,1,1,1,
     &        0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0/
       integer num5(8,11)    /1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,1,1,0,
     &        0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer num6(8,11)    /0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,1,
     &        0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,
     &        1,1,1,0,0,1,1,0,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer num7(8,11)    /1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,
     &        0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,
     &        0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,
     &        0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0/
       integer num8(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,
     &        0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer num9(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        0,1,1,0,0,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,
     &        1,1,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
*
       integer capa(8,11)    /0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,
     &        0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,
     &        0,1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,1,1,1,1,1,
     &        0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1/
       integer capb(8,11)    /1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,0,
     &        1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,0,1,1,1,1,1,1,0,0/
       integer capc(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer capd(8,11)    /1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,0,1,1,1,1,1,1,0,0/
       integer cape(8,11)    /1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
     &        1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1/
       integer capf(8,11)    /1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
     &        1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0/
       integer capg(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,
     &        1,1,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,1,0,0,1,1,1,1,1,1/
       integer caph(8,11)    /1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1/
       integer capi(8,11)    /0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0/
       integer capj(8,11)    /0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,
     &        0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,
     &        0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer capk(8,11)    /1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,0,
     &        1,1,0,0,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,
     &        1,1,1,1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,1,1,0,0,
     &        1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,1/
       integer capl(8,11)    /1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1/
       integer capm(8,11)    /1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,
     &        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     &        1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1/
       integer capn(8,11)    /1,1,0,0,0,1,1,0,1,1,1,0,0,1,1,0,
     &        1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,0,
     &        1,1,1,1,0,1,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,0,
     &        1,1,0,0,1,1,1,0,1,1,0,0,1,1,1,0,1,1,0,0,0,1,1,0/
       integer capo(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer capp(8,11)    /1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,1,1,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,
     &        1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0/
       integer capq(8,11)    /0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,1,
     &        1,1,0,0,1,1,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1/
       integer capr(8,11)    /1,1,1,1,1,1,1,0,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,
     &        1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,1/
       integer caps(8,11)    /0,1,1,1,1,1,1,0,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,
     &        0,0,0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,1,1,1,1,1,0/
       integer capt(8,11)    /1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0/
       integer capu(8,11)    /1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0/
       integer capv(8,11)    /1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,
     &        0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,
     &        0,0,1,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0/
       integer capw(8,11)    /1,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,
     &        1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     &        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,
     &        0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0/
       integer capx(8,11)    /1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,
     &        0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0,
     &        0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,0,
     &        0,1,1,0,0,1,1,0,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1/
       integer capy(8,11)    /1,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,
     &        0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,
     &        0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0/
       integer capz(8,11)    /1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,
     &        0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,
     &        0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,
     &        0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1/

       do j = 1, 11
       do i = 1, 8
         num(i,j,0) = num0(i,j)
         num(i,j,1) = num1(i,j)
         num(i,j,2) = num2(i,j)
         num(i,j,3) = num3(i,j)
         num(i,j,4) = num4(i,j)
         num(i,j,5) = num5(i,j)
         num(i,j,6) = num6(i,j)
         num(i,j,7) = num7(i,j)
         num(i,j,8) = num8(i,j)
         num(i,j,9) = num9(i,j)
         cap(i,j, 1) = capa(i,j)
         cap(i,j, 2) = capb(i,j)
         cap(i,j, 3) = capc(i,j)
         cap(i,j, 4) = capd(i,j)
         cap(i,j, 5) = cape(i,j)
         cap(i,j, 6) = capf(i,j)
         cap(i,j, 7) = capg(i,j)
         cap(i,j, 8) = caph(i,j)
         cap(i,j, 9) = capi(i,j)
         cap(i,j,10) = capj(i,j)
         cap(i,j,11) = capk(i,j)
         cap(i,j,12) = capl(i,j)
         cap(i,j,13) = capm(i,j)
         cap(i,j,14) = capn(i,j)
         cap(i,j,15) = capo(i,j)
         cap(i,j,16) = capp(i,j)
         cap(i,j,17) = capq(i,j)
         cap(i,j,18) = capr(i,j)
         cap(i,j,19) = caps(i,j)
         cap(i,j,20) = capt(i,j)
         cap(i,j,21) = capu(i,j)
         cap(i,j,22) = capv(i,j)
         cap(i,j,23) = capw(i,j)
         cap(i,j,24) = capx(i,j)
         cap(i,j,25) = capy(i,j)
         cap(i,j,26) = capz(i,j)
       enddo
       enddo

       do 110 m = 1, charlen
         do 110 j = 1, jch
         do 110 i = 1, ich
           l = 0
           cdummy = charout(m:m)
           idummy1 = ichar(cdummy)
           if ((idummy1 .GE. 48) .AND. (idummy1 .LE. 57)) then
             n = idummy1 - 48
             l = num(i,j,n)
           endif
           if ((idummy1 .GE. 65) .AND. (idummy1 .LE. 90)) then
             n = idummy1 - 64
             l = cap(i,j,n)
           endif
           if ((idummy1 .GE. 97) .AND. (idummy1 .LE. 122)) then
             n = idummy1 - 96 ! small -> cap.
             l = cap(i,j,n)
           endif

           if (l .EQ. 1) then
             i0 = ix + ich * m + i
             j0 = jy + jch - j
             if (itrape .EQ. 0) then !        not-transparent letter
               idummy1 = 0 ! 255
               idummy2 = 0 ! 255
               idummy3 = 0 ! 255
               rgb0(1,i0,j0) = char(idummy1)
               rgb0(2,i0,j0) = char(idummy2)
               rgb0(3,i0,j0) = char(idummy3)
             else if (itrape .EQ. 1) then !     transparent letter
               idummy1 = (ichar(rgb0(1,i0,j0)) + 255) / 2
               idummy2 = (ichar(rgb0(2,i0,j0)) + 255) / 2
               idummy3 = (ichar(rgb0(3,i0,j0)) + 255) / 2
               if (idummy1 .GT. 255) idummy1 = 255
               if (idummy2 .GT. 255) idummy2 = 255
               if (idummy3 .GT. 255) idummy3 = 255
               if (idummy1 .LT.   0) idummy1 =   0
               if (idummy2 .LT.   0) idummy2 =   0
               if (idummy3 .LT.   0) idummy3 =   0
               rgb0(1,i0,j0) = char(idummy1)
               rgb0(2,i0,j0) = char(idummy2)
               rgb0(3,i0,j0) = char(idummy3)
             endif
           endif
 110   continue
       return
      end



*
*** ----------------------------------------------------------------------
*
* degree and radian
*
*** ----------------------------------------------------------------------
*
       subroutine deg2rad2(deg1,deg2,rad1,rad2)
       implicit none
* interface
       real*8   deg1,deg2,rad1,rad2
* local
       real*8   pi
       parameter(pi = 3.14159265358979D+00)

       rad1 = deg1 / 180.0D+00 * pi
       rad2 = deg2 / 180.0D+00 * pi

       return
       end


* ------------
*
       subroutine rad2deg2(rad1,rad2,deg1,deg2)
       implicit none
* interface
       real*8   rad1,rad2,deg1,deg2
* local
       real*8   pi
       parameter(pi = 3.14159265358979D+00)

       deg1 = rad1 * 180.0D+00 / pi
       deg2 = rad2 * 180.0D+00 / pi

       return
       end


*
* ------------------------------------------------------------------
*
* vector operators
*
* ------------------------------------------------------------------
*
*
* -----------------
*
       real*8 function naiseki3(vec1, vec2)
       implicit none
* interface
       real*8 vec1(3), vec2(3)
* local
       real*8 ans

       ans = vec1(1) * vec2(1)
     &     + vec1(2) * vec2(2)
     &     + vec1(3) * vec2(3)
       naiseki3 = ans

       return
       end


*
* -----------------
*
       subroutine gaiseki3(xyz1, xyz2, xyzout)
       implicit none
* interface
       real*8  xyz1(3), xyz2(3), xyzout(3)

       xyzout(1) = xyz1(2)*xyz2(3) - xyz1(3)*xyz2(2)
       xyzout(2) = xyz1(3)*xyz2(1) - xyz1(1)*xyz2(3)
       xyzout(3) = xyz1(1)*xyz2(2) - xyz1(2)*xyz2(1)

       return
       end


*
* ----------------
*
       subroutine copyvec3(vec1,vec2)
       implicit none
* interface
       real*8 vec1(3), vec2(3)
* local
       integer m

       do m = 1, 3
         vec2(m) = vec1(m)
       enddo

       return
       end

*
* ----------------
*
       real*8 function absvec3(vec)
       implicit none
* interface
       real*8 vec(3)
* local
       real*8 rdummy

       rdummy = vec(1)**2 + vec(2)**2 + vec(3)**2
       rdummy = dsqrt(rdummy)
       absvec3 = rdummy

       return
       end

*
*** -----------------------------------------------------
*
       subroutine initrgb(rgb2,dstnt,bkrgb)
       implicit none
* interface
       integer ihpix0, jvpix0, ihpix2, jvpix2, ifine
       parameter(ifine = 2)
       parameter(ihpix0 = 256,          jvpix0 = 256)
       parameter(ihpix2 = ihpix0*ifine, jvpix2 = jvpix0*ifine)
       character*1 rgb2(3,ihpix2, jvpix2)
       real*8  dstnt(ihpix2, jvpix2)
       integer bkrgb(3)
* local
       logical lgrdbak
       parameter(lgrdbak = .false.)
       integer i0, j0
       integer bkrgb2(3)
       real    aa, bb

       do 100 j0 = 1, jvpix2
       do 100 i0 = 1, ihpix2
* gradation...
         if (lgrdbak) then
* partern 1 vertical
           bb = dfloat(j0) / dfloat(jvpix2) - 0.5D+00
* partern 2 diagonal
*           bb =(dfloat(j0) / dfloat(jvpix2)
*     &        + dfloat(i0) / dfloat(ihpix2)) * 0.5 - 0.5
* partern 3 horizontal
*           bb = dfloat(i0) / dfloat(ihpix2) - 0.5D+00

           aa = dfloat(bkrgb(1)) * (1.0D+00 + 0.2D+00 * bb)
           if (aa .LT.   0.00D+00) aa = 0.0D+00
           if (aa .GT. 255.01D+00) aa = 255.01D+00
           bkrgb2(1) = int(aa)
           aa = dfloat(bkrgb(2)) * (1.0D+00 + 0.2D+00 * bb)
           if (aa .LT.   0.00D+00) aa = 0.0D+00
           if (aa .GT. 255.01D+00) aa = 255.01D+00
           bkrgb2(2) = int(aa)
           aa = dfloat(bkrgb(3)) * (1.0D+00 + 0.2D+00 * bb)
           if (aa .LT.   0.00D+00) aa = 0.0D+00
           if (aa .GT. 255.01D+00) aa = 255.01D+00
           bkrgb2(3) = int(aa)
         else
* uniform
           bkrgb2(1) = bkrgb(1)
           bkrgb2(2) = bkrgb(2)
           bkrgb2(3) = bkrgb(3)
         endif

         rgb2(1,i0,j0) = char(bkrgb2(1))
         rgb2(2,i0,j0) = char(bkrgb2(2))
         rgb2(3,i0,j0) = char(bkrgb2(3))
         dstnt(i0, j0)= 1.00D+20
 100   continue

       return
       end

*
* -----------------------------------------
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
       parameter(ii = 72, jj = 32, kk = 64)
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
     &         (rtp(2).LE.theta(m+1))) ijk(2) = m
         enddo
       endif

       if (rtp(3) .LE. 0.0D+00)      rtp(3) = rtp(3) + 2.0D+00 * pi
       if (rtp(3) .GT. 2.0D+00 * pi) rtp(3) = rtp(3) - 2.0D+00 * pi
       ijk(3) = 0
       do m = 0, kk
         if ((rtp(3).GE.phi(m)).AND.
     &       (rtp(3).LE.phi(m+1))) ijk(3) = m
       enddo

       return
       end


*
** -----------------------------
*
       subroutine pxy2ij(xp,yp,ul,ix,jy,ic,jc)
       implicit none
* interface
       real*8  xp, yp, ul
       integer ix, jy, ic, jc
* local
       real*8  xp2, yp2, aa

       xp2 = xp / ul
       yp2 = yp / ul
       ix = ic + int(xp2)
       aa = dmod(xp2, 1.0D+00)
       if (aa .GE. 0.5D+00) ix = ix + 1
       jy = jc + int(yp2)
       aa = dmod(yp2, 1.0D+00)
       if (aa .GE. 0.5D+00) jy = jy + 1

       return
       end



* -------------------------------
*
       subroutine prjpnt1(xyz,xpt,ypt,dstntt)
       implicit none
* interface
       real*8  xyz(3), xpt, ypt, dstntt
* local variables
       integer nn
       parameter(nn = 3)
       real*8  vecma(nn), kmat(nn,nn), vecmb(nn)
       integer i
       real*8  nrml(3), xyzeye(3)
       real*8  dd
* common block
       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul

* normal of projection plane
       call gaiseki3(xyzxp, xyzyp, nrml)
* position of eye
       do i = 1, 3
         xyzeye(i) = xyzcp(i) + nrml(i) * dproj
       enddo

       dd = 0.0D+00
       do i = 1, 3
         kmat(i,1) =   xyzxp(i)
         kmat(i,2) =   xyzyp(i)
         kmat(i,3) = -(xyz(i) - xyzeye(i))
         vecma(i)  =   nrml(i) * dproj
         dd = dd + (xyz(i) - xyzeye(i))**2
       enddo
       dd = dsqrt(dd)
       call linear3(vecma, kmat, vecmb)

* location on the projection plane
       xpt = vecmb(1)
       ypt = vecmb(2)
* distance between object and projection plane
       dstntt = dd * (1.0D+00 - vecmb(3))

       return
       end




* -----------------------------
*
*       subroutine prjpnt3(xyz,xg,yg,dstnt2) ! now not used
*       implicit none
*       integer numpoint
*       parameter(numpoint = 3)
*       real*8  xyz(3,numpoint)
*       real*8  xg(numpoint), yg(numpoint), dstnt2(numpoint)
** local variables
*       integer nn
*       parameter(nn = 3)
*       real*8  vecma(nn), kmat(nn,nn), vecmb(nn)
*       integer i, l
*       real*8  nrml(3), xyzeye(3)
*       real*8  dd
** common block
*       real*8  xyzcp(3), dproj, xyzxp(3), xyzyp(3), xyzlght(3), ul
*       common  /varfm0/ xyzcp, dproj, xyzxp, xyzyp, xyzlght, ul
*
*       nrml(1) = xyzxp(2)*xyzyp(3)-xyzxp(3)*xyzyp(2)
*       nrml(2) = xyzxp(3)*xyzyp(1)-xyzxp(1)*xyzyp(3)
*       nrml(3) = xyzxp(1)*xyzyp(2)-xyzxp(2)*xyzyp(1)
*       do i = 1, 3
*         xyzeye(i) = xyzcp(i) + nrml(i) * dproj
*       enddo
*
*       do l = 1, numpoint
*         dd = 0.0D+00
*         do i = 1, 3
*           kmat(i,1) =   xyzxp(i)
*           kmat(i,2) =   xyzyp(i)
*           kmat(i,3) = -(xyz(i,l) - xyzeye(i))
*           vecma(i)  =   nrml(i) * dproj
*           dd = dd + (xyz(i,l) - xyzeye(i))**2
*         enddo
*         dd = dsqrt(dd)
*         call linear3(vecma, kmat, vecmb)
*         xg(l) = vecmb(1)
*         yg(l) = vecmb(2)
*         dstnt2(l) = dd * (1.0D+00 - vecmb(3))
*       enddo
*       return
*       end
*


*
* solve linear equation  using Jacob-method --------------------------
*
       subroutine linear3(vecma, kmat, vecmb) ! V_a = Mat_k * V_b
       implicit none
       integer nn
       parameter(nn = 3)
       real*8  vecma(nn), kmat(nn,nn), vecmb(nn)
       real*8  kmat2(nn,nn), kmati(nn,nn)
       integer i, j
       do i = 1, nn
         do j = 1, nn
           kmat2(i,j) = kmat(i,j)
         enddo
       enddo
       call invmat3(kmat2, kmati)
       call mulmtvc3(kmati, vecma, vecmb)
       return
       end



*
* get inversion matrix ----------------------------------------------
*
       subroutine invmat3(matin, matout)
       implicit none
       integer nn
       parameter(nn = 3)
       real*8   matin(nn,nn), matout(nn,nn)
       integer  pivo(nn)
       integer  i, j, k, lh, l, ii
       real*8   h, max
       real*8   e(nn,nn), a(nn,nn)
       do i = 1, nn
         pivo(i) = 0
       enddo
       do i = 1, nn
         do j = 1, nn
           if (i.EQ.j) e(i,j) = 1.0D+00
           if (i.NE.j) e(i,j) = 0.0D+00
           a(i,j) = matin(i,j)
         enddo
       enddo
       do k = 1, nn
         max = 0.0D+00
         do i = 1, nn
           h = abs(a(i,k))
           if (h .GT. max) then
             ii = 1
             if (pivo(i) .EQ. 0) then
               max = h
               j = i
             endif
           endif
         enddo
         pivo(j) = k
         h = 1.0D+00 / a(j,k)
         do lh = 1, nn
           a(j,lh) = a(j,lh) * h
           e(j,lh) = e(j,lh) * h
         enddo
         do i = 1, nn
           if (i .NE. j) then
             h = a(i,k)
             do lh = 1, nn
               a(i,lh) = a(i,lh) - a(j,lh) * h
               e(i,lh) = e(i,lh) - e(j,lh) * h
             enddo
           endif
         enddo
       enddo
       do k = 1, nn
         do i = 1, nn
*           if (pivo(i) .EQ. k) l = i
*           do lh = 1, nn
*             a(k, lh) = e(l,lh)
*           enddo
           if (pivo(i) .EQ. k) then
             l = i
             do lh = 1, nn
               a(k, lh) = e(l,lh)
             enddo
           endif
         enddo
       enddo
       do i = 1, nn
         do j = 1, nn
           matout(i,j) = a(i,j)
         enddo
       enddo
       return
       end



*
* multiple matrix and vector  ------------------------------
*
       subroutine mulmtvc3(mat, vec1, vecout)
       implicit none
       integer nn
       parameter(nn = 3)
       real*8  mat(nn,nn), vec1(nn), vecout(nn)
       real*8  sum
       integer i, j
       do i = 1, nn
         sum = 0.0D+00
         do j = 1, nn
           sum = sum + mat(i,j) * vec1(j)
         enddo
         vecout(i) = sum
       enddo
       return
       end


*
** which side -----------------------------------------------
*
       logical function sameside(xyz0,norm,xyz1,xyz2)
       implicit none
* interface
       real*8  xyz0(3), norm(3), xyz1(3), xyz2(3)
* local
       integer m
       real*8  aa, bb, cc

       aa = 0.0D+00
       bb = 0.0D+00
       do m = 1, 3
         aa = aa + norm(m) * (xyz1(m) - xyz0(m))
         bb = bb + norm(m) * (xyz2(m) - xyz0(m))
       enddo

       cc = aa * bb
       if (cc .LT. 0.0D+00) then
         sameside = .false. ! different side
       else
         sameside = .true.  ! same side
       endif

       return
       end
*
** --------
*
       logical function lnotdrwn(xyz0,norm,xyz1,xyz2)
       implicit none
* interface
       real*8  xyz0(3), norm(3), xyz1(3), xyz2(3)
* local
       integer m
       real*8  aa, bb, cc
       logical ldummy1
       logical ldummy2

       aa = 0.0D+00
       bb = 0.0D+00
       do m = 1, 3
         aa = aa + norm(m) * (xyz1(m) - xyz0(m))
         bb = bb + norm(m) * (xyz2(m) - xyz0(m))
       enddo

       cc = aa * bb
       if (cc .LT. 0.0D+00) then
         ldummy1 = .false. ! different side of the projection plain
       else
         ldummy1 = .true.  ! same side
       endif

       ldummy2 = .false. ! second condition, here always turned-off
* cut the region(s) : not used now....

*      ldummy2 =(xyz1(3) .GT. 0.0D+00) ! draw only south

** one choice ...
*       aa = 0.0D+00
*       bb = 0.0D+00
*       cc = 0.0D+00
*       do m = 1, 3
*         aa = aa + xyz1(m) * xyz2(m)
*         bb = bb + xyz1(m)**2
*         cc = cc + xyz2(m)**2
*       enddo
*       aa = aa / dsqrt(bb * cc) ! normalized inner product of LoS and Position

*       ldummy2 = (aa .GT. 0.75D+00) ! NOT draw near LoS
*       ldummy2 = (aa .GT. 0.00D+00) ! NOT draw near-side : draw only far-side
*       ldummy2 = (.NOT. ((aa .GT. -0.15) .AND. (aa .LT. 0.15)))  ! draw only nearby plane of ski

* 2nd choice
*       ldummy2 = ((xyz1(1) .LT. 0.0D+00) .AND.
*     +            (xyz1(2) .GT. 0.0D+00) .AND.
*     +            (xyz1(3) .GT. 0.0D+00))

* 3rd choice
*       aa = datan2(xyz1(1),xyz1(2))
*       bb = datan2(xyz2(1),xyz2(2))
*       cc = bb - aa
*       aa = dcos(cc)
*       bb = dsin(cc)
*       ldummy2 =((xyz1(3) .GT. 0.0D+00) .AND.
*     &           ((aa .GT.  0.00D+00) .AND.
*     &            (bb .GT. -0.50D+00) .AND.
*     &            (bb .LT.  0.90D+00)))

       lnotdrwn = (ldummy1 .OR. ldummy2)

       return
       end


******** end of this file ********************

