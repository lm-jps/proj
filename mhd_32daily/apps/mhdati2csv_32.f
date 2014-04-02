       program mhdati
       implicit none
* number of grid for each three direction.
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  ro3(0:ii+1,0:jj,0:kk+1),pg3(0:ii+1,0:jj,0:kk+1)
       real*8  ur3(0:ii+1,0:jj,0:kk+1),ut3(0:ii+1,0:jj,0:kk+1)
       real*8  up3(0:ii+1,0:jj,0:kk+1),br3(0:ii+1,0:jj,0:kk+1)
       real*8  bt3(0:ii+1,0:jj,0:kk+1),bp3(0:ii+1,0:jj,0:kk+1)
       integer ipickup, ncr, nstep, j, k
       real*8  brbot(0:jj,0:kk+1),vrbot(0:jj,0:kk+1)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       real*8  theta2(0:jj),phi2(0:kk+1)
       logical lmhd
       real*8  rpickup
       real*8  gamma,omega,v0,b0,n0,tmp0,r0
*       character*8 strdate1 ! fujitsu at nagoya
       character*9 strdate1
       character*7 strdate2
       character*3  mmchr
       integer imm, iyy, idd

** substitute calling date() at Stan.
*       open(12,file='date.txt',status='old')
*       read(12,*) strdate1
*       close(12)
*       read(strdate1,'(i2,x,a3,x,i2)') idd, mmchr, iyy
*       call chr2mon(mmchr,imm)
*       write(strdate2,'(2i2.2,''.'',i2.2)') mod(iyy,100), imm, idd
**
       open(12,file='carr.txt',status='old')
       read(12,*) strdate1
       close(12)
      
** at Nagoya
*       call date(strdate1)
*       write(strdate2,'(a2,a2,''.'',a2)')
*     &     strdate1(1:2),strdate1(4:5),strdate1(7:8)

       call readinit(gamma,omega,v0,b0,n0,tmp0,r0,rr,theta,phi)

       write(*,*) 'Input I for pickup map, Ncr and Nstep'
       read(*,*) ipickup, ncr, nstep

       rpickup = rr(ipickup) * r0 / 7.0D+10 ! in Rs
*       rpickup = rr(ipickup) / rr(0) ! in Rs

       write(*,*) 'Pick up at r (Rs) = ', rpickup

       call readmhd(ncr,nstep,ro3,pg3,ur3,ut3,up3,br3,bt3,bp3,lmhd)
       do j = 0, jj
         do k = 0, kk + 1
           brbot(j,k) = br3(ipickup,j,k) * b0 * 1.0D+05 ! in nT
           vrbot(j,k) = ur3(ipickup,j,k) * v0 / 1.0D+05 ! in km/s
         enddo
       enddo
       do j = 0, jj
         theta2(j) = 90.0D+00
     &             - theta(j) * 180.0D+00 / 3.14159265358797D+00
       enddo
       do k = 0, kk + 1
         phi2(k)   = phi(k)   * 180.0D+00 / 3.14159265358797D+00
       enddo
       call csvout(brbot,vrbot,phi2,theta2,rpickup,strdate2)

       stop
       end


*
* --------------------------------------------------------------
*
      subroutine chr2mon(mmchr,mm1)
      implicit none
      character*3  mmchr
      integer mm1
      mm1 = -1
      if (mmchr .EQ. 'Jan') mm1 =  1
      if (mmchr .EQ. 'Feb') mm1 =  2
      if (mmchr .EQ. 'Mar') mm1 =  3
      if (mmchr .EQ. 'Apr') mm1 =  4
      if (mmchr .EQ. 'May') mm1 =  5
      if (mmchr .EQ. 'Jun') mm1 =  6
      if (mmchr .EQ. 'Jul') mm1 =  7
      if (mmchr .EQ. 'Aug') mm1 =  8
      if (mmchr .EQ. 'Sep') mm1 =  9
      if (mmchr .EQ. 'Oct') mm1 = 10
      if (mmchr .EQ. 'Nov') mm1 = 11
      if (mmchr .EQ. 'Dec') mm1 = 12
      if (mm1 .LT. 0) write(*,*) ' ERR at chr2mon : ',mmchr
      return
      end

*
* --------------------------------------------------------------------
*
       subroutine csvout(brbot,vrbot,phi,theta,rpickup,strdate)
       implicit none
* number of grid for each three direction.
       integer jj, kk
       parameter(jj = 32, kk = 64)
       real*8  brbot(0:jj,0:kk+1),vrbot(0:jj,0:kk+1)
       real*8  theta(0:jj),phi(0:kk+1) ! this is in degree
       real*8  rpickup
       character*7 strdate
*
       integer j, k
       character*30 fmtstr
       real*8  aa

* supress warining
       aa = rpickup

       write(fmtstr,'(''f8.3'',i3,''('','',f8.3)'')') kk+1

       open(unit=11,file='solwind_vr.csv',status='unknown')
*       write(11,'(f8.3,'','',$)') rpickup
*       write(11,'(x,a7,'','',$)') strdate
       write(11,'(a7,''0,'',$)') strdate
       do k = 0, kk-1
         write(11,'(f8.3,'','',$)') phi(k)
       enddo
       write(11,'(f8.3)') phi(kk)
       do j = 0, jj
         write(11,'(f8.3,'','',$)') theta(j)
         do k = 0, kk-1
           write(11,'(f8.3,'','',$)') vrbot(j,k)
         enddo
         write(11,'(f8.3)') vrbot(j,kk)
       enddo
       close(11)

       open(unit=11,file='solwind_br.csv',status='unknown')
*       write(11,'(f8.3,'','',$)') rpickup
*       write(11,'(x,a7,'','',$)') strdate
       write(11,'(a7,''0,'',$)') strdate
       do k = 0, kk-1
         write(11,'(f8.3,'','',$)') phi(k)
       enddo
       write(11,'(f8.3)') phi(kk)
       do j = 0, jj
         write(11,'(f8.3,'','',$)') theta(j)
         do k = 0, kk-1
           write(11,'(f8.3,'','',$)') brbot(j,k)
         enddo
         write(11,'(f8.3)') brbot(j,kk)
       enddo
       close(11)

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
       real*8  a1, a2, aa, t0
       integer idummy
       integer iorien, igrid
       character*50 strdummy
       character*50 strdumm2
       integer i1, i2
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
       real*8  rr(0:ii+1),theta(0:jj),phi(0:kk+1)
       character*50  cmnt(30)
       integer imax, jmax, kmax, cmax

       open(unit=1,file='init.dat',status='old')
* read coordinate
       imax = 0
       jmax = 0
 200   continue
         read(1,*,END=299) idummy, a1, a2
         if (idummy .EQ. 0) goto 299
         aa = dfloat(idummy) / 1.0D+03 + 1.0D-05
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
       write(*,'('' rr   (  0) : '',f10.5)') rr(0)
       write(*,'('' rr   ('',i3,'') : '',f10.5)') imax, rr(imax)
       write(*,'('' theta(  0) : '',f10.5)') theta(0)
       write(*,'('' theta('',i3,'') : '',f10.5)') jmax, theta(jmax)
       write(*,'('' phi  (  0) : '',f10.5)') phi(0)
       write(*,'('' phi  ('',i3,'') : '',f10.5)') kmax, phi(kmax)
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
       subroutine readmhd(rt,nnn,ro3,pg3,ur3,ut3,up3,br3,bt3,bp3,lmhd)
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
           write(*,*) ' File not found ',flname
           lmhd = .false.
           do k = 0, kk + 1
             do j = 0, jj
               do i = 0, ii + 1
                 ro3(i,j,k) = 1.0D+00
                 pg3(i,j,k) = 1.0D+00
                 ur3(i,j,k) = 1.0D+00
                 ut3(i,j,k) = 1.0D+00
                 up3(i,j,k) = 1.0D+00
                 br3(i,j,k) = 1.0D+00
                 bt3(i,j,k) = 1.0D+00
                 bp3(i,j,k) = 1.0D+00
               enddo
             enddo
           enddo
         endif
       endif

       do 201 k = 0, kk + 1
       do 201 j = 0, jj
       do 201 i = 0, ii + 1
         if (.NOT. ((ro3(i,j,k) .LT. 1.0D-30) .OR.
     &              (ro3(i,j,k) .GE.-1.0D-30))) ro3(i,j,k) = 1.0D-30
         if (.NOT. ((pg3(i,j,k) .LT. 1.0D-30) .OR.
     &              (pg3(i,j,k) .GE.-1.0D-30))) pg3(i,j,k) = 1.0D-30
         if (.NOT. ((ur3(i,j,k) .LT. 1.0D-30) .OR.
     &              (ur3(i,j,k) .GE.-1.0D-30))) ur3(i,j,k) = 1.0D-30
         if (.NOT. ((ut3(i,j,k) .LT. 1.0D-30) .OR.
     &              (ut3(i,j,k) .GE.-1.0D-30))) ut3(i,j,k) = 1.0D-30
         if (.NOT. ((up3(i,j,k) .LT. 1.0D-30) .OR.
     &              (up3(i,j,k) .GE.-1.0D-30))) up3(i,j,k) = 1.0D-30
         if (.NOT. ((br3(i,j,k) .LT. 1.0D-30) .OR.
     &              (br3(i,j,k) .GE.-1.0D-30))) br3(i,j,k) = 1.0D-30
         if (.NOT. ((bt3(i,j,k) .LT. 1.0D-30) .OR.
     &              (bt3(i,j,k) .GE.-1.0D-30))) bt3(i,j,k) = 1.0D-30
         if (.NOT. ((bp3(i,j,k) .LT. 1.0D-30) .OR.
     &              (bp3(i,j,k) .GE.-1.0D-30))) bp3(i,j,k) = 1.0D-30
 201   continue

       return
       end

