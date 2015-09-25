       program dat2cr
       implicit none
       integer yyyy, mm, cr, iday
       integer yy1, mm1
       real    rlon, rdd
       integer imenu
       real    aa, b0
       character*9 sdate ! for Intel
       character*3 smonth
       integer itm(6)
       integer mdiday, idshift, idoy, ihh
       real    rday, rhh
       integer idatetime(8)
       character*10 bigben(3)
* function
       integer mon2char, dyofmth, mmdd2doy
*
       character*120 flname1
       integer rtstart, rtend
       parameter(rtstart = -9, rtend = 2999)
       integer yystrt(rtstart:rtend)
       integer mmstrt(rtstart:rtend)
       real    ddstrt(rtstart:rtend)
       common /crrot1/ yystrt, mmstrt
       common /crrot2/ ddstrt

       read(*,'(a)') flname1
       call readdate(flname1)

       call date_and_time(bigben(1),bigben(2),bigben(3),idatetime)
       itm(3) = idatetime(1) ! YYYY
       itm(2) = idatetime(2) ! MM
       itm(1) = idatetime(3) ! Day
       itm(4) = idatetime(4) ! diff to UT
       itm(5) = idatetime(5) ! local HH
       itm(6) = idatetime(6) ! local MM

** YYYY.MM.DD_HH:00:00 form
*       write(*,*) '--------------------------------------------'
       yy1 =  itm(3)
       mm1 =  itm(2)
       iday = itm(1)
*       idshift = -6
       idshift = 0

       rhh = float(itm(5)) + float(-itm(4)+itm(6)) / 60.0 + 0.00001
       ihh = int(rhh)
       if (ihh .LT. 0) then
         ihh = ihh + 24
         idshift = idshift - 1
       endif
       if (ihh .GT. 23) then
         ihh = ihh - 24
         idshift = idshift + 1
       endif
       call shiftday(yy1,mm1,iday,idshift)
       write(*,'(i4.4,i2.2,i2.2,''_120000'')')
     &     yy1, mm1, iday
**
*       write(*,'(''['',i4.4,''.'',i2.2,''.'',i2.2,''_'',
*     &              i2.2,'':00:00_TAI]'')')
*     &     yy1, mm1, iday, ihh
**
*       write(*,'(i4.4,''.'',i2.2,''.'',i2.2,''_'',i2.2)')
*     &     yy1, mm1, iday, ihh
**


*       call date(sdate)
*       read(sdate,'(i2.2,1x,a3,1x,i2.2)') itm(1),smonth,itm(3)
*       call chr2mon(smonth, mm)
*       itm(2) = mm
*       itm(3) = itm(3) + 2000 ! YY->CCYY

* 2000 Jan 1st was 2556th MDI day

*       mdiday = 2556 - 1
*       if (itm(3) .GT. 2000) then
*         do yyyy = 2000, itm(3) - 1
*           mdiday = mdiday + 365
*           if (mod(yyyy,4) .EQ. 0) mdiday = mdiday + 1
*         enddo
*       endif
*       
*       yyyy = itm(3)
*       mm   = itm(2)
*       iday = itm(1)
*       idoy = mmdd2doy(yyyy, mm, iday)
*       mdiday = mdiday + idoy

*       write(*,*) '--------------------------------------------'
*       yy1 = yyyy
*       mm1 = mm
*       iday = itm(1)
*       write(*,'('' Today (YYYYMMDD)='',i4,2i2.2,'' is ...'')')
*     &                                yy1,mm1,iday
*       write(*,'(i5,''th day of year '',i4)') idoy, yy1
*       write(*,'(i5,''th MDI day'')') mdiday
*       rdd = float(iday) + 0.5 - 1.0 ! MIND 1 day must be offset
*       call date2lng(cr,rlon,yy1,mm1,rdd)
*       aa = float(mm1-1) * 30.0 + rdd ! now one
*       aa = aa - 65.0 
*       aa = aa / 365.0 * 3.141592 * 2.0
*       b0 = cos(aa) * (-7.25)
*       write(*,'('' UT Noon is at CR'',i5,'' Lon'',f7.2,$)') cr, rlon
*       write(*,'('' approx B0 '',f5.2)') b0

** yesterday
*       write(*,*) '--------------------------------------------'
*       yy1 = yyyy
*       mm1 = mm
*       iday = itm(1)
*       idshift = -1
*       call shiftday(yy1,mm1,iday,idshift)
*       write(*,'('' Yesterday was (YYYYMMDD)='',i4,2i3,$)')
*     &                                yy1,mm1,iday
*       write(*,'('', '',i4,''th MDI day'')') mdiday + idshift
*       yy1 = yyyy
*       mm1 = mm
*       rdd = float(iday) + 0.5 - 1.0 ! MIND 1 day must be offset
*       call date2lng(cr,rlon,yy1,mm1,rdd)
*       aa = float(mm1-1) * 30.0 + rdd ! very very rough B0
*       aa = aa - 65.0 ! very rough offset 
*       aa = aa / 365.0 * 3.141592 * 2.0
*       b0 = cos(aa) * (-7.25)
*       write(*,'('' UT Noon is at CR'',i5,'' Lon'',f7.2,$)') cr, rlon
*       write(*,'('' approx B0 '',f5.2)') b0

** tomorrow
*       write(*,*) '--------------------------------------------'
*       yy1 = yyyy
*       mm1 = mm
*       iday = itm(1)
*       idshift = 1
*       call shiftday(yy1,mm1,iday,idshift)
*       write(*,'('' Tomorrow will be (YYYYMMDD)='',i4,2i3,$)')
*     &                                yy1,mm1,iday
*       write(*,'('', '',i4,''th MDI day'')') mdiday + idshift
*       yy1 = yyyy
*       mm1 = mm
*       rdd = float(iday) + 0.5 - 1.0 ! MIND 1 day must be offset
*       call date2lng(cr,rlon,yy1,mm1,rdd)
*       aa = float(mm1-1) * 30.0 + rdd ! very very rough B0
*       aa = aa - 65.0 ! very rough offset 
*       aa = aa / 365.0 * 3.141592 * 2.0
*       b0 = cos(aa) * (-7.25)
*       write(*,'('' UT Noon is at CR'',i5,'' Lon'',f7.2,$)') cr, rlon
*       write(*,'('' approx B0 '',f5.2)') b0

*       write(*,*) '--------------------------------------------'

       stop
       end

*
*** ---------------------------------------------
*
       subroutine date2lng(nrot,lon,yy1,mm1,dd1)
       implicit none
       integer nrot, yy1, mm1
       real    lon, dd1
       integer rt, rt1, rt2, rt0
       integer y1, m1, y2, m2, y0, m0, yy, mm
       real    d1, d2, d0, dd
       integer dyofmth
       integer month1, month2, monthm
       integer rtstart, rtend
       parameter(rtstart = -9, rtend = 2999)
       integer yystrt(rtstart:rtend)
       integer mmstrt(rtstart:rtend)
       real    ddstrt(rtstart:rtend)
       common /crrot1/ yystrt, mmstrt
       common /crrot2/ ddstrt

       yy = yy1
       mm = mm1
       dd = dd1

       nrot = -1000
       lon  = -1000.0

       do rt = rtstart + 1, rtend - 1
         rt0 = rt
         y0 = yystrt(rt0)
         m0 = mmstrt(rt0)
         d0 = ddstrt(rt0)
         if ((yy .EQ. y0) .AND. (mm .EQ. m0)) then
           if (d0 .GE. dd) then
             rt1 = rt0 - 1
             y1 = yystrt(rt1)
             m1 = mmstrt(rt1)
             d1 = ddstrt(rt1)
             rt2 = rt0
             y2  = y0
             m2  = m0
             d2  = d0
             if (m1 .EQ. m2) goto 300
           endif
           if (d0 .LT. dd) then
             rt1 = rt0
             y1  = y0
             m1  = m0
             d1  = d0
             rt2 = rt0 + 1
             y2 = yystrt(rt2)
             m2 = mmstrt(rt2)
             d2 = ddstrt(rt2)
           endif
           month1 = m1
           month2 = m2 + 12 * (y2 - y1)
           monthm = mm + 12 * (yy - y1)
           if (month1 .LT. month2) d2 = d2 + dyofmth(y1, m1)
           if (month1 .LT. monthm) dd = dd + dyofmth(y1, m1)
           lon = (1.00 - (dd - d1) / (d2 - d1)) * 360.0
           nrot = rt1
           if (lon .LE. 0.0) then
             lon = lon + 360.0
             nrot = rt2
           endif
           goto 200
         endif
 300     continue
       enddo
 200   continue

       if (nrot .LT. 0) then
         write(*,'('' ERR at date2lng 2 '',$)')
         write(*,'(2i5,f6.2,'':'',2(3i5,f6.2),'':'',i5,x,f8.2)')
     +     yy1, mm1, dd1, rt1, y1, m1, d1, rt2, y2, m2, d2, nrot, lon
       endif
       if ((lon .LT. 0.0) .OR. (lon .GT. 360.0)) then
         write(*,'('' ERR at date2lng 3 '',$)')
         write(*,'(2i5,f6.2,'':'',2(3i5,f6.2),'':'',i5,x,f8.2)')
     +     yy1, mm1, dd1, rt1, y1, m1, d1, rt2, y2, m2, d2, nrot, lon
       endif
       return
       end

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


**----------------------
*
       subroutine shiftday(yy1,mm1,idd1,idshift)
       implicit none
* interface
       integer yy1, mm1, idd1, idshift
* local
       integer dyofmth
       integer idummy

       if (abs(idshift) .GT. 28) write(*,*) ' Dshift is too big !!'

       if (idshift .GE. 0) then
         idd1 = idd1 + idshift
         idummy = dyofmth(yy1, mm1)
         if (idd1 .GT. idummy) then
           idd1 = idd1 - idummy
           mm1  = mm1 + 1
         endif
         if (mm1 .GT. 12) then
           mm1 = mm1 - 12
           yy1 = yy1 + 1
         endif
       else
         idd1 = idd1 + idshift
         if (idd1 .LE. 0) then
           mm1 = mm1 - 1
           if (mm1 .GE. 1) then
             idd1 = idd1 + dyofmth(yy1,mm1)
           else
             mm1 = 12
             yy1 = yy1 - 1
             idd1 = idd1 + dyofmth(yy1,mm1)
           endif
         endif
       endif

       return
       end


* ---------------------------------------------
* Day of Year
* ---------------------------------------------
*
       integer function mmdd2doy(yy, mm, idd)
       implicit none
       integer yy, mm, idd
       integer dyofmth
       integer mm2, idd2, idoy

       if (mm .EQ. 1) then
         idd2 = idd
       else if (mm .GT. 1) then
         idd2 = 0
         do mm2 = 1, mm - 1
           idd2 = idd2 + dyofmth(yy,mm2)
         enddo
         idd2 = idd2 + idd
       else
         write(*,*) ' MM is wrong '
       endif
       idoy = idd2
       mmdd2doy = idoy
       return
       end


*
*** ---------------------------------------------
*
       integer function dyofmth(yy, mm)
       implicit none
       integer yy, mm
       integer ans, uruu
       if ((mm .LT. 1) .OR. (mm .GT. 12))
     +    write(*,'('' Something is wrong with the Month '', 2i5)')
     +     yy, mm
       if (mm .EQ. 1) ans = 31
       if (mm .EQ. 2) then
         ans = 28
         uruu = mod(yy, 4)
         if (uruu .EQ. 0) ans = 29
       endif
       if (mm .EQ. 3) ans = 31
       if (mm .EQ. 4) ans = 30
       if (mm .EQ. 5) ans = 31
       if (mm .EQ. 6) ans = 30
       if (mm .EQ. 7) ans = 31
       if (mm .EQ. 8) ans = 31
       if (mm .EQ. 9) ans = 30
       if (mm .EQ. 10) ans = 31
       if (mm .EQ. 11) ans = 30
       if (mm .EQ. 12) ans = 31
       dyofmth = ans
       return
       end

*
*** ---------------------------------------------
*
       character*3 function mon2char(mm)
       implicit none
       integer mm
       character*3 cdummy

       if (mm .EQ.  1) cdummy = 'Jan'
       if (mm .EQ.  2) cdummy = 'Feb'
       if (mm .EQ.  3) cdummy = 'Mar'
       if (mm .EQ.  4) cdummy = 'Apr'
       if (mm .EQ.  5) cdummy = 'May'
       if (mm .EQ.  6) cdummy = 'Jun'
       if (mm .EQ.  7) cdummy = 'Jul'
       if (mm .EQ.  8) cdummy = 'Aug'
       if (mm .EQ.  9) cdummy = 'Sep'
       if (mm .EQ. 10) cdummy = 'Oct'
       if (mm .EQ. 11) cdummy = 'Nov'
       if (mm .EQ. 12) cdummy = 'Dec'
       if ((mm .LT. 1) .OR. (mm .GT. 12))
     +   write(*,*) 'MM is wrong at mm2char', mm

       mon2char = cdummy

       return
       end
*
*** ---------------------------------------------
*
* Note that 'day = 0.0' means 'just 0:00 AM of the first day of month'
*
       subroutine readdate(flname1)
       implicit none
       character*120 flname1
       integer rtstart, rtend
       parameter(rtstart = -9, rtend = 2999)
       integer yystrt(rtstart:rtend)
       integer mmstrt(rtstart:rtend)
       real    ddstrt(rtstart:rtend)
       integer nr1, y1, m1
       real    d1
       common /crrot1/ yystrt, mmstrt
       common /crrot2/ ddstrt

!       flname1 = '/home/jsoc/cvs/Development/JSOC/proj' //
!     &           '/mhd_32daily/apps/datelst4.dat'   
       open(unit = 1, file = flname1, status = 'old')
!       write(*,*) ' Now opened : ' , flname1
 100   continue
         read(1,*,END=200) nr1, y1, m1, d1
         if ((nr1 .GE. rtstart) .AND. (nr1 .LE. rtend)) then
           yystrt(nr1) = y1
           mmstrt(nr1) = m1
           ddstrt(nr1) = d1 - 1
           if (d1 - 1.0 .LT. 0.0) write(*,*) ' Wrong 2'
         endif
         goto 100
 200   continue
       close(1)
!       write(*,*) ' Now closed : ', flname1
       return
       end

*
*
* ----------------------------------------------------------------------------
* End of this file
* ----------------------------------------------------------------------------
