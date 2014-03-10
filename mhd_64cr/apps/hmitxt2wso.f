*
* make dummy STN data
*
      program hmitxt2wso
      implicit none
      character*48 flname
      integer i, j, ilon, jlat, cr
      real*8  aa
      real*8  magwso(72,30)
      real*8  magwpart(30)

* this does not be used .....
      open(unit=11,file='hmi2wso.txt',status='old')
      do i = 1, 72
      do j = 1, 30
        read(11,*) ilon, jlat, aa
        ilon = ilon + 1
        jlat = jlat + 1
        magwso(ilon,jlat) = aa ! at this point, latitude must run North to South and unit be in micro T
      enddo
      enddo
      close(11)

* write out in WSO format
      cr = 1000 ! dummy CR number
      write(flname,'(''stn'',i4,''.dat'')') cr
      open(unit = 12, file = flname, status = 'unknown')
      write(*,*) ' Now open file : ', flname

      do i = 72, 1, -1
        ilon = i * 5
        do jlat = 1, 30
          magwpart(jlat) = magwso(i,jlat)
        enddo
        write(12,'(''CT'',i4,'':'',i3.3,''        '',6f9.3)')
     +    cr, ilon, (magwpart(jlat),jlat=1,6)
        write(12,'(8f9.3)') (magwpart(jlat),jlat=  7,14)
        write(12,'(8f9.3)') (magwpart(jlat),jlat= 15,22)
        write(12,'(8f9.3)') (magwpart(jlat),jlat= 23,30)
      enddo
      close(12)

      stop
      end
