*-------------------------------------------------------------------------
*  this program will read fits of MDI/EOF and MDI/Xudong's map,
*  combine the two, then
*  make plain-text file
*-------------------------------------------------------------------------

       program hmi2txt
       implicit none
*       real    bpoles(0:3599,0:62,0:1) ! xudong's polar field
       integer ieof, jeof
       parameter(ieof = 3600, jeof = 1440) ! sinlat ?
       real    hmimag(0:ieof-1,0:jeof-1)
       real    hmimag2(0:ieof-1,0:jeof-1)
       integer iwso, jwso
       parameter(iwso = 72, jwso = 30) ! sinlat.
       real    brn2(0:iwso-1,0:1), brs2(0:iwso-1,0:1)
       real    magwso(0:iwso-1,0:jwso-1)
*       character*80 ftshead(200) ! 200 is chosen as a large number.
       character*24 fitsname
       real    llongio, llongi
       integer ncar, slongi
       integer i, j, k, j2, i2, k2, j3, k3
       integer is, ie, js, je
       integer nn
       real    aa, bb, cc, dd, zz, sinlat

* read renamed MDI/EOF file
       write(*,*) 'Input FITS name (HMI synoptic)'
       read(*,*) fitsname

       call readhmisynop(fitsname,hmimag,llongio)

       ncar  = int(llongio/360.0) + 2
       llongi=llongio-int(llongio/360.0)*360.0
       slongi=int(llongi) ! this is longitude at the leftmost part from lon=360

       open(unit=11,file='carr.txt',status='unknown')
       write(11,*) ncar !              carrington rot.num
       write(11,*) 360 - slongi ! + 60 ! carrington long. non-NRT
       close(11)
*       write(*,*) llongi, slongi

* shift map so that new one will have central longitude of 180 deg in Carrington system.
       do j = 0, jeof-1
       do i = 0, ieof-1
         hmimag2(i,j) = hmimag(i,j)
       enddo
       enddo
       do j = 0, jeof-1
       do i = 0, ieof-1
         aa = float(i) + float(slongi) * float(ieof) / 360.0
         i2 = int(aa)
         if (i2 .LT.    0) i2 = i2 + ieof
         if (i2 .GE. ieof) i2 = i2 - ieof
         hmimag(i,j) = hmimag2(i2,j)
       enddo
       enddo

* this part assume iwso = ieof/10, and jwso = jeof/12
       is =  - (ieof / iwso) / 2 + 1
       ie =    (ieof / iwso) / 2
       js =    0
       je =    (jeof / jwso) - 1
       do j = 0, jwso - 1
       do k = 0, iwso - 1
         aa = 0.0
         nn = 0
         do j2 =  js, je
         do k2 =  is, ie
           k3 = (float(k)-0.5)*(ieof / iwso)+k2 ! half-grid shift
           j3 =              j*(jeof / jwso)+j2
           if (k3 .GE. ieof) k3 = k3 - ieof
           if (k3 .LT.    0) k3 = k3 + ieof
           bb = hmimag(k3,j3)
           if ((bb .GE. -3000.0) .AND. (bb .LE. 3000.0)) then ! exclude NaN and some abnormal pixles
             aa = aa + bb
             nn = nn + 1
           endif
         enddo
         enddo
         if (nn .GT. 0) then
           magwso(k,j) = aa / float(nn)
         else
           magwso(k,j) = 0.0
         endif
       enddo
       enddo

* gauss -> micro T
       do j = 0, jwso-1
       do k = 0, iwso-1
         magwso(k,j) =  magwso(k,j) * 100.0
         if (magwso(k,j) .GE.  9999.0) magwso(k,j) = 9999.0
         if (magwso(k,j) .LE. -9999.0) magwso(k,j) =-9999.0
       enddo
       enddo

* convert Br to Blos by multiplying sin(theta)
       do j = 0, jwso-1
       do k = 0, iwso-1
         zz =  1.0 - 2.0 * (float(j)+0.5) / float(jwso)
         sinlat = sqrt(1.0 - zz*zz)
         magwso(k,j) =  magwso(k,j) * sinlat
       enddo
       enddo

* writeout text-ized file
       open(unit=12,file='hminew.txt',status='unknown')
       do k = 0, ieof-1
       do j = 0, jeof-1
         write(12,*) k,j,hmimag(k,(jeof-1)-j)
       enddo
       enddo
       close(12)

* add tweak.... until authentic polar field correction available...
       do k = 0, iwso-1
         magwso(k,     0) = magwso(k,     1)
         magwso(k,jwso-1) = magwso(k,jwso-2)
       enddo

       open(unit=13,file='hmi2wso.txt',status='unknown')
       do k = 0, iwso-1
       do j = 0, jwso-1
         write(13,*) k,j,magwso(k,(jwso-1)-j)
       enddo
       enddo
       close(13)

       stop
       end
*
**  end of main layer  --------------------------------------------------------
*

****  --------------------------------------------------------
* read synoptic map
****  --------------------------------------------------------
*
       subroutine readhmisynop(fitsname,hmimag,llongio)
       implicit none
       character*24 fitsname
*       character*80 ftshead(200) ! 200 is chosen as a large number.
       integer ieof, jeof
       parameter(ieof = 3600, jeof = 1440) ! sinlat ?
       real    hmimag(0:ieof-1,0:jeof-1)
       real    fitdat(ieof,jeof)
       real    llongio
*
       integer istatus,iunit,ireadwrite,blocksize
       integer nfound, naxes(2) ! data will be 2-D
       integer group,firstpix,npixels ! , nbuffer
       integer j, k
       real    datamin,datamax,nullval ! ,buffer(100)
       logical anynull
       character*72 hcomment
       character*21 partcmnt
       character*24 filename

       filename = fitsname
       write(*,*) 'READ_FITS : '//filename
       istatus=0
* Get an unused Logical Unit Number to use to open the FITS file.
       call ftgiou(iunit,istatus)

* open the FITS file previously created by WRITE_IMAGE.
       ireadwrite=0 ! to read
       call ftopen(iunit,filename,ireadwrite,blocksize,istatus)
* determine the size of the image
       call ftgknj(iunit,'NAXIS',1,2,naxes,nfound,istatus)
* check that it found both NAXIS1 and NAXIS2 keywords
       if (nfound .NE. 2)then
         print *,'READ_IMAGE failed to read the NAXISn keywords.'
         call ftclos(iunit, istatus)
         call ftfiou(iunit, istatus)
         return
       endif
* check if the sizes are as expected or not
       if ((naxes(1) .NE. ieof) .OR. (naxes(2) .NE. jeof)) then
         print *, 'Image size is not expected'
         call ftclos(iunit, istatus)
         call ftfiou(iunit, istatus)
         return
       endif

* get keyword
       call ftgkey(iunit,'LON_FRST',partcmnt,hcomment,istatus)
       read(partcmnt,*) llongio
*       llongio = 360.0
       write(*,*) 'LON_FRST=',llongio

* initialize variables.
       npixels=naxes(1)*naxes(2)
       group=1
       firstpix=1
       nullval= 0 ! leave NaN as NaN
*       nullval  = -999  ! substitute NaN to some number
       datamin=1.0E30
       datamax=-1.0E30

       call ftgpve(iunit,group,firstpix,npixels,
     &             nullval,fitdat,anynull,istatus)
       do j = 0, jeof - 1
       do k = 0, ieof - 1
         hmimag(k,j) = fitdat(k+1,j+1)
         if ((fitdat(k+1,j+1) .LT. 10.0) .OR.
     &       (fitdat(k+1,j+1) .GT. -10.0)) then ! exclude NaN
           datamin=min(datamin,fitdat(k+1,j+1))
           datamax=max(datamax,fitdat(k+1,j+1))
         endif
       enddo
       enddo

       print *,'Min and max values in the image are:',datamin,datamax
* close the file and free the unit number
       call ftclos(iunit, istatus)
       call ftfiou(iunit, istatus)
* Check for any error, and if so print out error messages
*       if (istatus .GT. 0) call print_error(istatus)

       return
       end

*
**  ---------------
*    reading polar field data, currently not used: as of 2014 Feb.
**  ---------------
*
       subroutine readpole(fitsname,bpoles)
       implicit none
       character*24 fitsname
       real    bpoles(0:3599,0:62,0:1) ! xudong's polar field
*
       real    fitdat(3600,63,2)
       integer istatus,iunit,ireadwrite,blocksize
       integer nfound, naxes(3) ! data will be 3-D
       integer group,firstpix,npixels !,nbuffer
       integer i, j, k
       real    datamin,datamax,nullval ! ,buffer(100)
       logical anynull
       character*24 filename
*
       filename = fitsname
       write(*,*) 'READ_FITS : '//filename
       istatus=0
* Get an unused Logical Unit Number to use to open the FITS file.
       call ftgiou(iunit,istatus)

* open the FITS file previously created by WRITE_IMAGE.
       ireadwrite=0 ! to read
       call ftopen(iunit,filename,ireadwrite,blocksize,istatus)
* determine the size of the image
       call ftgknj(iunit,'NAXIS',1,3,naxes,nfound,istatus)
* check that it found NAXIS1---NAXIS3 keywords
       if (nfound .NE. 3)then
         print *,'READ_IMAGE failed to read the NAXISn keywords.'
         call ftclos(iunit, istatus)
         call ftfiou(iunit, istatus)
         return
       endif

* check if the sizes are as expected or not
       if ((naxes(1) .NE. 3600) .OR.
     &     (naxes(2) .NE.   63) .OR.
     &     (naxes(3) .NE.    2)) then
         print *, 'Image size is not expected'
         call ftclos(iunit, istatus)
         call ftfiou(iunit, istatus)
         return
       endif

* initialize variables.
       npixels=naxes(1)*naxes(2)*naxes(3)
       group=1
       firstpix=1
       nullval= 0 ! leave NaN as NaN
*       nullval  = -999  ! substitute NaN to some number
       datamin=1.0E30
       datamax=-1.0E30

       call ftgpve(iunit,group,firstpix,npixels,
     &             nullval,fitdat,anynull,istatus)
       do i = 0, 1
       do j = 0, 62
       do k = 0, 3599
         bpoles(k,j,i) = fitdat(k+1,j+1,i+1)
         if ((fitdat(k+1,j+1,i+1) .LT. 10.0) .OR.
     &       (fitdat(k+1,j+1,i+1) .GT. -10.0)) then ! exclude NaN
           datamin=min(datamin,fitdat(k+1,j+1,i+1))
           datamax=max(datamax,fitdat(k+1,j+1,i+1))
         endif
       enddo
       enddo
       enddo

* read out with buffer
*       do while (npixels .GT. 0)
*         nbuffer=min(100,npixels)
*         call ftgpve(iunit,group,firstpix,nbuffer,
*     &               nullval,buffer,anynull,istatus)
*         do i=1,nbuffer
*           datamin=min(datamin,buffer(i))
*           datamax=max(datamax,buffer(i))
*         enddo
*         npixels=npixels-nbuffer
*         firstpix=firstpix+nbuffer
*       enddo

       print *,'Min and max values in the image are:',datamin,datamax
* close the file and free the unit number
       call ftclos(iunit, istatus)
       call ftfiou(iunit, istatus)
* Check for any error, and if so print out error messages
*       if (istatus .GT. 0) call print_error(istatus)

       return
       end

*
**  end of this file --------------------------------------------------------
*
