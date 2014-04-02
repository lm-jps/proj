*
       program u3_to_d3
       implicit none
       integer n, rt, ihead
       integer i, j, k
       character*48 flname
       logical ldummy
* number of grid for each three direction.
       integer ii, jj, kk
       parameter(ii = 72, jj = 32, kk = 64)
*       real    ro1(-2:jj+2,-1:kk+2,-2:ii+3),pg1(-2:jj+2,-1:kk+2,-2:ii+3)! single
*       real    ur1(-2:jj+2,-1:kk+2,-2:ii+3),ut1(-2:jj+2,-1:kk+2,-2:ii+3)
*       real    up1(-2:jj+2,-1:kk+2,-2:ii+3),br1(-2:jj+2,-1:kk+2,-2:ii+3)
*       real    bt1(-2:jj+2,-1:kk+2,-2:ii+3),bp1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  ro1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  pg1(-2:jj+2,-1:kk+2,-2:ii+3)! double
       real*8  ur1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  ut1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  up1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  br1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  bt1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  bp1(-2:jj+2,-1:kk+2,-2:ii+3)
       real*8  ro2(0:ii+1,0:jj,0:kk+1),pg2(0:ii+1,0:jj,0:kk+1)
       real*8  ur2(0:ii+1,0:jj,0:kk+1),ut2(0:ii+1,0:jj,0:kk+1)
       real*8  up2(0:ii+1,0:jj,0:kk+1),br2(0:ii+1,0:jj,0:kk+1)
       real*8  bt2(0:ii+1,0:jj,0:kk+1),bp2(0:ii+1,0:jj,0:kk+1)

* initialize all
       do i = -2, ii + 3
         do k = -1, kk + 2
           do j = -2, jj + 2
             ro1(j,k,i) = 0.0D+00
             pg1(j,k,i) = 0.0D+00
             ur1(j,k,i) = 0.0D+00
             ut1(j,k,i) = 0.0D+00
             up1(j,k,i) = 0.0D+00
             br1(j,k,i) = 0.0D+00
             bt1(j,k,i) = 0.0D+00
             bp1(j,k,i) = 0.0D+00
           enddo
         enddo
       enddo
       do k = 0, kk + 1
         do j = 0, jj
           do i = 0, ii + 1
             ro2(i,j,k) = 0.0D+00
             pg2(i,j,k) = 0.0D+00
             ur2(i,j,k) = 0.0D+00
             ut2(i,j,k) = 0.0D+00
             up2(i,j,k) = 0.0D+00
             br2(i,j,k) = 0.0D+00
             bt2(i,j,k) = 0.0D+00
             bp2(i,j,k) = 0.0D+00
           enddo
         enddo
       enddo

 100   continue

 111     continue
         write(*,*) ' INPUT 3 integers for CR Nstep'
         write(*,*) '     and Index of InputData(1-4=uxvy)'
         write(*,*) ' Any negative integer for quit'
         read(*,*) rt, n, ihead
         if (ihead .GT. 4) goto 111
         if ((rt .LT. 0) .OR. (n .LT. 0) .OR. (ihead .LT. 0)) goto 999

         if ((ihead .EQ. 1) .OR. (ihead .EQ. 3)) then
           if (ihead .EQ. 1) then
             write(flname,'(''u'',i6,''.'',i4)') n + 300000, rt
           else if (ihead .EQ. 3) then
             write(flname,'(''v'',i6,''.'',i4)') n + 300000, rt
           endif
           inquire(file = flname, exist = ldummy)
           if (ldummy) then
             write(*,*) ' File found : ',flname
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

             do k = 0, kk + 1
             do j = 0, jj
             do i = 0, ii + 1
               ro1(j,k,i) = ro2(i,j,k)
               pg1(j,k,i) = pg2(i,j,k)
               ur1(j,k,i) = ur2(i,j,k)
               ut1(j,k,i) = ut2(i,j,k)
               up1(j,k,i) = up2(i,j,k)
               br1(j,k,i) = br2(i,j,k)
               bt1(j,k,i) = bt2(i,j,k)
               bp1(j,k,i) = bp2(i,j,k)
             enddo
             enddo
             enddo

             if (ihead .EQ. 1) then
               write(flname,'(''x'',i6,''.'',i4)') n + 300000, rt
             else if (ihead .EQ. 3) then
               write(flname,'(''y'',i6,''.'',i4)') n + 300000, rt
             endif
             write(*,*) ' Output file  : ',flname
             open(unit=2,file=flname,
     +            status='unknown',form='unformatted')
             write(2) ro1
             write(2) pg1
             write(2) ur1
             write(2) ut1
             write(2) up1
             write(2) br1
             write(2) bt1
             write(2) bp1
             close(2)
           else
             write(*,*) ' Cannot find file : ', flname
           endif

         else if ((ihead .EQ. 2) .OR. (ihead .EQ. 4)) then

           if (ihead .EQ. 2) then
             write(flname,'(''x'',i6,''.'',i4)') n + 300000, rt
           else if (ihead .EQ. 4) then
             write(flname,'(''y'',i6,''.'',i4)') n + 300000, rt
           endif
           inquire(file = flname, exist = ldummy)

           if (ldummy) then
             write(*,*) ' File found : ',flname
             open(unit=1,file=flname,status='old',form='unformatted')
             read(1) ro1
             read(1) pg1
             read(1) ur1
             read(1) ut1
             read(1) up1
             read(1) br1
             read(1) bt1
             read(1) bp1

             do k = 0, kk + 1
             do j = 0, jj
             do i = 0, ii + 1
               ro2(i,j,k) = ro1(j,k,i)
               pg2(i,j,k) = pg1(j,k,i)
               ur2(i,j,k) = ur1(j,k,i)
               ut2(i,j,k) = ut1(j,k,i)
               up2(i,j,k) = up1(j,k,i)
               br2(i,j,k) = br1(j,k,i)
               bt2(i,j,k) = bt1(j,k,i)
               bp2(i,j,k) = bp1(j,k,i)
             enddo
             enddo
             enddo
  
             if (ihead .EQ. 2) then
               write(flname,'(''u'',i6,''.'',i4)') n + 300000, rt
             else if (ihead .EQ. 4) then
               write(flname,'(''v'',i6,''.'',i4)') n + 300000, rt
             endif
             write(*,*) ' Output file  : ',flname
             open(unit=2,file=flname,
     +            status='unknown',form='unformatted')
             write(2) ro2
             write(2) pg2
             write(2) ur2
             write(2) ut2
             write(2) up2
             write(2) br2
             write(2) bt2
             write(2) bp2
             close(2)

           else
             write(*,*) ' Cannot find file : ', flname
           endif

         endif

       goto 100
 999   continue

       stop
       end


