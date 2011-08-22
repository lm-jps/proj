   open(5,file='azimuth.dat')

! --> Recompute azimuth from ambiguity-resolved Cartesian components, 
! --> unless transverse field is zero, in which case use original angle.

   phase=0.5*float(iaflag)*pi
   do i=1,nx
      do j=1,ny
         if(bt(i,j).ne.0.) ba(i,j)=atan2(By(i,j),Bx(i,j))-phase
      enddo
   enddo

! --> Convert back to degrees, if necessary.
   if(iaunit.eq.1) then
      rtod=180./pi
      ba=ba*rtod
   endif

! --> Write out updated azimuth angle.
   do j=1,ny
      write(5,*) (ba(i,j),i=1,nx)
   enddo

