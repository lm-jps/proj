!***********************************************************************
! subroutine boxit                                                     *
!    Determine the corners of a single box which encloses all of the   *
! pixels contained in the mask.                                        *
!***********************************************************************

subroutine boxit

!-----------------------------------------------------------------------
   use sizes
   use bounds
   use maskvec
   use verbose

   implicit none

   integer :: i,j
!-----------------------------------------------------------------------
! --> Initialize.
   nxa=nx+1
   nxb=-1
   nya=ny+1
   nyb=-1
    
! --> Loop over pixels, checking which are above threshold. 
   do i=1,nx
      do j=1,ny

         if(mask(i,j).eq.1) then

            if(nxb.lt.i) nxb=i
            if(nxa.gt.i) nxa=i
            if(nyb.lt.j) nyb=j
            if(nya.gt.j) nya=j

         endif

      enddo
   enddo

! --> Expand by 1 in all directions (if possible) to allow for a 1 pixel border.
   if(nxa.gt.1) nxa=nxa-1
   if(nxb.lt.nx) nxb=nxb+1
   if(nya.gt.1) nya=nya-1
   if(nyb.lt.ny) nyb=nyb+1

! --> Determine the size of the region of interest.
   dnx=nxb-nxa+1
   dny=nyb-nya+1

   if(iverb.eq.2) write(*,*) 'box:',nxa,nxb,nya,nyb

end subroutine boxit
