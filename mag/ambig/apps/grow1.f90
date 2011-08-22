!***********************************************************************
! subroutine grow1                                                     *
!    Given an image and a threshold, construct an array labelling      *
! which pixels are either themselves above the threshold or have a     *
! neighbor which is above the threshold.                               *
!***********************************************************************

subroutine grow1(image,bthresh)

!-----------------------------------------------------------------------
   use sizes
   use maskvec

   implicit none

   integer :: i,j
   real :: bthresh
   real,dimension(nx,ny) :: image
!-----------------------------------------------------------------------
! --> Initialize mask array to be 1 for pixels above threshold.
   do i=1,nx
      do j=1,ny
         if(image(i,j).ge.bthresh) then
            mask(i,j)=1
         else
            mask(i,j)=0
         endif
      enddo
   enddo

! --> Grow the region by shifting the binary image by one pixel in each 
! --> direction. Need to treat edges and corners separately, since not
! --> all directions are available.

! --> Interior pixels.
   do i=2,nx-1
      do j=2,ny-1
         if(image(i,j).ge.bthresh) then
            mask(i,j-1)=1
            mask(i,j+1)=1
            mask(i-1,j)=1
            mask(i+1,j)=1
            mask(i-1,j-1)=1
            mask(i+1,j-1)=1
            mask(i-1,j+1)=1
            mask(i+1,j+1)=1
         endif
      enddo
   enddo

! --> Top and bottom edges.
   do i=1,nx-1
      if(image(i,1).ge.bthresh) then
         mask(i,2)=1
         mask(i-1,1)=1
         mask(i+1,1)=1
         mask(i-1,2)=1
         mask(i+1,2)=1
      endif
      if(image(i,ny).ge.bthresh) then
         mask(i,ny-1)=1
         mask(i-1,ny)=1
         mask(i+1,ny)=1
         mask(i-1,ny-1)=1
         mask(i+1,ny-1)=1
      endif
   enddo

! --> Left and right edges.
   do j=2,ny-1
      if(image(1,j).ge.bthresh) then
         mask(2,j)=1
         mask(1,j-1)=1
         mask(1,j+1)=1
         mask(2,j-1)=1
         mask(2,j+1)=1
      endif
      if(image(nx,j).ge.bthresh) then
         mask(nx-1,j)=1
         mask(nx,j-1)=1
         mask(nx,j+1)=1
         mask(nx-1,j-1)=1
         mask(nx-1,j+1)=1
      endif
   enddo

! --> Four corners.
   if(image(1,1).ge.bthresh) then
      mask(1,2)=1
      mask(2,2)=1
      mask(2,1)=1
   endif
   if(image(1,ny).ge.bthresh) then
      mask(1,ny-1)=1
      mask(2,ny-1)=1
      mask(2,ny)=1
   endif
   if(image(nx,1).ge.bthresh) then
      mask(nx,2)=1
      mask(nx-1,2)=1
      mask(nx-1,1)=1
   endif
   if(image(nx,ny).ge.bthresh) then
      mask(nx,ny-1)=1
      mask(nx-1,ny-1)=1
      mask(nx-1,ny)=1
   endif

end subroutine grow1

