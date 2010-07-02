!***********************************************************************
! subroutine grown                                                     *
!    Given an image and a distance, construct an array labelling which *
! pixels are either themselves above the threshold or have a neighbor  *
! above the threshold which is within the given distance.              *
!***********************************************************************

subroutine grown(n)

!-----------------------------------------------------------------------
   use sizes
   use maskvec

   implicit none

   integer :: i,j,k,l,n
   integer(kind=1),dimension(:,:),allocatable :: imagepad
!   real,dimension(nx,ny) :: image
!-----------------------------------------------------------------------
!open(4,file='mask.dat')
! --> If distance is not positive definite, return the original mask. 
   if(n.gt.0) then

! --> Allocate memory for a padded copy of the image.
      allocate(imagepad(nx+2*n,ny+2*n))
      imagepad=0

! --> Grow the region by shifting the binary image by n pixels in each 
! --> direction. 
      do i=1,nx
         do j=1,ny
            if(mask(i,j).eq.1) then
               do k=0,2*n
                  do l=0,2*n
                     imagepad(i+k,j+l)=1
                  enddo
               enddo
            endif
         enddo
      enddo

! --> Replace input image with grown image.

      do i=1,nx
         do j=1,ny
            mask(i,j)=imagepad(i+n,j+n)
         enddo
      enddo

! --> Deallocate memory.
      deallocate(imagepad)

   endif
!do j=1,ny
!   write(4,*) (mask(i,j),i=1,nx)
!enddo

end subroutine grown

