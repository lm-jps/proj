!***********************************************************************
! subroutine buffer                                                    *
!    Pad the input array with zeros.                                   *
!***********************************************************************

subroutine buffer(bl,nxp,nyp)

!-----------------------------------------------------------------------
   use sizes
   use pad

   implicit none

   integer :: i,j,nxp,nyp
   real :: bl(nx,ny)
!-----------------------------------------------------------------------
! --> non-zero data will be centered

   ixmin=npad+1
   ixmax=ixmin+nx-1
   iymin=npad+1
   iymax=iymin+ny-1

   nxp=2*npad+nx
   nyp=2*npad+ny

! --> Allocate memory.
   allocate(blpad(nxp,nyp))

   do i=1,ixmin-1
      do j=iymin,iymax
         blpad(i,j)=0.
      enddo
   enddo

   do i=ixmax+1,nxp
      do j=iymin,iymax
         blpad(i,j)=0.
      enddo
   enddo

   do i=1,nxp
      do j=1,iymin-1
         blpad(i,j)=0.
      enddo
   enddo

   do i=1,nxp
      do j=iymax+1,nyp
         blpad(i,j)=0.
      enddo
   enddo

   do i=1,nx
      do j=1,ny
         blpad(i+ixmin-1,j+iymin-1)=bl(i,j)
      enddo
   enddo

end subroutine buffer
