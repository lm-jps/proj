!***********************************************************************
! subroutine tukey                                                     *
!    Apply a Tukey windowing function to reduce ringing in the         *
! potential field calculation.                                         *
!***********************************************************************

subroutine tukey(b,nx,ny,nap)

!-----------------------------------------------------------------------
   use constant
   use pad

   implicit none

   real :: b(nx,ny)
   integer,intent(in) :: nx,ny,nap
   integer :: i,j
!-----------------------------------------------------------------------
! --> Allocate memory for output array
   allocate(blpad(nx,ny))

   do i=1,nap-1
      do j=1,nap-1
         blpad(i,j)=b(i,j)*0.25*(1.-cos(float(i-1)/float(nap-1)*pi))&
                               *(1.-cos(float(j-1)/float(nap-1)*pi))
      enddo

      do j=nap,ny-nap
         blpad(i,j)=b(i,j)*0.5*(1.-cos(float(i-1)/float(nap-1)*pi))
      enddo

      do j=ny-nap+1,ny
         blpad(i,j)=b(i,j)*0.25*(1.-cos(float(i-1)/float(nap-1)*pi))&
                               *(1.-cos(float(ny-j)/float(nap-1)*pi))
      enddo
   enddo

   do i=nap,nx-nap
      do j=1,nap-1
         blpad(i,j)=b(i,j)*0.5*(1.-cos(float(j-1)/float(nap-1)*pi))
      enddo

      do j=nap,ny-nap
         blpad(i,j)=b(i,j)
      enddo

      do j=ny-nap+1,ny
         blpad(i,j)=b(i,j)*0.5*(1.-cos(float(ny-j)/float(nap-1)*pi))
      enddo
   enddo

   do i=nx-nap+1,nx
      do j=1,nap-1
         blpad(i,j)=b(i,j)*0.25*(1.-cos(float(nx-i)/float(nap-1)*pi))&
                               *(1.-cos(float(j-1)/float(nap-1)*pi))
      enddo

      do j=nap,ny-nap
         blpad(i,j)=b(i,j)*0.5*(1.-cos(float(nx-i)/float(nap-1)*pi))
      enddo

      do j=ny-nap+1,ny
         blpad(i,j)=b(i,j)*0.25*(1.-cos(float(nx-i)/float(nap-1)*pi))&
                               *(1.-cos(float(ny-j)/float(nap-1)*pi))
      enddo
   enddo

end subroutine tukey

