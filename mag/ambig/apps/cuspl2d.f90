!***********************************************************************
! subroutine cuspl2d                                                   *
!    Computes 2-D spline coefficients.                                 *
!***********************************************************************

subroutine cuspl2d(y,f,a,n1,n2)

!-----------------------------------------------------------------------
   use sizes

   implicit none

   integer :: i,j,n1,n2
   !real,dimension(nmax,nmax) :: f,a
   !real,dimension(nmax) :: y
   real,dimension(n1,n2) :: f,a
   real,dimension(n2) :: y
   real,dimension(n2) :: fdum,adum
!-----------------------------------------------------------------------
   do i=1,n1
      do j=1,n2
         fdum(j)=f(i,j)
      enddo

      call cuspl(y,fdum,adum,n2)

      do j=1,n2
         a(i,j)=adum(j)
      enddo
   enddo

end subroutine cuspl2d
