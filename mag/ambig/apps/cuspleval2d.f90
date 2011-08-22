!***********************************************************************
! subroutine cuspleval2d                                               *
!    Evaluates 2-D cubic spline interpolation with coefficients        *
! computed in cuspl2d.f                                                *
!***********************************************************************

subroutine cuspleval2d(x,y,x1,y1,f,a,fi,n1,n2,n3)

!-----------------------------------------------------------------------
   use sizes

   implicit none

   integer :: i,j,n1,n2,n3
   real :: y1,fdum
   real,dimension(n1,n2) :: f,a
   real,dimension(n1) :: x,ffdum
   real,dimension(n2) :: y,f1d,a1d
   real,dimension(n3) :: x1,fi
!-----------------------------------------------------------------------
! --> Evaluate the row splines.
   do i=1,n1
      do j=1,n2
         f1d(j)=f(i,j)
         a1d(j)=a(i,j)
      enddo
      call cuspleval(y,y1,f1d,a1d,ffdum(i),n2)
   enddo

! --> Construct and evaluate the column spline.
   call cuspl(x,ffdum,a1d,n1)

! --> Loop over values of x1, evaluating the interpolation.
   do i=1,n3
      call cuspleval(x,x1(i),ffdum,a1d,fdum,n1)
      fi(i)=fdum
   enddo

end subroutine cuspleval2d
