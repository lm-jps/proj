!***********************************************************************
! subroutine cuspl                                                     *
!    Computes one dimensional cubic spline coefficients.               *
!***********************************************************************

subroutine cuspl(x,f,a,n)

!-----------------------------------------------------------------------
   use sizes

   implicit none

   integer :: i,n
   real :: dxp,dxm
   real,dimension(n) :: x,f,a,g
!-----------------------------------------------------------------------
! --> Natural boundary condition. 
   a(1)=0.
   g(1)=0.
! --> Require first derivative to be 0. at the first boundary. 
!   a(1)=-0.5
!   g(1)=3.*(f(2)-f(1))/(x(2)-x(1))**2
   do i=2,n-1
      dxp=x(i+1)-x(i)
      dxm=x(i)-x(i-1)
      a(i)=-dxp/(a(i-1)*dxm+2.*(dxm+dxp))
      g(i)=(6.*((f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm)-dxm*g(i-1))/(a(i-1)*dxm+2.*(dxm+dxp))
   enddo

! --> Require first derivative to be 0. at the second boundary. 
   g(n)=-3.*(f(n)-f(n-1))/(x(n)-x(1))**2
   a(n)=(2.*g(n)-g(n-1))/(a(n-1)+2.)
! --> Natural boundary condition.
!   a(n)=0.
   do i=n-1,1,-1
      a(i)=a(i)*a(i+1)+g(i)
   enddo

end subroutine cuspl
