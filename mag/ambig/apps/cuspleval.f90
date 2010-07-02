!***********************************************************************
! subroutine cuspleval                                                 *
!    Evaluates cubic spline interpolation with coefficients computed   *
! in cuspl.f                                                           *
!***********************************************************************

subroutine cuspleval(x,x1,f,a,fi,n)

!-----------------------------------------------------------------------
   use sizes

   implicit none

   integer :: i,im,ip,n
   real :: g1,g2,dx,x1,fi
   real,dimension(n) :: f,a,x
!-----------------------------------------------------------------------
! --> Use bisection to find the indices of the points bracketing the 
! --> interpolation point. 
   im=1
   ip=n
   do while(ip-im.gt.1)
      i=(ip+im)/2
      if(x(i).gt.x1)then
         ip=i
      else
         im=i
      endif
   enddo
   dx=x(ip)-x(im)

! --> Evaluate the interpolation. 
   g1=(x(ip)-x1)/dx
   g2=(x1-x(im))/dx
   fi=g1*f(im)+g2*f(ip)+((g1**3-g1)*a(im)+(g2**3-g2)*a(ip))*(dx**2)/6.

end subroutine cuspleval
