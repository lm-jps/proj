    module polyval_f
      integer, parameter :: d=8
      contains
      REAL(kind=d) FUNCTION POLYVAL(n,c,x)
!
!     Evaluate the polynomial
!
!         c(n+1)*x^n + c(n)*x^(n-1) ... c(2)*x + c(1)
!
!     using Horner's scheme.
!
      implicit none
      integer, parameter :: double=8
      integer :: n,i
      real(kind=double), dimension(n+1) :: c
      real(kind=double) :: x

      polyval = c(n+1)
      do i = n,1,-1
         polyval = (x*polyval + c(i))
      end do
      return
      end function
    end module
