!***********************************************************************
! subroutine mollweide                                                 *
!    Compute coordinate transforms for the Mollweide equal area        *
! projection of a sphere onto a plane.  With (x,y) the standard        *
! Cartesian coordinates in the plane, and (theta,phi) the longitude    *
! (central meridian angle) and latitude, the transform is given by     *
!     x = \frac{2 \sqrt 2}{\pi} \theta \cos\left(\lambda \right)       *
!     y = \sqrt 2 \sin\left(\lambda \right)                            *
! where \lambda is an auxiliary angle defined by                       *
!     2 \lambda + \sin(2 \lambda ) = \pi \sin(\phi)                    *
!***********************************************************************

subroutine mollweide(x,y,theta,phi)

!-----------------------------------------------------------------------
   use constant

   real :: x,y,theta,phi,lambda
!-----------------------------------------------------------------------
   if(y.gt.sqrt(2.)) then
      write(*,*) 'y coordinate outside projection range',y
         phi=-pi
         theta=-pi
      return
   endif
   lambda=asin(y/sqrt(2.))
   if(abs(2.*lambda+sin(2.*lambda))/pi.gt.1.) then
      write(*,*) 'coordinates outside projection range',x,y
         phi=-pi
         theta=-pi
      return
   endif
   phi=asin((2.*lambda+sin(2.*lambda))/pi)
   theta=pi*x/(2.*sqrt(2.)*cos(lambda))

end subroutine mollweide
