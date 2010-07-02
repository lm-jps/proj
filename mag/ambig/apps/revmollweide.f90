!***********************************************************************
! subroutine revmollweide                                              *
!    Compute reverse coordinate transforms for the Mollweide equal     *
! area projection of a sphere onto a plane.  With (x,y) the standard   *
! Cartesian coordinates in the plane, and (theta,phi) the longitude    *
! (central meridian angle) and latitude, the (forward) transform is    *
! given by                                                             *
!     x = \frac{2 \sqrt 2}{\pi} \theta \cos\left(\lambda \right)       *
!     y = \sqrt 2 \sin\left(\lambda \right)                            *
! where \lambda is an auxiliary angle defined by                       *
!     2 \lambda + \sin(2 \lambda ) = \pi \sin(\phi)                    *
!***********************************************************************

subroutine revmollweide(theta,phi,x,y)

!-----------------------------------------------------------------------
   use constant

   integer :: iter,maxiter
   real :: x,y,theta,phi,lambda,lambdap,dlambdap
!-----------------------------------------------------------------------
   lambdap=phi
   dlambdap=1.
   iter=0
   maxiter=50
   do while(abs(dlambdap).gt.1.e-4.and.iter.le.maxiter)
      dlambdap=-(lambdap+sin(lambdap)-pi*sin(phi))/(1.+cos(lambdap))
      lambdap=lambdap+dlambdap
      iter=iter+1
   enddo
   lambda=lambdap/2.
   x=2.*sqrt(2.)*theta*cos(lambda)/pi
   y=sqrt(2.)*sin(lambda)

end subroutine revmollweide
