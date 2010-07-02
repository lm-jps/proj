!***********************************************************************
! subroutine transform                                                 *
!    Construct coordinate transformation matrices from Gary and        *
! Hagyard (1990).                                                      *
! Note that the matrix C is the transpose of A, which is also the      *
! inverse of A, since A is orthogonal.                                 *
!***********************************************************************

subroutine transform(theta,phi)

!-----------------------------------------------------------------------
   use ephemeris
   use trnsfrm
!-----------------------------------------------------------------------
   real :: theta,phi
!-----------------------------------------------------------------------

   c11=-sin(b)*sin(p0)*sin(theta) + cos(p0)*cos(theta)
   c12=-sin(phi)*(sin(b)*sin(p0)*cos(theta) + cos(p0)*sin(theta))-cos(phi)*cos(b)*sin(p0)
   c21= sin(b)*cos(p0)*sin(theta) + sin(p0)*cos(theta)
   c22= sin(phi)*(sin(b)*cos(p0)*cos(theta) - sin(p0)*sin(theta))+cos(phi)*cos(b)*cos(p0)

   a13=-cos(b)*sin(theta)
   a23=-cos(b)*sin(phi)*cos(theta) + sin(b)*cos(phi)
   a31= cos(phi)*(sin(b)*sin(p0)*cos(theta) + cos(p0)*sin(theta))-sin(phi)*cos(b)*sin(p0)
   a32=-cos(phi)*(sin(b)*cos(p0)*cos(theta) - sin(p0)*sin(theta))+sin(phi)*cos(b)*cos(p0)
   a33= cos(phi)*cos(b)*cos(theta) + sin(phi)*sin(b)

end subroutine transform
