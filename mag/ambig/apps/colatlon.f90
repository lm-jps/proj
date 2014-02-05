!***********************************************************************
! subroutine colatlon                                                  *
!    Generates arrays containing colatitude and longitude of each      *
! pixel, as well as sine and cosine of colatitude and longitude.       *
!***********************************************************************

subroutine colatlon()

!-----------------------------------------------------------------------
   use sizes
   use constant
   use disk_center
   use ephemeris
   use pix_size
   use spherical_position

   implicit none

   integer :: i,j
   real :: xx,yy,y2
!-----------------------------------------------------------------------
! --> Allocate memory for colat/lon arrays.
   allocate(t(nx,ny),p(ny))
   allocate(sint(nx,ny),cost(nx,ny),sinp(ny),cosp(ny))

! --> Generate colat/lon for each pixel.
   do j=1,ny
      yy=(float(j)*dyi-ycen)/radius

! --> Check this pixel is between the north and south poles.
      if(abs(yy).lt.1.) then
         y2=yy*yy
         p(j)=0.5*pi-asin(yy)
         cosp(j)=yy
         sinp(j)=sqrt(1.-y2)
         do i=1,nx
            xx=(float(i)*dxi-xcen)/radius

! --> Check this pixel is on the disk.
            if(xx*xx+y2.lt.0.999) then
               t(i,j)=asin(xx/sqrt(1.-y2))
               sint(i,j)=xx/sqrt(1.-y2)
               cost(i,j)=sqrt(1.-sint(i,j)**2)
            else
               t(i,j)=-999.
               sint(i,j)=-999.
               cost(i,j)=-999.
            endif
         enddo
      else
         p(j)=-999.
         sinp(j)=-999.
         cosp(j)=-999.
         do i=1,nx
            t(i,j)=-999.
            sint(i,j)=-999.
            cost(i,j)=-999.
         enddo
      endif
   enddo

end subroutine colatlon

