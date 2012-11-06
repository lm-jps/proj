!***********************************************************************************************************************************
subroutine setup_spherical_PF_4p
!===================================================================================================================================
! This subroutine computes the coefficients used to approximate derivatives in the theta- and phi-directions given the magnetic
! field Bxi, Byi, Bzi and the radial derivative of the radial component of the field, r*dBrdr, at discrete locations on a spherical
! surface. All derivatives are approximated at a point which is located between the centres of four neighbouring pixels: 1=(i,j),
! 2=(i+1,j), 3=(i,j+1), and 4=(i+1,j+1). It is assumed that the positions theta and phi are supplied for each pixel (arrays t and p)
! and that phi is a function of yi only (this assumption could be relaxed but it will change the coefficients for the derivatives). 
! It is further assumed that r*dBr/dr is supplied at one height at the same location as Bxi, Byi, Bzi and that r*dBr/dr is not
! affected by the ambiguity.
! Note that r*dBr/dr is stored in the array dBzdz.
!
! 2009, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use bounds
   use maskvec
   use mgram_data
   use spherical_deriv_coefficients
   use spherical_position

   implicit none

   integer :: i,j,k
!
! Determine the size of the region of interest
!
   write (*,*) nxa, nxb
   dnx=nxb-nxa+1
   dny=nyb-nya+1
   if (dnx.lt.2) then
      print*,'setup: dnx too small'
      stop
   endif
   if (dny.lt.2) then
      print*,'setup: dny too small'
      stop
   endif
   if (dnx.gt.nx) then
      print*,'setup: dnx too big'
      stop
   endif
   if (dny.gt.ny) then
      print*,'setup: dny too big'
      stop
   endif
!
! Allocate memory for derivative coefficient arrays, and zero the arrays.
!
   allocate(ddt1(nx,ny),ddt3(nx,ny))
   allocate(ddp(nx,ny,4),ddr(nx,ny,4))
   ddt1=0.
   ddt3=0.
   ddp=0.
   ddr=0.

!
! Coefficients for the derivative d/dtheta (multiplied by 2)
!
   do i=1,(nx-1)
      do j=1,(ny-1)
         if(mask(i,j).ge.1.and.mask(i+1,j).ge.1.and.mask(i,j+1).ge.1.and.mask(i+1,j+1).ge.1) then
            ddt1(i,j)=1./(t(i,j)-t(i+1,j))
            ddt3(i,j)=1./(t(i,j+1)-t(i+1,j+1))
         endif
      enddo
   enddo
!
! Coefficients for the derivative d/dphi (multiplied by 2)
!
   do i=1,(nx-1)
      do j=1,(ny-1)
         if(mask(i,j).ge.1.and.mask(i+1,j).ge.1.and.mask(i,j+1).ge.1.and.mask(i+1,j+1).ge.1) then
            ddp(i,j,1)=(t(i,j)+t(i,j+1)-3.*t(i+1,j)+t(i+1,j+1))/(2.*(p(j)-p(j+1))*(t(i,j)-t(i+1,j)))
            ddp(i,j,2)=-(-3.*t(i,j)+t(i,j+1)+t(i+1,j)+t(i+1,j+1))/(2.*(p(j)-p(j+1))*(t(i,j)-t(i+1,j)))
            ddp(i,j,3)=-(t(i,j)+t(i,j+1)+t(i+1,j)-3.*t(i+1,j+1))/(2.*(p(j)-p(j+1))*(t(i,j+1)-t(i+1,j+1)))
            ddp(i,j,4)=(t(i,j)-3.*t(i,j+1)+t(i+1,j)+t(i+1,j+1))/(2.*(p(j)-p(j+1))*(t(i,j+1)-t(i+1,j+1)))
         endif
      enddo
   enddo
!
! Coefficients for the derivative d/dr (multiplied by 2)
!
   do i=1,(nx-1)
      do j=1,(ny-1)
         if(mask(i,j).ge.1.and.mask(i+1,j).ge.1.and.mask(i,j+1).ge.1.and.mask(i+1,j+1).ge.1) then
            ddr(i,j,1)=(t(i,j)+t(i,j+1)-3.*t(i+1,j)+t(i+1,j+1))/(4.*t(i,j)-4.*t(i+1,j))
            ddr(i,j,2)=-(-3.*t(i,j)+t(i,j+1)+t(i+1,j)+t(i+1,j+1))/(4.*t(i,j)-4.*t(i+1,j))
            ddr(i,j,3)=(t(i,j)+t(i,j+1)+t(i+1,j)-3.*t(i+1,j+1))/(4.*t(i,j+1)-4.*t(i+1,j+1))
            ddr(i,j,4)=-((t(i,j)-3.*t (i,j+1)+t(i+1,j)+t(i+1,j+1))/(4.*t(i,j+1)-4.*t(i+1,j+1)))
         endif
      enddo
   enddo
!
! For convenience and speed replace ddp with sin(phi_0)*ddp+cos(phi_0)*ddr
! and ddr with 2*sin(phi_0)*ddr
!
   do i=1,(nx-1)
      do j=1,(ny-1)
         do k=1,4
            ddp(i,j,k)=ddp(i,j,k)*sin(0.5*(p(j)+p(j+1)))+ddr(i,j,k)*cos(0.5*(p(j)+p(j+1)))
            ddr(i,j,k)=2.*ddr(i,j,k)*sin(0.5*(p(j)+p(j+1)))
         enddo
      enddo
   enddo
!
end subroutine setup_spherical_PF_4p
!***********************************************************************************************************************************
