!***********************************************************************************************************************************
subroutine CalcE_spherical_PF_4p(E)
!===================================================================================================================================
! This subroutine calculates the divergence of the field div(B) and the radial current density Jr in spherical coordinates. 
! Derivatives in the theta- and phi-directions are approximated with a 4-point finite differencing stencil. To improve the order of 
! accuracy div(B) and Jr are evaluated at a point which is located between the centres of four neighbouring pixels: 1=(i,j), 
! 2=(i+1,j), 3=(i,j+1), and 4=(i+1,j+1). The derivative r*dBrdr is approximated at this point by interpolation from the tabulated 
! values at each pixel. It is assumed that the position phi (array p) is a function of yi only (this assumption could be relaxed but
! it will change the coefficients for the derivatives).
! Note that r*dBr/dr is stored in the array dBzdz.
!
! This subroutine also calculates the "energy", which is |div(B)|+lambda*|Jr| summed over the region of interest.
!
! 2009, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use bounds
   use energy_arrays
   use mgram_data
   use spherical_deriv_coefficients
   use spherical_position
   use WeightingFactor

   implicit none

   integer :: i,j
   real :: E

   E=0.
   do i=nxa,(nxb-1)
      do j=nya,(nyb-1)
!
! r*sin(phi)*div(B)
!
         DivB(i,j)=&
    0.5*ddr(i,j,1)*(sinp(j  )*sint(i  ,j  )*Bx(i  ,j  )+cosp(j  )*By(i  ,j  )+sinp(j  )*cost(i  ,j  )*Bz(i  ,j  ))&
   +0.5*ddr(i,j,2)*(sinp(j  )*sint(i+1,j  )*Bx(i+1,j  )+cosp(j  )*By(i+1,j  )+sinp(j  )*cost(i+1,j  )*Bz(i+1,j  ))&
   +0.5*ddr(i,j,3)*(sinp(j+1)*sint(i  ,j+1)*Bx(i  ,j+1)+cosp(j+1)*By(i  ,j+1)+sinp(j+1)*cost(i  ,j+1)*Bz(i  ,j+1))&
   +0.5*ddr(i,j,4)*(sinp(j+1)*sint(i+1,j+1)*Bx(i+1,j+1)+cosp(j+1)*By(i+1,j+1)+sinp(j+1)*cost(i+1,j+1)*Bz(i+1,j+1))&
  +0.25*ddr(i,j,1)*dBzdz(i,j)+0.25*ddr(i,j,2)*dBzdz(i+1,j)+0.25*ddr(i,j,3)*dBzdz(i,j+1)+0.25*ddr(i,j,4)*dBzdz(i+1,j+1)&
   +0.5*ddp(i,j,1)*(cosp(j  )*sint(i  ,j  )*Bx(i  ,j  )-sinp(j  )*By(i  ,j  )+cosp(j  )*cost(i  ,j  )*Bz(i  ,j  ))&
   +0.5*ddp(i,j,2)*(cosp(j  )*sint(i+1,j  )*Bx(i+1,j  )-sinp(j  )*By(i+1,j  )+cosp(j  )*cost(i+1,j  )*Bz(i+1,j  ))&
   +0.5*ddp(i,j,3)*(cosp(j+1)*sint(i  ,j+1)*Bx(i  ,j+1)-sinp(j+1)*By(i  ,j+1)+cosp(j+1)*cost(i  ,j+1)*Bz(i  ,j+1))&
   +0.5*ddp(i,j,4)*(cosp(j+1)*sint(i+1,j+1)*Bx(i+1,j+1)-sinp(j+1)*By(i+1,j+1)+cosp(j+1)*cost(i+1,j+1)*Bz(i+1,j+1))&
  +0.5*ddt1(i,j)*(cost(i  ,j  )*Bx(i  ,j  )-sint(i  ,j  )*Bz(i  ,j  ))&
  -0.5*ddt1(i,j)*(cost(i+1,j  )*Bx(i+1,j  )-sint(i+1,j  )*Bz(i+1,j  ))&
  +0.5*ddt3(i,j)*(cost(i  ,j+1)*Bx(i  ,j+1)-sint(i  ,j+1)*Bz(i  ,j+1))&
  -0.5*ddt3(i,j)*(cost(i+1,j+1)*Bx(i+1,j+1)-sint(i+1,j+1)*Bz(i+1,j+1))
!
! r*sin(phi)*Jr
!
         Jz(i,j)=&
    0.5*ddp(i,j,1)*(cost(i  ,j  )*Bx(i  ,j  )-sint(i  ,j  )*Bz(i  ,j  ))&
   +0.5*ddp(i,j,2)*(cost(i+1,j  )*Bx(i+1,j  )-sint(i+1,j  )*Bz(i+1,j  ))&
   +0.5*ddp(i,j,3)*(cost(i  ,j+1)*Bx(i  ,j+1)-sint(i  ,j+1)*Bz(i  ,j+1))&
   +0.5*ddp(i,j,4)*(cost(i+1,j+1)*Bx(i+1,j+1)-sint(i+1,j+1)*Bz(i+1,j+1))&
  -0.5*ddt1(i,j)*(cosp(j  )*sint(i  ,j  )*Bx(i  ,j  )-sinp(j  )*By(i  ,j  )+cosp(j  )*cost(i  ,j  )*Bz(i  ,j  ))&
  +0.5*ddt1(i,j)*(cosp(j  )*sint(i+1,j  )*Bx(i+1,j  )-sinp(j  )*By(i+1,j  )+cosp(j  )*cost(i+1,j  )*Bz(i+1,j  ))&
  -0.5*ddt3(i,j)*(cosp(j+1)*sint(i  ,j+1)*Bx(i  ,j+1)-sinp(j+1)*By(i  ,j+1)+cosp(j+1)*cost(i  ,j+1)*Bz(i  ,j+1))&
  +0.5*ddt3(i,j)*(cosp(j+1)*sint(i+1,j+1)*Bx(i+1,j+1)-sinp(j+1)*By(i+1,j+1)+cosp(j+1)*cost(i+1,j+1)*Bz(i+1,j+1))

         E=E+abs(DivB(i,j))+lambda*abs(Jz(i,j))
      enddo
   enddo

end subroutine CalcE_spherical_PF_4p
!***********************************************************************************************************************************

