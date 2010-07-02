!***********************************************************************************************************************************
subroutine CalcDE_reconfig_spherical_PF_4p(i,j)
!===================================================================================================================================
! For a reconfiguration of the azimuthal angle at pixel (i,j) this subroutine computes Div(B) and Jz at each point with a finite 
! differencing stencil that references pixel (i,j).
!
! 2009,2010, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use bounds
   use energy_arrays
   use mgram_data
   use recon
   use spherical_deriv_coefficients
   use spherical_position

   implicit none

   integer :: i,j
   real :: Br_ambig,Bp_ambig,Bx_cost

   Br_ambig=sinp(j)*sint(i,j)*Bx(i,j)+cosp(j)*By(i,j)
   Bp_ambig=cosp(j)*sint(i,j)*Bx(i,j)-sinp(j)*By(i,j)
   Bx_cost=cost(i,j)*Bx(i,j)

   i_reconfig(1)=i
   j_reconfig(1)=j
   DivB_reconfig(1)=DivB(i,j)-ddr(i,j,1)*Br_ambig-ddp(i,j,1)*Bp_ambig-ddt1(i,j)*Bx_cost
   Jz_reconfig(1)=Jz(i,j)-Bx_cost*ddp(i,j,1)+ddt1(i,j)*Bp_ambig

   i_reconfig(2)=i
   j_reconfig(2)=j-1
   DivB_reconfig(2)=DivB(i,j-1)-ddr(i,j-1,3)*Br_ambig-ddp(i,j-1,3)*Bp_ambig-ddt3(i,j-1)*Bx_cost
   Jz_reconfig(2)=Jz(i,j-1)-Bx_cost*ddp(i,j-1,3)+ddt3(i,j-1)*Bp_ambig

   i_reconfig(3)=i-1
   j_reconfig(3)=j
   DivB_reconfig(3)=DivB(i-1,j)-ddr(i-1,j,2)*Br_ambig-ddp(i-1,j,2)*Bp_ambig+ddt1(i-1,j)*Bx_cost
   Jz_reconfig(3)=Jz(i-1,j)-Bx_cost*ddp(i-1,j,2)-ddt1(i-1,j)*Bp_ambig

   i_reconfig(4)=i-1
   j_reconfig(4)=j-1
   DivB_reconfig(4)=DivB(i-1,j-1)-ddr(i-1,j-1,4)*Br_ambig-ddp(i-1,j-1,4)*Bp_ambig+ddt3(i-1,j-1)*Bx_cost
   Jz_reconfig(4)=Jz(i-1,j-1)-Bx_cost*ddp(i-1,j-1,4)-ddt3(i-1,j-1)*Bp_ambig

end subroutine CalcDE_reconfig_spherical_PF_4p
!***********************************************************************************************************************************
