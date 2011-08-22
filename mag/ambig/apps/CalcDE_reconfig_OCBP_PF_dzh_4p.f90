!***********************************************************************************************************************************
subroutine CalcDE_reconfig_OCBP_PF_dzh_4p(i,j)
!===================================================================================================================================
! For a reconfiguration of the azimuthal angle at pixel (i,j) this subroutine computes Div(B) and Jz at each point with a finite 
! differencing stencil that references pixel (i,j). For most points this involves four neighbouring pixels. At the corners and edges 
! it involves less than four pixels.
!
! 2009,2010, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use bounds
   use deriv_coefficients
   use energy_arrays
   use mgram_data
   use recon
   use trnsfrm

   implicit none

   integer :: i,j
   real :: Bx_ambig,By_ambig

   Bx_ambig=c11*Bx(i,j)+c21*By(i,j)
   By_ambig=c12*Bx(i,j)+c22*By(i,j)

   i_reconfig(1)=i
   j_reconfig(1)=j
   DivB_reconfig(1)=DivB(i,j)-ddxc(1)*Bx_ambig-ddyc(1)*By_ambig
   Jz_reconfig(1)=Jz(i,j)-ddxc(1)*By_ambig+ddyc(1)*Bx_ambig

   i_reconfig(2)=i
   j_reconfig(2)=j-1
   DivB_reconfig(2)=DivB(i,j-1)-ddxc(3)*Bx_ambig-ddyc(3)*By_ambig
   Jz_reconfig(2)=Jz(i,j-1)-ddxc(3)*By_ambig+ddyc(3)*Bx_ambig

   i_reconfig(3)=i-1
   j_reconfig(3)=j
   DivB_reconfig(3)=DivB(i-1,j)-ddxc(2)*Bx_ambig-ddyc(2)*By_ambig
   Jz_reconfig(3)=Jz(i-1,j)-ddxc(2)*By_ambig+ddyc(2)*Bx_ambig

   i_reconfig(4)=i-1
   j_reconfig(4)=j-1
   DivB_reconfig(4)=DivB(i-1,j-1)-ddxc(4)*Bx_ambig-ddyc(4)*By_ambig
   Jz_reconfig(4)=Jz(i-1,j-1)-ddxc(4)*By_ambig+ddyc(4)*Bx_ambig

end subroutine CalcDE_reconfig_OCBP_PF_dzh_4p
!***********************************************************************************************************************************
