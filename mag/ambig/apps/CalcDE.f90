!***********************************************************************************************************************************
subroutine CalcDE(i,j,de,CalcDE_reconfig)
!===================================================================================================================================
! This subroutine calculates the change in the "energy" caused by reconfiguring the azimuthal angle at pixel (i,j). Given Div(B) and
! Jz at the relevant finite differencing stencils for the reconfigured state, the change in the energy is determined by examining
! the appropriate sums over the affected stencils
!
! 2009,2010, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use energy_arrays
   use recon
   use WeightingFactor

   implicit none

   integer :: i,j,ir
   real :: E,Eprev,de
   external :: CalcDE_reconfig
!
! Compute Div(B) and Jz after the reconfiguration
!
   call CalcDE_reconfig(i,j)
!
! Calculate the sum of the contributions to the energy before and after the reconfiguration
!
   Eprev=0.
   E=0.
   do ir=1,4
      Eprev=Eprev+abs(DivB(i_reconfig(ir),j_reconfig(ir)))+lambda*abs(Jz(i_reconfig(ir),j_reconfig(ir)))
      E=E+abs(DivB_reconfig(ir))+lambda*abs(Jz_reconfig(ir))
   enddo
!
! Calculate the change in energy due to the reconfiguration
!
   de=E-Eprev

end subroutine CalcDE
!***********************************************************************************************************************************
