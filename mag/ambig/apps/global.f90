subroutine global(bthresh,setup,CalcE,CalcDE_reconfig)
!===================================================================================================================================
! This subroutine calls the subroutines for the global minimisation.
!
! 2009, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   implicit none

   real :: bthresh
   external :: setup,CalcE,CalcDE_reconfig

   call setup
   call minimise_energy(CalcE,CalcDE_reconfig)

end subroutine global
