!***********************************************************************************************************************************
subroutine reconfig(i,j)
!===================================================================================================================================
! This subroutine implements the reconfiguration caused by changing the azimuthal angle at pixel (i,j). This involves:
! (1) reversing the sign of the transverse components of the field at pixel (i,j); and 
! (2) updating div(B) and Jz for each finite differencing stencil that is affected by the reconfiguration.
!
! 2009,2010, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use energy_arrays
   use mgram_data
   use recon

   implicit none

   integer :: i,j,ir
!
! Reverse the sign of the transvserse components of the field
!
   Bx(i,j)=-Bx(i,j)
   By(i,j)=-By(i,j)
!
! Update div(B) and Jz
!
   do ir=1,4
      DivB(i_reconfig(ir),j_reconfig(ir))=DivB_reconfig(ir)
      Jz(i_reconfig(ir),j_reconfig(ir))=Jz_reconfig(ir)
   enddo

end subroutine reconfig
!***********************************************************************************************************************************
