!***********************************************************************************************************************************
subroutine CalcE_OCBP_PF_dzh_4p(E)
!===================================================================================================================================
! This subroutine calculates the divergence of the field div(B) and the associated vertical current density Jz. Derivatives in the
! heliographic directions, xh and yh, are approximated with a 4-point finite differencing stencil. Note that div(B) and Jz are
! evaluated at a point which is located between the centres of four neighbouring pixels: 1=(i,j), 2=(i+1,j), 3=(i,j+1), and 
! 4=(i+1,j+1) where the approximations for d/dx and d/dy are second order accurate. The derivative dBzdz is also approximated at
! this point (by interpolation) as required.
! This subroutine also calculates the "energy", which is |div(B)|+lambda*|Jz| summed over the region of interest.
!
! 2009, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use bounds
   use deriv_coefficients
   use energy_arrays
   use mgram_data
   use trnsfrm
   use WeightingFactor

   implicit none

   integer :: i,j
   real :: E

!
! Allocate memory for energy arrays
!
   if(.not.allocated(DivB)) allocate(DivB(nx,ny))
   if(.not.allocated(Jz)) allocate(Jz(nx,ny))

   E=0.
   do i=nxa,(nxb-1)
      do j=nya,(nyb-1)
         DivB(i,j)=ddxb(1)*(c11*Bx(i,j)+c21*By(i,j)+a13*Bz(i,j))&
                  +ddxb(2)*(c11*Bx(i+1,j)+c21*By(i+1,j)+a13*Bz(i+1,j))&
                  +ddxb(3)*(c11*Bx(i,j+1)+c21*By(i,j+1)+a13*Bz(i,j+1))&
                  +ddxb(4)*(c11*Bx(i+1,j+1)+c21*By(i+1,j+1)+a13*Bz(i+1,j+1))&
                  +ddyb(1)*(c12*Bx(i,j)+c22*By(i,j)+a23*Bz(i,j))&
                  +ddyb(2)*(c12*Bx(i+1,j)+c22*By(i+1,j)+a23*Bz(i+1,j))&
                  +ddyb(3)*(c12*Bx(i,j+1)+c22*By(i,j+1)+a23*Bz(i,j+1))&
                  +ddyb(4)*(c12*Bx(i+1,j+1)+c22*By(i+1,j+1)+a23*Bz(i+1,j+1))&
                  +ddzb(1)*dBzdz(i,j)+ddzb(2)*dBzdz(i+1,j)+ddzb(3)*dBzdz(i,j+1)+ddzb(4)*dBzdz(i+1,j+1)

         Jz(i,j)=ddxb(1)*(c12*Bx(i,j)+c22*By(i,j)+a23*Bz(i,j))&
                +ddxb(2)*(c12*Bx(i+1,j)+c22*By(i+1,j)+a23*Bz(i+1,j))&
                +ddxb(3)*(c12*Bx(i,j+1)+c22*By(i,j+1)+a23*Bz(i,j+1))&
                +ddxb(4)*(c12*Bx(i+1,j+1)+c22*By(i+1,j+1)+a23*Bz(i+1,j+1))&
                -ddyb(1)*(c11*Bx(i,j)+c21*By(i,j)+a13*Bz(i,j))&
                -ddyb(2)*(c11*Bx(i+1,j)+c21*By(i+1,j)+a13*Bz(i+1,j))&
                -ddyb(3)*(c11*Bx(i,j+1)+c21*By(i,j+1)+a13*Bz(i,j+1))&
                -ddyb(4)*(c11*Bx(i+1,j+1)+c21*By(i+1,j+1)+a13*Bz(i+1,j+1))

         E=E+abs(DivB(i,j))+lambda*abs(Jz(i,j))
      enddo
   enddo

end subroutine CalcE_OCBP_PF_dzh_4p
!***********************************************************************************************************************************

