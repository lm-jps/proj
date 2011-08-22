!***********************************************************************************************************************************
subroutine setup_OCBP_PF_dzh_4p
!===================================================================================================================================
! This subroutine computes the coefficients used to approximate heliographic derivatives given the magnetic field and the vertical
! derivative of the vertical component of the field, dBzdz, at discrete locations. All heliographic derivatives are approximated at
! a point which is located between the centres of four neighbouring pixels: 1=(i,j), 2=(i+1,j), 3=(i,j+1), and 4=(i+1,j+1). The 
! location of this point varies with theta, phi, p, b in order to ensure that the approximations for d/dx and d/dy are second order
! accurate. The location of this point affects the coefficients for d/dz because it must be interpolated, but does not affect the
! coefficients for d/dx and d/dy.
!
! 2009, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use bounds
   use deriv_coefficients
   use mgram_data
   use pix_size
   use trnsfrm

   implicit none

   integer :: i
!
! Coefficients for the heliographic derivative d/dx
!
   ddxb(1)=0.5*((a23*a32 - a33*c22)/dxi + (a33*c12 - a23*a31)/dyi)
   ddxb(2)=0.5*((a33*c22 - a23*a32)/dxi + (a33*c12 - a23*a31)/dyi)
   ddxb(3)=-ddxb(2)
   ddxb(4)=-ddxb(1)
!
! Coefficients for the heliographic derivative d/dy
!
   ddyb(1)=0.5*((a33*c21 - a13*a32)/dxi + (a13*a31 - a33*c11)/dyi)
   ddyb(2)=0.5*((a13*a32 - a33*c21)/dxi + (a13*a31 - a33*c11)/dyi)
   ddyb(3)=-ddyb(2)
   ddyb(4)=-ddyb(1)
!
! Coefficients for the heliographic derivative d/dz
!
   ddzb(4)=4*(a33*(-a23*a31*c21 - a23*a32*c11 + a33*c12*c21 + a33*c11*c22) + &
      a13*(2*a23*a31*a32 - a33*(a32*c12 + a31*c22)))*dxi*dyi
   ddzb(1)=((a13*a31*dxi - a33*c11*dxi - a13*a32*dyi + a33*c21*dyi)*&
      (-a23*a31*dxi + a33*c12*dxi + a23*a32*dyi - a33*c22*dyi))/ddzb(4)
   ddzb(2)=((a33*(c11*dxi + c21*dyi) - a13*(a31*dxi + a32*dyi))*&
      (a33*(c12*dxi + c22*dyi) - a23*(a31*dxi + a32*dyi)))/ddzb(4)
   ddzb(3)=ddzb(2)
   ddzb(4)=ddzb(1)
!
! The following quantities are used repeatedly in CalcDE 
!
   do i=1,4
      ddxc(i)=2.*ddxb(i)
      ddyc(i)=2.*ddyb(i)
   enddo

end subroutine setup_OCBP_PF_dzh_4p
!***********************************************************************************************************************************
