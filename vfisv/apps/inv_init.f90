SUBROUTINE INV_INIT (NUM_ITERATIONS, SVD_TOLERANCE, CHI2_STOP, POLARIZATION_THRESHOLD, &
     INTENSITY_THRESHOLD, PERCENTAGE_JUMP)
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE CONS_PARAM
  USE INV_PARAM
  IMPLICIT NONE
  REAL(DP), INTENT(IN), OPTIONAL :: SVD_TOLERANCE, CHI2_STOP, POLARIZATION_THRESHOLD
  REAL(DP), INTENT(IN), OPTIONAL :: INTENSITY_THRESHOLD, PERCENTAGE_JUMP
  INTEGER,  INTENT(IN)           :: NUM_ITERATIONS
  !
  ITER=NUM_ITERATIONS
  SVDTOL = SVD_TOLERANCE
  GOODFIT = CHI2_STOP
  TREPOL =  POLARIZATION_THRESHOLD
  TREIC = INTENSITY_THRESHOLD
  RANDOM_JUMP = PERCENTAGE_JUMP
  !
END SUBROUTINE INV_INIT
!CVSVERSIONINFO "$Id: inv_init.f90,v 1.4 2012/04/09 22:21:23 keiji Exp $"
