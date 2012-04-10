SUBROUTINE FILT_INIT (NUM_LAMBDA_FILTER,  WSPACING, NUM_LAMBDA_LONG)
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  ! RCE Apr 21, 2010: Commented out almost everything (all that's 
  ! related to computing the actual filter profiles because they are
  ! computed in the C wrapper (by Sebastien Couvidat).
  !
  ! RCE May 5, 2010: Uncommented the line that computes TUNEPOS calculation:        
  ! TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
  ! because WFA_GUESS.F90 needs it to estimate the velocity.


  USE FILT_PARAM
  USE LINE_PARAM
  USE CONS_PARAM
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)           :: NUM_LAMBDA_FILTER, NUM_LAMBDA_LONG
  REAL(DP), INTENT(IN)           :: WSPACING
  INTEGER                        :: I
  !
  NUMW_LONG = NUM_LAMBDA_LONG
  NBINS = NUM_LAMBDA_FILTER

   ! TUNEPOS is used in wfa_guess to calculate velocity guess
   DO I=1,NBINS
      TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
   ENDDO

END SUBROUTINE FILT_INIT
!CVSVERSIONINFO "$Id: filt_init.f90,v 1.6 2012/04/10 22:16:24 keiji Exp $"
