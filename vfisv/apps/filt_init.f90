SUBROUTINE FILT_INIT (NUM_LAMBDA_FILTER, NUM_TUNNING, WSPACING, NUM_LAMBDA_LONG)
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  ! RCE Apr 21, 2010: Commented out almost everything (all that's 
  ! related to computing the actual filter profiles because Sebastien 
  ! is computing them in the C wrapper.
  !
  ! RCE May 5, 2010: Uncommented the line that computes TUNEPOS calculation:        
  ! TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
  ! because WFA_GUESS.F90 needs it to estimate the velocity.


  USE FILT_PARAM
  USE LINE_PARAM
  USE CONS_PARAM
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)           :: NUM_LAMBDA_FILTER, NUM_TUNNING, NUM_LAMBDA_LONG
  REAL(DP), INTENT(IN)           :: WSPACING
  INTEGER                        :: I
  !
  NUMW_LONG = NUM_LAMBDA_LONG
  NBINS = NUM_LAMBDA_FILTER
  NTUNE = NUM_TUNNING

   !-------------------------------------------------------------
   ! Allocating arrays
   !-------------------------------------------------------------
!   ALLOCATE (FILTER(NUMW,NBINS),TUNEPOS(NBINS))
   !-------------------------------------------------------------

   ! TUNEPOS is used in wfa_guess to calculate velocity guess
   DO I=1,NBINS
      TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
   ENDDO

END SUBROUTINE FILT_INIT
!CVSVERSIONINFO "$Id: filt_init.f90,v 1.5 2012/04/09 22:20:44 keiji Exp $"
