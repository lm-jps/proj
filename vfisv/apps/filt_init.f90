SUBROUTINE FILT_INIT (NUM_LAMBDA_FILTER, NUM_TUNNING, CONTINUUM, LYOTFWHM, WNARROW, WSPACING, NUM_LAMBDA_LONG)
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  ! RCE Apr 21, 2010: Commented out almost everything (all that's related to computing the 
  ! actual filter profiles because Sebastien is computing them in the C wrapper.
  !
  ! RCE May 5, 2010: Uncommented the line that computes TUNEPOS calculation:        
  ! TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
  ! because WFA_GUESS.F90 needs it to estimate the velocity.


  USE FILT_PARAM
  USE LINE_PARAM
  USE CONS_PARAM
  IMPLICIT NONE
  !
  INTEGER,  PARAMETER                        :: NMICH = 3
  INTEGER,  INTENT(IN)                       :: NUM_LAMBDA_FILTER, NUM_TUNNING, NUM_LAMBDA_LONG
  INTEGER,  INTENT(IN)                       :: CONTINUUM
  REAL(DP), INTENT(IN)                       :: LYOTFWHM, WNARROW, WSPACING
  REAL(DP), DIMENSION(NMICH)                 :: WMICH
  REAL(DP)                                   :: LYOTW
  INTEGER                                    :: I, J
  !
  !-------------------------------------------------------
  ! Check for optional input parameters
  !-------------------------------------------------------
  ! 5L:   NBINS=5; NTUNE=5, CONT=.FALSE.
  ! 6L:   NBINS=6, NTUNE=6, CONT=.FALSE.
  ! 5L+C: NBINS=6, NTUNE=11, CONT=.TRUE.
  !--------------------------------------------------------
  NUMW_LONG = NUM_LAMBDA_LONG
  NBINS = NUM_LAMBDA_FILTER
  NTUNE = NUM_TUNNING
  CONT=.FALSE.
  IF (CONTINUUM.EQ.0) CONT=.FALSE.
  IF (CONTINUUM.EQ.1) CONT=.TRUE.

  ! Consistency check
!  CHECK = 0
!  IF (NBINS.EQ.5.AND.NTUNE.EQ.5.AND.CONT.EQ..FALSE.) CHECK=CHECK+1
!  IF (NBINS.EQ.6.AND.NTUNE.EQ.6.AND.CONT.EQ..FALSE.) CHECK=CHECK+1
!  IF (NBINS.EQ.6.AND.NTUNE.EQ.11.AND.CONT.EQ..TRUE.) CHECK=CHECK+1
!  IF (CHECK.EQ.0) THEN
!     PRINT*,'Illegal combination of Filter parameters: NBINS, NTUNE, CONT'
!     PRINT*,'NBIBS:',NBINS
!     PRINT*,'NTUNE:',NTUNE
!     PRINT*,'CONT:',CONT
!     STOP
!  ENDIF

     !-------------------------------------------------------------
     ! Allocating arrays
     !-------------------------------------------------------------
     ALLOCATE (FILTER(NUMW,NBINS),TUNEPOS(NBINS))
     !-------------------------------------------------------------

     DO I=1,NBINS
        TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
     ENDDO

END SUBROUTINE FILT_INIT
