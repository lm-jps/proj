PURE SUBROUTINE WFA_GUESS(OBS,LAMBDA_0,VELOCITY, AZIMUTH, ETA0, DOPPLERW, S0, S1)
  ! J M Borrero
  ! Oct 28, 2007
  ! HAO-NCAR for HMI-Stanford
  USE CONS_PARAM
  USE INV_PARAM
  USE FILT_PARAM
  USE LINE_PARAM
  IMPLICIT NONE
 
  REAL(8), INTENT(IN)                       :: OBS(NBINS,4)
  REAL(8), INTENT(IN)                       :: LAMBDA_0
  REAL(8), INTENT(OUT)                      :: VELOCITY, AZIMUTH, ETA0, DOPPLERW, S0, S1
  INTEGER                                   :: L
  REAL(DP), DIMENSION(NBINS)                :: AVPROF

! By RCE: April 2012: moved to this routine all initializations that override the wrapper's initial guess.

!-- ETA0
     ETA0 = 5
!-- DOPPLER WIDTH --
     DOPPLERW = 20
!-- AZIMUTH
     AZIMUTH = atan2(sum(obs(:,3)),sum(obs(:,2)))*90/DPI
! The following is slightly better
    S0=0.15D0*0.98D0*ICONT
    S1=0.85D0*0.98D0*ICONT

  !--------------------------------------------
  ! Velocity: we use the COG method in Stokes I
  !--------------------------------------------
  AVPROF(:)=0D0
  DO L=1,NBINS
     AVPROF(L)=MAXVAL(OBS(:,1))-OBS(L,1)
  ENDDO

  VELOCITY=SUM(AVPROF*TUNEPOS)/SUM(AVPROF)
  IF (MAXVAL(ABS(OBS(:,4))) .GT. 2D0*NOISE(1)) &
       VELOCITY = SUM(ABS(OBS(:,4))*TUNEPOS)/SUM(ABS(OBS(:,4)))
  VELOCITY=1E-3*(2.998E10/LAMBDA_0)*VELOCITY    ! cm/s

END SUBROUTINE WFA_GUESS
!CVSVERSIONINFO "$Id: wfa_guess.f90,v 1.5 2012/04/10 22:18:01 keiji Exp $"
