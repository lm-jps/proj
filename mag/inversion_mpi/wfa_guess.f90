PURE SUBROUTINE WFA_GUESS(OBS,LAMBDA_0,CAMPO,INCLI,AZIM,VELOCITY, TOTPOL)
  ! J M Borrero
  ! Oct 28, 2007
  ! HAO-NCAR for HMI-Stanford
  USE CONS_PARAM
  USE INV_PARAM
  USE FILT_PARAM
  IMPLICIT NONE
 
  REAL(8), INTENT(IN)                       :: OBS(NBINS,4)
  REAL(8), INTENT(IN)                       :: LAMBDA_0
 ! REAL(DP), INTENT(IN), DIMENSION(NTUNE)   :: TUNEPOS
  !LOGICAL, INTENT(IN)                      :: QUICKLOOK
  REAL(8), INTENT(OUT)                      :: CAMPO, INCLI, AZIM, VELOCITY, TOTPOL
  INTEGER                                   :: L, K, J, W
  REAL(DP)                                  :: Q_INT, U_INT, V_INT, VDIFF, VSIGN, BLOS, BT
  REAL(DP), DIMENSION(NBINS)                :: POP, DPOP, DDPOP, AVPROF
  !
!  Q_INT=SUM(ABS(OBS(:,2)/MAXVAL(OBS(:,1))))/NBINS
!  U_INT=SUM(ABS(OBS(:,3)/MAXVAL(OBS(:,1))))/NBINS
!  V_INT=SUM(ABS(OBS(:,4)/MAXVAL(OBS(:,1))))/NBINS
  Q_INT=SUM(ABS(OBS(:,2)))/NBINS
  U_INT=SUM(ABS(OBS(:,3)))/NBINS
  V_INT=SUM(ABS(OBS(:,4)))/NBINS
  TOTPOL=Q_INT**2D0+U_INT**2D0+V_INT**2D0
!PRINT*, "STOKES V = ", OBS(:,4)
!PRINT*, "V_INT = ", ABS(OBS(:,4)), ABS(OBS(:,1)), SUM(ABS(OBS(:,4)/MAXVAL(OBS(:,1)))),V_INT
!PRINT*, "OBS = ", OBS(1,1)
  !--------------------------------------------
  ! Velocity: we use the COG method in Stokes I
  !--------------------------------------------
  AVPROF(:)=0D0
  DO L=1,NTUNE
     AVPROF(L)=MAXVAL(OBS(:,1))-OBS(L,1)
  ENDDO
  VELOCITY=SUM(AVPROF*TUNEPOS)/SUM(AVPROF)
  VELOCITY=1E-3*(2.998E10/LAMBDA_0)*VELOCITY    ! cm/s
  !----------------------------------------------------------
  ! Field Strength using WFA: constants obtained by S.Tomczyk
  !----------------------------------------------------------
  BLOS=V_INT
  BLOS=-3.8924328D0+7112.1008D0*BLOS
  BT=(Q_INT**2D0+U_INT**2D0)**0.25D0
  BT=-20.492311D0+5134.6516D0*BT
  CAMPO=SQRT(BT**2D0+BLOS**2D0)
!PRINT*, "CAMPO = ", CAMPO
  !-------------------------------------------------------------
  ! Field Inclination using WFA: constants obtained by S.Tomczyk
  !-------------------------------------------------------------
  IF (ABS(VELOCITY).LT.2E5) VDIFF=SUM(OBS(1:3,4))-SUM(OBS(4:6,4))
  IF (VELOCITY.GT.2E5) VDIFF=SUM(OBS(2:4,4))-SUM(OBS(5:6,4))
  IF (VELOCITY.LT.-2E5) VDIFF=SUM(OBS(1:2,4))-SUM(OBS(3:4,4))
  VSIGN=VDIFF/ABS(VDIFF)
  BLOS=V_INT
  BLOS=-3.8924328D0+7112.1008D0*BLOS
  BT=(Q_INT**2D0+U_INT**2D0)**0.25D0
  BT=-20.492311D0+5134.6516D0*BT
  INCLI=-VSIGN*ATAN2(BLOS,BT)/D2R+90D0
  !-------------------------------------------------
  ! Regardless of the polarization level
  ! Field Azimuth: Auer's SoPh 1977, 55, 47
  !-------------------------------------------------
  AZIM=0.25D0*ATAN2(2D0*SUM(OBS(:,2)*OBS(:,3)),SUM(OBS(:,2)**2D0-OBS(:,3)**2D0))/D2R
  POP=OBS(:,2)*COS(2D0*AZIM*D2R)+OBS(:,3)*SIN(2D0*AZIM*D2R)
  DO W=1,NBINS-1
     DPOP(W)=POP(W+1)-POP(W)
  ENDDO
  DO W=2,NBINS-2
     DDPOP(W)=DPOP(W+1)-DPOP(W)
  ENDDO
  IF (DDPOP(4).GT.0.AND.INCLI.LT.90) AZIM=AZIM+90D0
  IF (DDPOP(4).GT.0.AND.INCLI.GT.90) AZIM=AZIM+90D0
  !
  IF (AZIM.GT.180D0) AZIM=AZIM-180D0
  IF (AZIM.LT.0) AZIM=AZIM+180D0
  !
  ! By RCE: commented out next sentence because QUICKLOOK is not defined in the new 
  ! version of the code
  !IF (QUICKLOOK.EQ..FALSE.) TOTPOL=1D0
  !

END SUBROUTINE WFA_GUESS
