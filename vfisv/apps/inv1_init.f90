SUBROUTINE INV1_INIT (ICONT,REGUL_FLAG)
! Initialize limits and norms

  USE INV_PARAM
  IMPLICIT NONE
  REAL(DP), INTENT(IN)           :: ICONT
  INTEGER, INTENT(IN)            :: REGUL_FLAG

! Used to be in inv_utils.f90
! Limits for the model parameters. Used in FINE_TUNE_MODEL
! Order: eta0, inclination (deg), azimuth (deg), damping, Doppler width, 
! field strength (gauss), line-of-sight velocity (cm/s), source function, 
! source function gradient, filling factor

  ICONT_REF=ICONT

  LOWER_LIMIT = (/1D0, 0D0 , 0D0 , 1D-4,1D0,5D0,-7D5, 1.5E-1*ICONT, 1.5E-1*ICONT, 0D0/)
  UPPER_LIMIT  = (/1D3, 180D0,180D0,5D0, 5D2,5D3, 7D5, 1.2D0*ICONT,  1.2D0*ICONT,  1D0/)

  !NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.5D0/)
  NORM(:)=(/0.2D0,7D0,50D0,1D0,0.4D0,50D0,5000D0,0.003D0*ICONT,0.004D0*ICONT,0.5D0/)

! Limits on DMODEL
! I have tried with more generous limits at this is the best combination I could find
  DLIMIT(:)=(/10D0,25D0,25D0,0.1D0,30D0,500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.25D0/)
! Relative limits. 0 means don't use
  RLIMIT(:)=(/0.0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/)
 
END SUBROUTINE INV1_INIT
!CVSVERSIONINFO "$Id: inv1_init.f90,v 1.1 2012/04/09 22:21:10 keiji Exp $"
