MODULE INV_PARAM
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  ! By RCE, April 2012: Most of the variables defined here are 
  ! initialized in calls to INV_INIT and LIM_INIT from
  ! the wrapper. 

  USE CONS_PARAM
  IMPLICIT NONE
  LOGICAL,  DIMENSION(10)            :: FREE
  INTEGER,     ALLOCATABLE           :: FREELOC(:)
  INTEGER                            :: ITER, NUMFREE_PARAM, NUMFREE_DEG
  REAL(DP)                           :: SVDTOL, GOODFIT, TREPOL, TREIC, ICONT
  REAL(DP), DIMENSION(10)            :: LOWER_LIMIT,UPPER_LIMIT
  REAL(DP)                           :: ICONT_REF
  REAL(DP), DIMENSION(10)            :: DLIMIT,RLIMIT
  REAL(DP), DIMENSION(10)            :: NORM
  REAL(DP)                           :: DELTACHIMIN, LAMBDA_UP, LAMBDA_DOWN, LAMBDA_MIN=1D-4, LAMBDA_MAX

END MODULE INV_PARAM
!CVSVERSIONINFO "$Id: inv_param.f90,v 1.5 2012/04/10 22:16:59 keiji Exp $"
