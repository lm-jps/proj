MODULE INV_PARAM
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE CONS_PARAM
  IMPLICIT NONE
  LOGICAL,  DIMENSION(10)            :: FREE
  INTEGER,     ALLOCATABLE           :: FREELOC(:)
  INTEGER                            :: ITER, NUMFREE_PARAM, NUMFREE_DEG
  REAL(DP)                           :: SVDTOL
  REAL(DP)                           :: GOODFIT
  REAL(DP)                           :: TREPOL
  REAL(DP)                           :: TREIC
  REAL(DP)                           :: RANDOM_JUMP
  REAL(DP), DIMENSION(10)            :: LOWER_LIMIT,UPPER_LIMIT
  REAL(DP)                           :: ICONT_REF
  REAL(DP), DIMENSION(10)            :: DLIMIT,RLIMIT
  REAL(DP), DIMENSION(10)            :: NORM
  REAL(DP)                           :: DELTACHIMIN, LAMBDA_UP, LAMBDA_DOWN, LAMBDA_MIN=1D-4, LAMBDA_MAX

END MODULE INV_PARAM
!CVSVERSIONINFO "$Id: inv_param.f90,v 1.4 2012/04/09 22:21:28 keiji Exp $"
