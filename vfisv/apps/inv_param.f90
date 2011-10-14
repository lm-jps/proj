MODULE INV_PARAM
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE CONS_PARAM
  IMPLICIT NONE
  LOGICAL,  DIMENSION(10)            :: FREE
  INTEGER                            :: ITER, NUMFREE_PARAM, NUMFREE_DEG
  REAL(DP)                           :: SVDTOL
  REAL(DP)                           :: GOODFIT
  REAL(DP)                           :: TREPOL
  REAL(DP)                           :: TREIC
  REAL(DP)                           :: RANDOM_JUMP
  !
END MODULE INV_PARAM
!CVSVERSIONINFO "$Id: inv_param.f90,v 1.3 2011/10/14 17:22:51 keiji Exp $"
