MODULE SVD_PARAM
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  USE CONS_PARAM
  IMPLICIT NONE
  REAL(DP),    ALLOCATABLE  :: HESS(:,:), DIVC(:)
  INTEGER,     ALLOCATABLE  :: FREELOC(:)

!  WRITE(*,*) 'done svd_param'
END MODULE SVD_PARAM
!CVSVERSIONINFO "$Id: svd_param.f90,v 1.2 2011/05/31 22:25:25 keiji Exp $"
