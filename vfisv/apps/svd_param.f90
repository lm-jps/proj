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
!CVSVERSIONINFO "$Id: svd_param.f90,v 1.3 2011/10/14 17:23:32 keiji Exp $"
