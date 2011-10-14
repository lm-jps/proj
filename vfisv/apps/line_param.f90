MODULE LINE_PARAM
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE CONS_PARAM
  !
  REAL(DP)                                   :: LANDA0, SHIFT, STEPW, NOISE
  INTEGER                                    :: NUMW
  REAL(DP),         ALLOCATABLE              :: WAVE(:)
  !
END MODULE LINE_PARAM
!CVSVERSIONINFO "$Id: line_param.f90,v 1.3 2011/10/14 17:23:06 keiji Exp $"
