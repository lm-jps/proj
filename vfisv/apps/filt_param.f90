MODULE FILT_PARAM
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE CONS_PARAM
  USE LINE_PARAM
  INTEGER                         :: NBINS, NTUNE, NUMW_LONG
  LOGICAL                         :: CONT
  REAL(DP), ALLOCATABLE           :: FILTER(:,:)
  REAL(DP), ALLOCATABLE           :: TUNEPOS(:)
END MODULE FILT_PARAM
!CVSVERSIONINFO "$Id: filt_param.f90,v 1.3 2011/05/31 22:24:06 keiji Exp $"
