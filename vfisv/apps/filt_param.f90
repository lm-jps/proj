MODULE FILT_PARAM
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE CONS_PARAM
  USE LINE_PARAM
  INTEGER                         :: NBINS, NUMW_LONG
  REAL(DP), ALLOCATABLE           :: FILTER(:,:)
  REAL(DP), ALLOCATABLE           :: TUNEPOS(:)
END MODULE FILT_PARAM
!CVSVERSIONINFO "$Id: filt_param.f90,v 1.6 2012/04/10 22:16:29 keiji Exp $"
