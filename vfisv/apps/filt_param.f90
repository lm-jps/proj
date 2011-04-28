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
