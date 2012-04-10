SUBROUTINE FREE_INIT (LIST_FREE_PARAMS)
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE INV_PARAM
  USE FILT_PARAM
  IMPLICIT NONE
  INTEGER, INTENT(IN), DIMENSION(10)         :: LIST_FREE_PARAMS
  INTEGER                                    :: I,J
  !
  FREE(:) = .FALSE.
  NUMFREE_PARAM = 0
  NUMFREE_DEG = 4*NBINS
  !
  DO I=1,10
     IF (LIST_FREE_PARAMS(I).EQ.1) THEN
        FREE(I) = .TRUE.
        NUMFREE_PARAM=NUMFREE_PARAM+1
        NUMFREE_DEG = NUMFREE_DEG - 1
     ENDIF
  ENDDO

  ALLOCATE (FREELOC(NUMFREE_PARAM))
  J=1
  DO I=1,10
     IF (FREE(I).EQ..TRUE.) THEN
        FREELOC(J)=I
        J=J+1
     ENDIF
  ENDDO

  !FREE(1)                                    ! ETA0
  !FREE(2)                                    ! FIELD INCLINATION
  !FREE(3)                                    ! FIELD AZIMUTH
  !FREE(4)                                    ! DAMPING
  !FREE(5)                                    ! DOPPLER WIDTH
  !FREE(6)                                    ! FIELD STREGNTH
  !FREE(7)                                    ! MAGNETIC LOS VELOCITY
  !FREE(8)                                    ! SOURCE FUNCTION CONTINUUM
  !FREE(9)                                    ! SOURCE FUNCTION GRADIENT
  !FREE(10)                                   ! MAGNETIC FILLING FACTOR
END SUBROUTINE FREE_INIT
!CVSVERSIONINFO "$Id: free_init.f90,v 1.5 2012/04/10 22:16:39 keiji Exp $"
