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
  INTEGER                                    :: I
  !
  FREE(:) = .FALSE.
  NUMFREE_DEG = 4*NBINS-10
  !
  DO I=1,10
     IF (LIST_FREE_PARAMS(I).EQ.1) THEN
        FREE(I) = .TRUE.
        NUMFREE_DEG = NUMFREE_DEG + 1
     ENDIF
  ENDDO
  WRITE(*,*) 'freed parameters'
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
!CVSVERSIONINFO "$Id: free_init.f90,v 1.2 2011/05/31 22:24:18 keiji Exp $"
