SUBROUTINE FILT_INIT (NUM_LAMBDA_FILTER, NUM_TUNNING, CONTINUUM, LYOTFWHM, WNARROW, WSPACING)
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  ! RCE Apr 21, 2010: Commented out almost everything (all that's related to computing the 
  ! actual filter profiles because Sebastien is computing them in the C wrapper.
  !
  ! RCE May 5, 2010: Uncommented the line that computes TUNEPOS calculation:        
  ! TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
  ! because WFA_GUESS.F90 needs it to estimate the velocity.


  USE FILT_PARAM
  USE LINE_PARAM
  USE CONS_PARAM
  IMPLICIT NONE
  !
  INTEGER,  PARAMETER                        :: NMICH = 3
  INTEGER,  INTENT(IN)                       :: NUM_LAMBDA_FILTER, NUM_TUNNING
  INTEGER,  INTENT(IN)                       :: CONTINUUM
  REAL(DP), INTENT(IN)                       :: LYOTFWHM, WNARROW, WSPACING
  REAL(DP), ALLOCATABLE                      :: T1(:), LYOT(:), TEMP_TUNE(:)
  REAL(DP), ALLOCATABLE                      :: FILTERR(:,:), TEMP1(:,:), TEMP2(:,:)
  REAL(DP), DIMENSION(NMICH)                 :: WMICH
  REAL(DP)                                   :: LYOTW
  INTEGER                                    :: I, J, CHECK
  !
  !-------------------------------------------------------
  ! Check for optional input parameters
  !-------------------------------------------------------
  ! 5L:   NBINS=5; NTUNE=5, CONT=.FALSE.
  ! 6L:   NBINS=6, NTUNE=6, CONT=.FALSE.
  ! 5L+C: NBINS=6, NTUNE=11, CONT=.TRUE.
  !--------------------------------------------------------
  NBINS = NUM_LAMBDA_FILTER
  NTUNE = NUM_TUNNING
  CONT=.FALSE.
  IF (CONTINUUM.EQ.0) CONT=.FALSE.
  IF (CONTINUUM.EQ.1) CONT=.TRUE.
  ! Consistency check
  CHECK = 0
  IF (NBINS.EQ.5.AND.NTUNE.EQ.5.AND.CONT.EQ..FALSE.) CHECK=CHECK+1
  IF (NBINS.EQ.6.AND.NTUNE.EQ.6.AND.CONT.EQ..FALSE.) CHECK=CHECK+1
  IF (NBINS.EQ.6.AND.NTUNE.EQ.11.AND.CONT.EQ..TRUE.) CHECK=CHECK+1
  IF (CHECK.EQ.0) THEN
     PRINT*,'Illegal combination of Filter parameters: NBINS, NTUNE, CONT'
     PRINT*,'NBIBS:',NBINS
     PRINT*,'NTUNE:',NTUNE
     PRINT*,'CONT:',CONT
     STOP
  ENDIF
  !-------------------------------------------------------------
  ! No continuum point needed
  !-------------------------------------------------------------
  IF (CONT.EQ..FALSE.) THEN
     !-------------------------------------------------------------
     ! Allocating arrays
     !-------------------------------------------------------------
     ALLOCATE (FILTER(NUMW,NBINS),FILTERR(NUMW,NBINS), TUNEPOS(NBINS))
     ALLOCATE (T1(NUMW), LYOT(NUMW))
     !-------------------------------------------------------------
     !-------------------------------------------------------------
     ! Here is where the actual filter profiles are calculated
     !-------------------------------------------------------------
     ! Calculating Lyot filter profile. Assume 1,1,2,2,4,8 periods
!     LYOTW=2D0*1.545*LYOTFWHM
!     T1 = DPI*WAVE/LYOTW
!     LYOT = (COS(T1)**2D0*COS(T1/2D0)**2D0*COS(T1/4D0)**2D0*COS(T1/8D0))**2D0
     ! For the Michelson
!     WMICH(1)=1D0*WNARROW
!     WMICH(2)=2D0*WNARROW
!     WMICH(3)=4D0*WNARROW
     ! Tunning positions
     DO I=1,NBINS
        TUNEPOS(I) = (-(NBINS-1D0)/2D0+DBLE(I-1))*WSPACING
     ENDDO
     ! Getting filter profile
!     DO I=1,NBINS
!        FILTER(:,I) = LYOT
!        DO J=1,NMICH
!           FILTER(:,I)=FILTER(:,I)*COS(DPI*(WAVE+TUNEPOS(I))/(WMICH(J)))**2D0
!        ENDDO
!     ENDDO
     ! Reversing filter profiles in order of increasing wavelengths

!     DO I=1,NBINS
!        DO J=1,NUMW
!           FILTERR(NUMW+1-J,I)=FILTER(J,I)
!        ENDDO
!     ENDDO

!By RCE: Reversing the ORDER of the filters (equivalent to reversing 
!        the spectral line)
!     DO I=1,NBINS
!           FILTERR(:,NBINS+1-I)=FILTER(:,I)
!     ENDDO


!     FILTER=FILTERR
     ! Normalizing filters so they conserve convolved area 
!     DO I=1,NBINS
!        FILTER(:,I)=FILTER(:,I)/SUM(FILTER(:,I))
!     ENDDO
  ENDIF
  !-------------------------------------------------------------
  ! Yes continuum point needed
  !-------------------------------------------------------------
  IF (CONT.EQ..TRUE.) THEN
     !-------------------------------------------------------------
     ! Allocating arrays
     !-------------------------------------------------------------
     ALLOCATE (FILTER(NUMW,NBINS),FILTERR(NUMW,NBINS), TUNEPOS(NBINS))
     ALLOCATE (TEMP1(NUMW,NTUNE), TEMP2(NUMW,NTUNE), TEMP_TUNE(NTUNE))
     ALLOCATE (T1(NUMW), LYOT(NUMW))
     !-------------------------------------------------------------
     !-------------------------------------------------------------
     ! Here is where the actual filter profiles are calculated
     !-------------------------------------------------------------
     ! Calculating Lyot filter profile. Assume 1,1,2,2,4,8 periods
!     LYOTW=2D0*1.545*LYOTFWHM
!     T1 = DPI*WAVE/LYOTW
!     LYOT = (COS(T1)**2D0*COS(T1/2D0)**2D0*COS(T1/4D0)**2D0*COS(T1/8D0))**2D0
     ! For the Michelson
!     WMICH(1)=1D0*WNARROW
!     WMICH(2)=2D0*WNARROW
!     WMICH(3)=4D0*WNARROW
     ! Tunning positions
!     DO I=1,NTUNE
!        TEMP_TUNE(I) = (-(NTUNE-1D0)/2D0+DBLE(I-1))*WSPACING
!     ENDDO
     !PRINT*,TEMP_TUNE
     !STOP
     ! Getting filter profile
!     DO I=1,NTUNE
!        TEMP1(:,I) = LYOT
!        DO J=1,NMICH
!           TEMP1(:,I)=TEMP1(:,I)*COS(DPI*(WAVE+TEMP_TUNE(I))/(WMICH(J)))**2D0
!        ENDDO
!     ENDDO
     ! Reversing filter profiles in order of increasing wavelengths
!     DO I=1,NTUNE
!        DO J=1,NUMW
!           TEMP2(NUMW+1-J,I)=TEMP1(J,I)
!        ENDDO
!     ENDDO

!By RCE: Reversing ORDER of the 6 filter profiles (equivalent 
!        to reversing spectral line)
!     DO I=1,NTUNE
!         !DO J=1,NUMW
!           TEMP2(:,NTUNE+1-I)=TEMP1(:,I)
!        !ENDDO
!     ENDDO
      
!     TEMP1=TEMP2
     ! Normalizing filters so they conserve convolved area 
!     DO I=1,NTUNE
!        TEMP1(:,I)=TEMP1(:,I)/SUM(TEMP1(:,I))
!     ENDDO
     ! Adding only part we are interested in
!     FILTER(:,1:5) = TEMP1(:,4:8)
!     FILTER(:,6) = TEMP1(:,1)
!     TUNEPOS(1:5) = TEMP_TUNE(4:8)
!     TUNEPOS(6) = TEMP_TUNE(NTUNE)
     ! Deallocating temporal arrays
     DEALLOCATE(TEMP1,TEMP2,TEMP_TUNE)
     !
  ENDIF
  !-----------------------------------------------------
  ! Deallocating unneeded arrays
  !-----------------------------------------------------
  DEALLOCATE (T1,LYOT,FILTERR)
  !-----------------------------------------------------
END SUBROUTINE FILT_INIT
