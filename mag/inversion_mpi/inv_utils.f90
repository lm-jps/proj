MODULE INV_UTILS
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
CONTAINS
  !---------------------------------------------------
  PURE SUBROUTINE GET_TOTPOL(OBS,TOTPOL)
    USE CONS_PARAM
    USE FILT_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(NBINS*4) :: OBS
    REAL(DP), INTENT(OUT)                    :: TOTPOL
    REAL(DP), DIMENSION(NBINS)               :: I, Q, U, V
    !
    I=OBS(1:6)
    Q=OBS(7:12)/I
    U=OBS(13:18)/I
    V=OBS(19:24)/I
    !
    TOTPOL=SUM(SQRT(Q**2D0+U**2D0+V**2D0))
  END SUBROUTINE GET_TOTPOL
  !----------------------------------------------------
  PURE SUBROUTINE SORT_OBS(OBS_LONG,OBS)
    USE FILT_PARAM
    IMPLICIT NONE
    REAL(DP),  INTENT(IN),  DIMENSION(NBINS*4)      :: OBS_LONG
    REAL(DP),  INTENT(OUT), DIMENSION(NBINS,4)      :: OBS
    !
    OBS(:,1) = OBS_LONG(1:6)
    OBS(:,2) = OBS_LONG(7:12)
    OBS(:,3) = OBS_LONG(13:18)
    OBS(:,4) = OBS_LONG(19:24)
    !
  END SUBROUTINE SORT_OBS

  !----------------------------------------------------
  PURE SUBROUTINE GET_WEIGHT (OBS,WEIGHTS)
    !
    USE LINE_PARAM
    USE CONS_PARAM
    USE FILT_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN),  DIMENSION(NBINS,4) :: OBS
    REAL(DP), INTENT(OUT), DIMENSION(4)       :: WEIGHTS
    !
    WEIGHTS(:)=1D0/NOISE
    WEIGHTS(1)=8E-1*WEIGHTS(1)
    WEIGHTS(2) = 2D0*WEIGHTS(2)
    WEIGHTS(3) = 2D0*WEIGHTS(3)
    !WEIGHTS(1)=(0.8D0*ABS(1.01D0-MAXVAL(OBS(:,1)))+0.2D0)/NOISE
    !
  END SUBROUTINE GET_WEIGHT
  !----------------------------------------------------
  PURE SUBROUTINE GET_CHI2(SYN,OBS,WEIGHTS,C2)
    !
    USE INV_PARAM
    USE CONS_PARAM
    USE FILT_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN),   DIMENSION(4)       :: WEIGHTS
    REAL(DP), INTENT(IN),   DIMENSION(NBINS,4) :: OBS, SYN
    REAL(DP), INTENT(OUT)                      :: C2
    INTEGER                                    :: I
    !
    C2=0D0
    
    DO I=1,4
       C2=C2+(1D0/NUMFREE_DEG)*WEIGHTS(I)**2D0*SUM((OBS(:,I)-SYN(:,I))**2D0)
    ENDDO
  END SUBROUTINE GET_CHI2
  !------------------------------------------------------
  PURE SUBROUTINE NORMALIZE_DSYN(DSYN)
    !
    USE FILT_PARAM
    USE CONS_PARAM
    USE INV_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT),     DIMENSION(10,NBINS,4) ::  DSYN
    REAL(DP),                    DIMENSION(10)         ::  NORM
    REAL(DP)                                           ::  ICONT
    INTEGER                                            ::  I
    !
    ICONT = 1D3*TREIC
    NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.5D0/)
    !
    
    DO I=1,10
       DSYN(I,:,:)=DSYN(I,:,:)*NORM(I)
    ENDDO
  END SUBROUTINE NORMALIZE_DSYN
  !-------------------------------------------------------
  PURE SUBROUTINE ZERO_DSYN (DSYN)
    USE CONS_PARAM
    USE FILT_PARAM
    USE INV_PARAM
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT),  DIMENSION(10,NBINS,4)  :: DSYN
    INTEGER                                           :: I
    DO I=1,10
       IF (FREE(I).EQ..FALSE.) THEN 
          DSYN(I,:,:)=0D0
       ENDIF
    ENDDO
    !
  END SUBROUTINE ZERO_DSYN
  !--------------------------------------------------------
  PURE SUBROUTINE GET_LAMBDA(LAMBDA_OLD,IMPROVE,LAMBDA_NEW)
    !
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN)           :: LAMBDA_OLD
    REAL(DP), INTENT(OUT)          :: LAMBDA_NEW
    LOGICAL, INTENT(IN)            :: IMPROVE
    !
    IF (IMPROVE.EQ..TRUE.) THEN
       ! We start from smallest lambda because we are going to reduce it
       ! and therefore it will not go through 2 ifs.
       IF (LAMBDA_OLD.LE.1E-4) LAMBDA_NEW=LAMBDA_OLD/2D0       
       IF (LAMBDA_OLD.LT.1E4.AND.LAMBDA_OLD.GT.1E-4) LAMBDA_NEW=LAMBDA_OLD/10D0
       IF (LAMBDA_OLD.GE.1E4) LAMBDA_NEW=LAMBDA_OLD/100D0
    ENDIF
    !
    IF (IMPROVE.EQ..FALSE.) THEN
       ! We start from largest lambda because we are going to reduce it
       ! and therefore it will not go through 2 ifs.
       IF (LAMBDA_OLD.GE.1E4) LAMBDA_NEW=2D0*LAMBDA_OLD       
       IF (LAMBDA_OLD.LT.1E4.AND.LAMBDA_OLD.GT.1E-4) LAMBDA_NEW=10D0*LAMBDA_OLD
       IF (LAMBDA_OLD.GE.1E-4) LAMBDA_NEW=100D0*LAMBDA_OLD
    ENDIF
  END SUBROUTINE GET_LAMBDA
  !--------------------------------------------------------
  SUBROUTINE GET_DIVC(SYN,OBS,DSYN,WEIGHTS)
    !
    USE CONS_PARAM
    USE FILT_PARAM 
    USE INV_PARAM
    USE SVD_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(NBINS,4)      :: OBS, SYN
    REAL(DP), INTENT(IN), DIMENSION(10,NBINS,4)   :: DSYN
    REAL(DP), INTENT(IN), DIMENSION(4)            :: WEIGHTS
    INTEGER                                       :: I, J
    !
    DIVC(:)=0D0
    DO I=1,NUMFREE_PARAM
       DO J=1,4
          DIVC(I)=DIVC(I)-(2D0/NUMFREE_DEG)*WEIGHTS(J)**2D0*SUM((OBS(:,J)-SYN(:,J))*DSYN(FREELOC(I),:,J))
       ENDDO
    ENDDO
  END SUBROUTINE GET_DIVC
  !--------------------------------------------------------
  SUBROUTINE RANDOM_MODEL_JUMP(JUMP,MODELR)
    USE CONS_PARAM
    USE INV_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN)                      :: JUMP
    REAL(DP), INTENT(INOUT),    DIMENSION(10) :: MODELR
    REAL(DP),                   DIMENSION(10) :: RAN
    INTEGER                                   :: I
    !
    CALL RANDOM_SEED
    CALL RANDOM_NUMBER(RAN)
    DO I=1,10
       IF (FREE(I).EQ..TRUE.) MODELR(I)=(1D0-JUMP/100D0)*MODELR(I)+(2D0*JUMP/100D0)*MODELR(I)*RAN(I)
    ENDDO
    ! Checking that perturbation did not go too far
    CALL FINE_TUNE_MODEL(MODELR)
  END SUBROUTINE RANDOM_MODEL_JUMP
  !--------------------------------------------------------
  PURE SUBROUTINE FINE_TUNE_MODEL(MODEL)
    USE INV_PARAM
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT), DIMENSION(10)    :: MODEL
    INTEGER                                   :: REV
    ! Avoid negative or too small values in thermodynamic parameters
    IF (MODEL(1).LE.1D0) MODEL(1)=1D0
    IF (MODEL(4).LE.1D-4) MODEL(4)=1D-4
    IF (MODEL(5).LE.5D0) MODEL(5)=5D0
    IF (MODEL(8).LE. 0.5D0*TREIC) MODEL(8)=0.5D0*TREIC
    IF (MODEL(9).LE. 0.5D0*TREIC) MODEL(9)=0.5D0*TREIC
    ! Avoid too large values in thermodynamic parameter
    IF (MODEL(1).GE.100D0) MODEL(1)=100D0
    IF (MODEL(4).GE.5D0) MODEL(4)=5D0
    IF (MODEL(5).GE.100D0) MODEL(5)=100D0
    IF (MODEL(8).GE.8E5) MODEL(8)=8E5
    IF (MODEL(9).GE.8E5) MODEL(9)=8E5
    ! All negatives Gamma treated here to be positive and between 0,PI
    IF (MODEL(2).LT.-360D0) THEN
       REV=ABS((INT(MODEL(2)/360D0)))
       MODEL(2)=MODEL(2)+360D0*REV
    ENDIF
    IF (MODEL(2).LT.0D0 .AND. MODEL(2).GT.-90D0) MODEL(2)=-MODEL(2)
    IF (MODEL(2).LT.-90D0 .AND. MODEL(2).GT.-180D0) MODEL(2)=-MODEL(2)
    IF (MODEL(2).LT.-180 .AND. MODEL(2).GT.-270D0) MODEL(2)=360D0+MODEL(2)
    IF (MODEL(2).LT.-270D0 .AND. MODEL(2).GT. -360D0) MODEL(2)=360D0+MODEL(2)
    ! All positive Gamma larger than PI or 2PI put between 0,PI
    IF (MODEL(2).GT.360D0) THEN
       REV=INT(MODEL(2)/360D0)
       MODEL(2)=MODEL(2)-360D0*REV
    ENDIF
    IF (MODEL(2).GT.180) MODEL(2)=360D0-MODEL(2)
    ! Shift azimuthal angles
    IF (MODEL(3).GT.360D0) THEN
       REV=INT(MODEL(3)/360D0)
       MODEL(3)=MODEL(3)-360D0*REV
    ENDIF
    IF (MODEL(3).GT.180) MODEL(3)=MODEL(3)-180D0
    IF (MODEL(3).LT.-360D0) THEN
       REV=ABS((INT(MODEL(3)/360D0)))
       MODEL(3)=MODEL(3)+360D0*REV
    ENDIF
    IF (MODEL(3).LT.-180) MODEL(3)=MODEL(3)+180D0
    IF (MODEL(3).LT.0) MODEL(3)=MODEL(3)+180D0
    IF (MODEL(6).LT.0D0) THEN
       MODEL(6)=-MODEL(6)
    ENDIF
    ! Too large magnetic field
    IF (MODEL(6).GT.5000D0) MODEL(6)=5000D0
    ! Too large velocities
    IF (MODEL(7).GT.7E5) MODEL(7)=7E5
    IF (MODEL(7).LT.-7E5) MODEL(7)=-7E5
    ! Filling factor negative or larger than 1 because of MODEL
    IF (FREE(10).EQ..TRUE.) THEN
       IF (MODEL(10).LT.0D0) MODEL(10)=0D0
       IF (MODEL(10).GT.1D0) MODEL(10)=1D0
    ENDIF
  END SUBROUTINE FINE_TUNE_MODEL
  !--------------------------------------------------------
  SUBROUTINE GET_HESS(DSYN,LAMBDA,WEIGHTS)
    !
    USE CONS_PARAM
    USE FILT_PARAM
    USE INV_PARAM
    USE SVD_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(10,NBINS,4) :: DSYN
    REAL(DP), INTENT(IN), DIMENSION(4)          :: WEIGHTS
    REAL(DP), INTENT(IN)                        :: LAMBDA
    INTEGER                                     :: I, J, K
    !
    HESS(:,:)=0D0
    !
    DO I=1,NUMFREE_PARAM
       DO J=1,NUMFREE_PARAM
          DO K=1,4
             HESS(I,J)=HESS(I,J)-(2D0/NUMFREE_DEG)*(WEIGHTS(K)**2D0)*SUM(DSYN(FREELOC(I),:,K) &
                  *DSYN(FREELOC(J),:,K))
          ENDDO
          IF (I.EQ.J) HESS(I,J)=HESS(I,J)*(1D0+LAMBDA)
       ENDDO
    ENDDO
  END SUBROUTINE GET_HESS
  !------------------------------------------------------------
  SUBROUTINE SVDSOL(DMODEL)
    !
    USE INV_PARAM
    USE CONS_PARAM
    USE FILT_PARAM
    USE SVD_PARAM
    IMPLICIT NONE
    REAL(DP),       DIMENSION(NUMFREE_PARAM,NUMFREE_PARAM)           :: VT, U, WINV, HINV, V
    REAL(DP),       DIMENSION(NUMFREE_PARAM)                         :: PLUSMODEL, W
    REAL(DP)                                                         :: WMAX, DIAG
    INTEGER                                                          :: I, J, INFO
    REAL(DP), DIMENSION(10)                                          :: DMODEL
    !
    ! Now call Numerical Recipes SVD
    CALL SVDCMP(HESS,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,W,V)
    !
    WMAX=MAXVAL(W)
    DO I=1,NUMFREE_PARAM
       IF (W(I).LT.SVDTOL*WMAX) W(I)=0D0
    ENDDO
    CALL SVBKSB(HESS,W,V,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,DIVC,PLUSMODEL)
    !
    DO I=1,NUMFREE_PARAM
       DMODEL(FREELOC(I))=PLUSMODEL(I)
    ENDDO
    !
  END SUBROUTINE SVDSOL
  !------------------------------------------------------------
  PURE SUBROUTINE NORMALIZE_DMODEL(DMODEL)
    !
    USE CONS_PARAM
    USE INV_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT), DIMENSION(10)         ::  DMODEL
    REAL(DP),     DIMENSION(10)                    ::  NORM
    REAL(DP)                                       ::  ICONT
    INTEGER                                        ::  I
    !
    ICONT = 1D3 * TREIC
    NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*ICONT,0.3D0*ICONT,0.5D0/)
    !
    DO I=1,10
       DMODEL(I)=DMODEL(I)*NORM(I)
    ENDDO
    !
  END SUBROUTINE NORMALIZE_DMODEL
  !------------------------------------------------------------
  PURE SUBROUTINE CUT_DMODEL(DMODEL)
    !
    USE CONS_PARAM
    USE INV_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT),   DIMENSION(10)         ::  DMODEL
    REAL(DP),     DIMENSION(10)                      ::  SIGNO, LIMIT
    REAL(DP)                                         ::  ICONT
    INTEGER                                          ::  I
    !
    ICONT = 1D3*TREIC
    ! I have tried with more generous limits at this is the best combination I could find
    LIMIT(:)=(/10D0,25D0,25D0,0.1D0,30D0,500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.25D0/)
    !
    DO I=1,10
       IF (DMODEL(I).GT.0D0) SIGNO(I)=1D0
       IF (DMODEL(I).LT.0D0) SIGNO(I)=-1D0
       IF (ABS(DMODEL(I)).GT.LIMIT(I)) DMODEL(I)=SIGNO(I)*LIMIT(I)
    ENDDO
    !   
  END SUBROUTINE CUT_DMODEL
  !-------------------------------------------
  SUBROUTINE GET_ERR(CHI2,ERR,SIGMA)
    !
    USE CONS_PARAM
    USE INV_PARAM
    USE SVD_PARAM
    IMPLICIT NONE
    REAL(DP), DIMENSION(NUMFREE_PARAM,NUMFREE_PARAM)       :: VT, U, WINV, HINV, V
    REAL(DP), DIMENSION(15*NUMFREE_PARAM)                  :: WORK
    REAL(DP), DIMENSION(NUMFREE_PARAM)                     :: W
    REAL(DP)                                               :: WMAX, CHI2
    REAL(DP), DIMENSION(10)                                :: SIGMA, DSIGMA, NORM
    REAL(DP), DIMENSION(11)                                 :: ERR
    INTEGER                                                :: I, INFO, NFREE
    !
    NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*1.5E6,0.3D0*1.5E6,0.5D0/)
    ! Covariance BFIELD-INCLINATION
    ERR(6) = SQRT(CHI2/DBLE(NUMFREE_DEG))*HESS(6,2)*SQRT(NORM(2))*SQRT(NORM(6))
    ! Covariance BFIELD-AZIMUTH
    ERR(7) = SQRT(CHI2/DBLE(NUMFREE_DEG))*HESS(6,3)*SQRT(NORM(3))*SQRT(NORM(6))
    ! Covariance INCLINATION-AZIMUTH
    ERR(8) = SQRT(CHI2/DBLE(NUMFREE_DEG))*HESS(2,3)*SQRT(NORM(2))*SQRT(NORM(3))
    ! Covariance FIELD-ALPHA
    ERR(9) = SQRT(CHI2/DBLE(NUMFREE_DEG))*HESS(6,10)*SQRT(NORM(6))*SQRT(NORM(10))
    ! Covariance INCLINATION-ALPHA
    ERR(10) = SQRT(CHI2/DBLE(NUMFREE_DEG))*HESS(2,10)*SQRT(NORM(2))*SQRT(NORM(10))
    ! Covariance AZIMUTH-ALPHA
    ERR(11) = SQRT(CHI2/DBLE(NUMFREE_DEG))*HESS(3,10)*SQRT(NORM(3))*SQRT(NORM(10))

    ! Use LAPACK for calculating SVD of the HESSIAN
    ! We cannot use SVD from Numerical Recipes because SVDCMP destroys the contents
    ! of the U and V matrixes.
    W(:)=0D0
    DSIGMA(:)=0D0
    !
    CALL DGESVD('A','A',NUMFREE_PARAM,NUMFREE_PARAM,HESS,NUMFREE_PARAM,W,U,NUMFREE_PARAM &
         ,VT,NUMFREE_PARAM,WORK,15*NUMFREE_PARAM,INFO)
    ! WRITE(*,*) 'INFO is equal to', INFO
    IF (INFO.NE.0) THEN
       PRINT*,'ERROR IN SVD'
       STOP
    ENDIF 
    WMAX=MAXVAL(W)
    WINV(:,:)=0D0
    HINV(:,:)=0D0
    DO I=1,NUMFREE_PARAM
       IF (W(I).LT.SVDTOL*WMAX) WINV(I,I)=0D0
       IF (W(I).GE.SVDTOL*WMAX) WINV(I,I)=1D0/W(I)
    ENDDO
    !
    HINV=-1D0*MATMUL(MATMUL(TRANSPOSE(VT),WINV),TRANSPOSE(U))
    ! 
    DO I=1,NUMFREE_PARAM
       ! This is basically Equation 11.29 in book by Jose Carlos del Toro
       SIGMA(FREELOC(I))=SQRT((CHI2/DBLE(NUMFREE_DEG))*HINV(I,I))
       SIGMA(FREELOC(I))=SIGMA(FREELOC(I))*NORM(FREELOC(I))
    ENDDO
    !
    ERR(1) = SIGMA(6) ! Field strength
    ERR(2) = SIGMA(2) ! Field inclination
    ERR(3) = SIGMA(3) ! Field azimuth
    ERR(4) = SIGMA(7) ! LOS-velocity
    ERR(5) = SIGMA(10)! Alpha_mag
  END SUBROUTINE GET_ERR
  !-------------------------------------------
 
END MODULE INV_UTILS
