MODULE INV_UTILS
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !

  ! By RCE (June 8, 2010): in  SUBROUTINE SVDSOL(DMODEL, CONV_FLAG) added conv_flag to flag the convergence of the SVD solving routine. If SVD doesn't converge, then CONV_FLAG takes value 0. Changes propagate into SVDCMP.f90. 
  ! This flag is an indication for INVERT to set the results of the inversions to NaN for a particular pixel.
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
    I=OBS(1:NBINS)
    Q=OBS(NBINS+1:2*NBINS)/I
    U=OBS(2*NBINS+1:3*NBINS)/I
    V=OBS(3*NBINS+1:4*NBINS)/I
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
    OBS(:,1) = OBS_LONG(1:NBINS)
    OBS(:,2) = OBS_LONG(NBINS+1:2*NBINS)
    OBS(:,3) = OBS_LONG(2*NBINS+1:3*NBINS)
    OBS(:,4) = OBS_LONG(3*NBINS+1:4*NBINS)
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
    WEIGHTS(:) = 1D0/7D0/NOISE
    WEIGHTS(2:3) = 7D0* WEIGHTS(2:3)
    WEIGHTS(4) = 3D0 * WEIGHTS(4)
    !-----------------------------------------------------------
    ! JM Borrero: Apr 15, 2010
    ! WEIGHTS(1)=(0.8D0*ABS(1.01D0-MAXVAL(OBS(:,1)))+0.2D0)/NOISE
    ! Change weights only if you know what you are doing !!!
    !-----------------------------------------------------------
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
    
    C2=0D0

    DO I=1,4    
    C2=C2+(1D0/NUMFREE_DEG)*(WEIGHTS(I)**2D0)*SUM(((OBS(:,I)-SYN(:,I)))**2D0)/SUM(WEIGHTS)
    ENDDO


  END SUBROUTINE GET_CHI2
  !------------------------------------------------------
  PURE SUBROUTINE NORMALIZE_DSYN(DSYN,ICONT)
    !
    USE FILT_PARAM
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT),     DIMENSION(10,NBINS,4) ::  DSYN
    REAL(DP), INTENT(IN)                               ::  ICONT
    REAL(DP),                    DIMENSION(10)         ::  NORM
    INTEGER                                            ::  I
    !
    NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.5D0/)
    !NORM(:)=(/1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0/)
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
          DIVC(I)=DIVC(I)-(2D0/NUMFREE_DEG)*(WEIGHTS(J)**2D0)*SUM((OBS(:,J)-SYN(:,J))*DSYN(FREELOC(I),:,J))
       ENDDO
    ENDDO
  END SUBROUTINE GET_DIVC
  !--------------------------------------------------------
  SUBROUTINE RANDOM_MODEL_JUMP(JUMP,MODELR,ICONT)
    USE CONS_PARAM
    USE INV_PARAM
    USE RAN_MOD
    IMPLICIT NONE
    REAL(DP), INTENT(IN)                      :: JUMP, ICONT
    REAL(DP), INTENT(INOUT),    DIMENSION(10) :: MODELR
    REAL(DP),                   DIMENSION(10) :: RAN
    INTEGER                                   :: I
    !--------------------------------------------
    ! JM Borrero: Apr 15, 2010
    ! Modified to include and normal distribution
    ! instead of uniform distribution
    !--------------------------------------------
    DO I=1,10
       RAN(I)=NORMAL(0D0,1D0)
       IF (FREE(I).EQ..TRUE.) MODELR(I)=MODELR(I)+RAN(I)*JUMP*MODELR(I)/100D0
    ENDDO
    ! Checking that perturbation did not go too far
    CALL FINE_TUNE_MODEL(MODELR,ICONT)
  END SUBROUTINE RANDOM_MODEL_JUMP
  !--------------------------------------------------------
  PURE SUBROUTINE FINE_TUNE_MODEL(MODEL,ICONT)
    USE INV_PARAM
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT), DIMENSION(10)    :: MODEL
    REAL(DP), INTENT(IN)                      :: ICONT
    INTEGER                                   :: REV
    ! Avoid negative or too small values in thermodynamic parameters
    IF (MODEL(1).LE.1D0) MODEL(1)=1D0
    IF (MODEL(4).LE.1D-4) MODEL(4)=1D-4
    IF (MODEL(5).LE.5D0) MODEL(5)=5D0
    IF (MODEL(8).LE.1.5E-1*ICONT) MODEL(8)=1.5E-1*ICONT
    IF (MODEL(9).LE.1.5E-1*ICONT) MODEL(9)=1.5E-1*ICONT
    ! Avoid too large values in thermodynamic parameter
    IF (MODEL(1).GE.100D0) MODEL(1)=100D0
    IF (MODEL(4).GE.5D0) MODEL(4)=5D0
    IF (MODEL(5).GE.50D0) MODEL(5)=50D0
    IF (MODEL(8).GE.1.5D0*ICONT) MODEL(8)=1.5D0*ICONT
    IF (MODEL(9).GE.1.5D0*ICONT) MODEL(9)=1.5D0*ICONT
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
    !   MODEL(6)=-MODEL(6)
	MODEL(6) = 10D0
    ENDIF
    ! Too large magnetic field
    IF (MODEL(6).GT.5000D0) MODEL(6)=5000D0
    ! Too large velocities
    IF (MODEL(7).GT.5E5) MODEL(7)=7E5
    IF (MODEL(7).LT.-5E5) MODEL(7)=-7E5
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
  SUBROUTINE SVDSOL(DMODEL, CONV_FLAG)
    !
    USE INV_PARAM
    USE CONS_PARAM
    USE FILT_PARAM
    USE SVD_PARAM
    IMPLICIT NONE

    REAL(DP), DIMENSION(15*NUMFREE_PARAM)                  :: WORK
    REAL(DP),       DIMENSION(NUMFREE_PARAM,NUMFREE_PARAM)           :: VT, U, WINV, HINV, V
    REAL(DP),       DIMENSION(NUMFREE_PARAM)                         :: PLUSMODEL, W
    REAL(DP)                                                         :: WMAX, DIAG
    INTEGER                                                          :: I, J, INFO
    REAL(DP), DIMENSION(10)                                          :: DMODEL
    INTEGER                                                          :: CONV_FLAG
    !
    ! Now call Numerical Recipes SVD
    !
!    CALL SVDCMP(HESS,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,W,V,CONV_FLAG)
    !
! By RCE: call LAPACK SVD
    CALL DGESVD('A','A',NUMFREE_PARAM,NUMFREE_PARAM,HESS,NUMFREE_PARAM,W,U,NUMFREE_PARAM &
         ,VT,NUMFREE_PARAM,WORK,15*NUMFREE_PARAM,INFO)

    IF (INFO.NE.0) THEN
       PRINT*,'ERROR IN SVD'
       CONV_FLAG = 1
    ENDIF

IF (CONV_FLAG.EQ.1) PRINT*, "CONV_FLAG EQ 1. ERROR IN DGESVD"
   
    IF (CONV_FLAG.NE.1) THEN
! By RCE: transpose VT matrix to obtain V
    DO I = 1, NUMFREE_PARAM
	DO J = 1, NUMFREE_PARAM
		V(i,j) = VT(j,i)
	ENDDO
    ENDDO

    WMAX=MAXVAL(W)
    DO I=1,NUMFREE_PARAM
       IF (W(I).LT.SVDTOL*WMAX) W(I)=0D0
    ENDDO
    CALL SVBKSB(U,W,V,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,DIVC,PLUSMODEL)
    !
    DO I=1,NUMFREE_PARAM
       DMODEL(FREELOC(I))=PLUSMODEL(I)
    ENDDO
    ENDIF
    !
  END SUBROUTINE SVDSOL
  !------------------------------------------------------------
  PURE SUBROUTINE NORMALIZE_DMODEL(DMODEL,ICONT)
    !
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT), DIMENSION(10)         ::  DMODEL
    REAL(DP), INTENT(IN)                           ::  ICONT
    REAL(DP),     DIMENSION(10)                    ::  NORM
    INTEGER                                        ::  I
    !
    NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.5D0/)
    !NORM(:)=(/1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0/)
    !
    DO I=1,10
       DMODEL(I)=DMODEL(I)*NORM(I)
    ENDDO
    !
  END SUBROUTINE NORMALIZE_DMODEL
  !------------------------------------------------------------
  PURE SUBROUTINE CUT_DMODEL(DMODEL,ICONT)
    !
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT),   DIMENSION(10)         ::  DMODEL
    REAL(DP), INTENT(IN)                             :: ICONT
    REAL(DP), DIMENSION(10)                          ::  SIGNO, LIMIT
    INTEGER                                          ::  I
    !
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
  SUBROUTINE GET_ERR(CHI2,ICONT,SIGMA)
    !
    USE CONS_PARAM
    USE INV_PARAM
    USE SVD_PARAM
    USE LINE_PARAM
    IMPLICIT NONE
    REAL(DP), DIMENSION(NUMFREE_PARAM,NUMFREE_PARAM)       :: VT, U, COV
    REAL(DP), DIMENSION(15*NUMFREE_PARAM)                  :: WORK
    REAL(DP), DIMENSION(NUMFREE_PARAM)                     :: W
    REAL(DP)                                               :: CHI2, ICONT
    REAL(DP), DIMENSION(10)                                :: NORM
    REAL(DP), DIMENSION(10,10)                             :: COV_DUMMY
    REAL(DP), DIMENSION(16)                                :: SIGMA
    INTEGER                                                :: I, J, INFO, NFREE, CONV_FLAG
    !
    NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.5D0/)
    !NORM(:)=(/1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0/)
    ! Use LAPACK for calculating SVD of the HESSIAN
    ! We cannot use SVD from Numerical Recipes because SVDCMP destroys the contents
    ! of the U and V matrixes.
    W(:)=0D0
    !
    CALL DGESVD('A','A',NUMFREE_PARAM,NUMFREE_PARAM,HESS,NUMFREE_PARAM,W,U,NUMFREE_PARAM &
         ,VT,NUMFREE_PARAM,WORK,15*NUMFREE_PARAM,INFO)
    IF (INFO.NE.0) THEN
       PRINT*,'ERROR IN SVD'
       STOP
    ENDIF
    !
    DO I=1,NUMFREE_PARAM
       DO J=1,NUMFREE_PARAM
          COV(I,J)=(CHI2/DBLE(NUMFREE_DEG))*SUM(U(:,I)*U(:,J)/W(:))
          COV(I,J)=2D0*COV(I,J)*(NORM(FREELOC(I))*NORM(FREELOC(J)))
       ENDDO
    ENDDO
    !
    DO I=1,NUMFREE_PARAM
       DO J=1,NUMFREE_PARAM
          COV_DUMMY(FREELOC(I),FREELOC(J))=COV(I,J)
       ENDDO
    ENDDO
    !----------
    ! Variances
    !----------
    DO I=1,10
       SIGMA(I)=COV_DUMMY(I,I)
    ENDDO
    !------------
    ! Covariances
    !------------
    ! B-Gamma
    SIGMA(11)=COV_DUMMY(6,2)
    ! B-Phi
    SIGMA(12)=COV_DUMMY(6,3)
    ! Gamma-Phi
    SIGMA(13)=COV_DUMMY(2,3)
    ! B-Alpha
    SIGMA(14)=COV_DUMMY(6,10)
    ! Gamma-Alpha
    SIGMA(15)=COV_DUMMY(2,10)
    ! Phi-Alpha
    SIGMA(16)=COV_DUMMY(3,10)
  END SUBROUTINE GET_ERR
  !-------------------------------------------
  
END MODULE INV_UTILS
