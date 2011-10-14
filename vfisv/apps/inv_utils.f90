MODULE INV_UTILS
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !

  ! By RCE (June 8, 2010): in  SUBROUTINE SVDSOL(DMODEL, CONV_FLAG) added 
  ! conv_flag to flag the convergence of the SVD solving routine. If SVD 
  ! doesn't converge, then CONV_FLAG takes value 0. Changes propagate into SVDCMP.f90. 
  ! This flag is an indication for INVERT to set the results of the inversions 
  ! to NaN for a particular pixel.
  !
  ! By RCE, Jan 2011: fixed a bug in the GET_ERROR routine.
  !
  ! By RCE, April 2011: Small fixes in GET_LAMBDA routine: Changed 
  ! increment/decrement sizes, limitted lower value of lambda and fixed a bug.
  !
  ! By RCE, May 20 2011: Define NORMALIZATION routine to define NORM vector in only 
  ! one place. Other routines make calls to NORMALIZATION.
  ! By RCE, May 23 2011: Define LIMITS routine that sets the LOWER and UPPER limits
  ! allowed for the model parameters. FINE_TUNE_MODEL will call this routine.
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
  SUBROUTINE GET_WEIGHT (OBS,WEIGHTS)
    !
    USE LINE_PARAM
    USE CONS_PARAM
    USE FILT_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN),  DIMENSION(NBINS,4) :: OBS
    REAL(DP), DIMENSION(4)       :: WEIGHTS
    !
    ! By RCE, Feb 2011: Now the weights are passed by the wrapper to the code. 
    ! We normalize them to the maximum 
     WEIGHTS(:) = WEIGHTS(:)/MAXVAL(WEIGHTS(:))
    ! By RCE: and we divide them by the noise
     WEIGHTS(:) = WEIGHTS(:)/NOISE

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
    C2=C2+(1D0/NUMFREE_DEG)*(WEIGHTS(I)**2D0)*SUM(((OBS(:,I)-SYN(:,I)))**2D0)
    ENDDO

  END SUBROUTINE GET_CHI2
  !------------------------------------------------------
  !By RCE, May 20, 2011: Defines Normalization vector for the derivatives
  !
  PURE SUBROUTINE NORMALIZATION(NORM, ICONT)    !!NEW!!
      
    USE CONS_PARAM

    IMPLICIT NONE
    REAL(DP), INTENT(OUT),         DIMENSION(10)       :: NORM
    REAL(DP), INTENT(IN)                               :: ICONT

    NORM(:)=(/25D0,90D0,90D0,1D0,50D0,1500D0,1E5_DP,0.5D0*ICONT,0.5D0*ICONT,0.5D0/)

  END SUBROUTINE NORMALIZATION
  !
  !------------------------------------------------------
  !
  !By RCE May 20 2011: We define the lower and upper limits for the model variables
  ! 
  PURE SUBROUTINE LIMITS(LOWER, UPPER, ICONT)     !!NEW!!
  !
    USE CONS_PARAM

    IMPLICIT NONE

    REAL(DP), INTENT(OUT),          DIMENSION(10)         :: LOWER, UPPER
    REAL(DP), INTENT(IN)                                  :: ICONT

    ! Limits for the model parameters. Used in FINE_TUNE_MODEL
    ! Order: eta0, inclination (deg), azimuth (deg), damping, Doppler width, 
    ! field strength (gauss), line-of-sight velocity (cm/s), source function, 
    ! source function gradient, filling factor

    LOWER = (/1D0, 0D0 , 0D0 , 1D-4,5D0,5D0,-7D5, 1.5E-1*ICONT, 1.5E-1*ICONT, 0D0/)
    UPPER  = (/1D3, 180D0,180D0,5D0, 5D2,5D3, 7D5, 1.2D0*ICONT,  1.2D0*ICONT,  1D0/)
 
  END SUBROUTINE LIMITS
  !
  !------------------------------------------------------

  SUBROUTINE NORMALIZE_DSYN(DSYN,ICONT)
    ! 
    ! By RCE, May 20, 2011: We use the NORM vector defined in 
    ! the NORMALIZATION subroutine, by calling NORMALIZATION(NORM)
    !
    USE FILT_PARAM
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT),     DIMENSION(10,NBINS,4) ::  DSYN
    REAL(DP), INTENT(IN)                               ::  ICONT
    REAL(DP),                    DIMENSION(10)         ::  NORM
    INTEGER                                            ::  I
    !
    CALL NORMALIZATION(NORM, ICONT)
    !
    DO I=1,10
       DSYN(I,:,:)=DSYN(I,:,:)*NORM(I)
    ENDDO
  END SUBROUTINE NORMALIZE_DSYN
  !
  !-------------------------------------------------------
  !
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
  !
  !--------------------------------------------------------
  !
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
       ! By RCE march 2011: If lambda gets too small, then we shouldn't decrease it
       ! IF (LAMBDA_OLD.LE.1E-4) LAMBDA_NEW=LAMBDA_OLD/2D0
       IF (LAMBDA_OLD.LE.1E-4) LAMBDA_NEW=1E-4
       IF (LAMBDA_OLD.LT.1E4.AND.LAMBDA_OLD.GT.1E-4) LAMBDA_NEW=LAMBDA_OLD/5D0
       IF (LAMBDA_OLD.GE.1E4) LAMBDA_NEW=LAMBDA_OLD/10D0
    ENDIF
    !
    IF (IMPROVE.EQ..FALSE.) THEN
       ! We start from largest lambda because we are going to reduce it
       ! and therefore it will not go through 2 ifs.
       IF (LAMBDA_OLD.GE.1E4) LAMBDA_NEW=1E4     
       IF (LAMBDA_OLD.LT.1E4.AND.LAMBDA_OLD.GT.1E-4) LAMBDA_NEW=5D0*LAMBDA_OLD
       ! By RCE: Bug in next line. It used to say .GE. rather than .LE.
       !IF (LAMBDA_OLD.GE.1E-4) LAMBDA_NEW=100D0*LAMBDA_OLD
       IF (LAMBDA_OLD.LE.1E-4) LAMBDA_NEW=10D0*LAMBDA_OLD
    ENDIF
  END SUBROUTINE GET_LAMBDA
  !
  !--------------------------------------------------------
  !
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
          DIVC(I)=DIVC(I)+(2D0/NUMFREE_DEG)*(WEIGHTS(J)**2D0) &
               *SUM((OBS(:,J)-SYN(:,J))*DSYN(FREELOC(I),:,J))
       ENDDO
    ENDDO
  END SUBROUTINE GET_DIVC
  !
  !--------------------------------------------------------
  !
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
    !CALL RANDOM_SEED()

    DO I=1,10
       RAN(I)=NORMAL(0D0,1D0)
       IF (FREE(I).EQ..TRUE.) MODELR(I)=MODELR(I)+RAN(I)*JUMP*MODELR(I)/100D0
    ENDDO
	
    ! By RCE, March 2011: Random jumps for the angles should be distributed 
    ! evenly between 0 and 180. CHECK THIS because I believe it's not right.
    MODELR(2) = 90D0 + 180D0*NORMAL(0D0,1D0)
    MODELR(3) = 90D0 + 180D0*NORMAL(0D0,1D0)

    ! Checking that perturbation did not go too far
    CALL FINE_TUNE_MODEL(MODELR,ICONT)

  END SUBROUTINE RANDOM_MODEL_JUMP
  !
  !--------------------------------------------------------
  !
   SUBROUTINE FINE_TUNE_MODEL(MODEL,ICONT)
    USE INV_PARAM
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT), DIMENSION(10)    :: MODEL
    REAL(DP), DIMENSION(10)                   :: LOWER, UPPER
    REAL(DP), INTENT(IN)                      :: ICONT
    INTEGER                                   :: REV, I

    ! Call routine that sets lower and upper limits for model parameters
    CALL LIMITS(LOWER, UPPER, ICONT)        !!NEW!!

    ! Check all the limits except the angles and non-free parameters
    DO I = 1, 10
      IF ((FREE(I).EQ..TRUE.) .AND. (I .NE. 2) .AND. (I .NE. 3)) THEN
        IF (MODEL(I) .LT. LOWER(I)) MODEL(I) = LOWER(I)        !!NEW!!
        IF (MODEL(I) .GT. UPPER(I)) MODEL(I) = UPPER(I)        !!NEW!!
      ENDIF
    ENDDO

   ! Additional constraint on the source function and its gradient. This
   ! is to ensure that the sum of both quantities stays close to the  
   ! value of the continuum intensity.

IF ((MODEL(8)+MODEL(9)) .LT. 0.8 * ICONT) MODEL(8) = 0.8*ICONT - MODEL(9)        !!NEW!!
IF ((MODEL(8)+MODEL(9)) .GT. 1.2*ICONT) MODEL(8) = 1.2*ICONT-MODEL(9)            !!NEW!!

    ! Check the angles independently so that they vary between 0 and 180.
 
    ! All negatives Inclination treated here to be positive and between 0,PI
    IF (MODEL(2).LT.-360D0) THEN
       REV=ABS((INT(MODEL(2)/360D0)))
       MODEL(2)=MODEL(2)+360D0*REV
    ENDIF
    IF (MODEL(2).LT.0D0 .AND. MODEL(2).GT.-90D0) MODEL(2)=-MODEL(2)
    IF (MODEL(2).LT.-90D0 .AND. MODEL(2).GT.-180D0) MODEL(2)=-MODEL(2)
    IF (MODEL(2).LT.-180 .AND. MODEL(2).GT.-270D0) MODEL(2)=360D0+MODEL(2)
    IF (MODEL(2).LT.-270D0 .AND. MODEL(2).GT. -360D0) MODEL(2)=360D0+MODEL(2)
    ! All positive inclination larger than PI or 2PI put between 0,PI
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
 
  END SUBROUTINE FINE_TUNE_MODEL
  !
  !--------------------------------------------------------
  !
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
             HESS(I,J)=HESS(I,J)+(2D0/NUMFREE_DEG)*(WEIGHTS(K)**2D0) &
                  *SUM(DSYN(FREELOC(I),:,K)*DSYN(FREELOC(J),:,K))
          ENDDO
          IF (I.EQ.J) HESS(I,J)=HESS(I,J)*(1D0+LAMBDA)
       ENDDO
    ENDDO

  END SUBROUTINE GET_HESS
  !
  !------------------------------------------------------------
  !
  SUBROUTINE SVDSOL(DMODEL, CONV_FLAG)
    !
    USE INV_PARAM
    USE CONS_PARAM
    USE FILT_PARAM
    USE SVD_PARAM
    IMPLICIT NONE

    REAL(DP), DIMENSION(15*NUMFREE_PARAM)                   :: WORK
    REAL(DP),       DIMENSION(NUMFREE_PARAM,NUMFREE_PARAM)  :: VT, U, WINV, HINV, V
    REAL(DP),       DIMENSION(NUMFREE_PARAM)                :: PLUSMODEL, W
    REAL(DP)                                                :: WMAX, DIAG
    INTEGER                                                 :: I, J, INFO
    REAL(DP), DIMENSION(10)                                 :: DMODEL
    INTEGER                                                 :: CONV_FLAG
    !
    ! By RCE: commented out all calls to Numerical Recipes routines because of 
    ! distribution issues. I will stick to LAPACK.
    ! Now call Numerical Recipes SVD
    !
    !    CALL SVDCMP(HESS,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM,W,V,CONV_FLAG)
    !
    ! By RCE: call LAPACK SVD
    CALL DGESVD('A','A',NUMFREE_PARAM,NUMFREE_PARAM,HESS,NUMFREE_PARAM,W,U &
   ,NUMFREE_PARAM,VT,NUMFREE_PARAM,WORK,15*NUMFREE_PARAM,INFO)

    IF (INFO.NE.0) THEN
       PRINT*,'ERROR IN SVD'
       CONV_FLAG = 1
    ENDIF

    IF (CONV_FLAG.EQ.1) PRINT*, "CONV_FLAG EQ 1. ERROR IN DGESVD"
   
    IF (CONV_FLAG.NE.1) THEN
    ! By RCE: transpose VT matrix to obtain V and use it in SVBKSB
    DO I = 1, NUMFREE_PARAM
	DO J = 1, NUMFREE_PARAM
		V(i,j) = VT(j,i)
	ENDDO
    ENDDO

    WMAX=MAXVAL(W)
    DO I=1,NUMFREE_PARAM
       IF (W(I).LT.SVDTOL*WMAX) W(I)=0D0
    ENDDO

    CALL SVBKSB(U,W,V,NUMFREE_PARAM,NUMFREE_PARAM,NUMFREE_PARAM, &
         NUMFREE_PARAM,DIVC,PLUSMODEL)
    !
    DO I=1,NUMFREE_PARAM
       DMODEL(FREELOC(I))=PLUSMODEL(I)
    ENDDO
    ENDIF
    !
  END SUBROUTINE SVDSOL
  !
  !------------------------------------------------------------
  !
  PURE SUBROUTINE NORMALIZE_DMODEL(DMODEL,ICONT)
    !
    ! By RCE: May 20 2011: Call NORMALIZATION subroutine for NORM vector.
    !
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT), DIMENSION(10)         ::  DMODEL
    REAL(DP), INTENT(IN)                           ::  ICONT
    REAL(DP),     DIMENSION(10)                    ::  NORM
    INTEGER                                        ::  I
    !
    CALL NORMALIZATION(NORM, ICONT)
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
    REAL(DP), INTENT(IN)                             ::  ICONT
    REAL(DP), DIMENSION(10)                          ::  SIGNO, LIMIT
    INTEGER                                          ::  I
    !
    ! I have tried with more generous limits at this is the best combination 
    ! I could find
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
    ! By RCE May 20 2011: Call NORMALIZATION routine to get NORM vector
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

    ! By RCE: call normalization vector to re-normalize derivatives in Hessian.
    CALL NORMALIZATION(NORM, ICONT)
 
    ! Use LAPACK for calculating SVD of the HESSIAN
    ! We cannot use SVD from Numerical Recipes because SVDCMP destroys the contents
    ! of the U and V matrixes.
    W(:)=0D0
    !
    CALL DGESVD('A','A',NUMFREE_PARAM,NUMFREE_PARAM,HESS,NUMFREE_PARAM, &
         W,U,NUMFREE_PARAM,VT,NUMFREE_PARAM,WORK,15*NUMFREE_PARAM,INFO)
    IF (INFO.NE.0) THEN
       PRINT*,'ERROR IN SVD'
       STOP
    ENDIF
    !

    DO I=1,NUMFREE_PARAM
       DO J=1,NUMFREE_PARAM

    ! By RCE, Jan 2011: Switching indices in what Juanma used 
    !          COV(I,J)=(CHI2/DBLE(NUMFREE_DEG))*SUM(U(:,I)*U(:,J)/W(:))
    !          COV(I,J)=2D0*COV(I,J)*(NORM(FREELOC(I))*NORM(FREELOC(J)))

          COV(I,J)=CHI2/DBLE(NUMFREE_DEG)*SUM(U(I,:)*U(J,:)/W(:))
    ! Re-normalizing the derivatives:
          COV(I,J)=COV(I,J)*(NORM(FREELOC(I))*NORM(FREELOC(J)))
       ENDDO
    ENDDO
    !
    DO I=1,NUMFREE_PARAM
       DO J=1,NUMFREE_PARAM
          COV_DUMMY(FREELOC(I),FREELOC(J))=COV(I,J)
       ENDDO
    ENDDO
    !--------------------------------
    ! Variances (square root of)
    !--------------------------------
    DO I=1,10
       SIGMA(I)=SQRT(COV_DUMMY(I,I))
    ENDDO
    !--------------------------------
    ! Covariances (normalized)
    !--------------------------------
    ! B-Gamma
    SIGMA(11)=COV_DUMMY(6,2)/SIGMA(6)/SIGMA(2)
    ! B-Phi
    SIGMA(12)=COV_DUMMY(6,3)/SIGMA(6)/SIGMA(3)
    ! Gamma-Phi
    SIGMA(13)=COV_DUMMY(2,3)/SIGMA(2)/SIGMA(3)
    ! B-Alpha
    SIGMA(14)=COV_DUMMY(6,10)/SIGMA(6)/SIGMA(10)
    ! Gamma-Alpha
    SIGMA(15)=COV_DUMMY(2,10)/SIGMA(10)/SIGMA(2)
    ! Phi-Alpha
    SIGMA(16)=COV_DUMMY(3,10)/SIGMA(3)/SIGMA(10)
  END SUBROUTINE GET_ERR
  !-------------------------------------------
  
END MODULE INV_UTILS
!CVSVERSIONINFO "$Id: inv_utils.f90,v 1.4 2011/10/14 17:22:55 keiji Exp $"
