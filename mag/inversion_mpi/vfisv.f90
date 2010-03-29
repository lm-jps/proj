SUBROUTINE VFISV (OBS_LONG,SCAT_LONG,GUESS,RES,ERR)
  !
  ! J M Borrero
  ! Dec 16, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  USE FILT_PARAM
  USE CONS_PARAM
  USE LINE_PARAM
  USE INV_PARAM
  USE SVD_PARAM
  USE INV_UTILS
  USE FORWARD
  IMPLICIT NONE
  !---------------------------------------------------------------------
  REAL(DP), INTENT(IN),  DIMENSION(NBINS*4)          :: OBS_LONG
  REAL(DP), INTENT(IN),  DIMENSION(NBINS*4)          :: SCAT_LONG
  REAL(DP), INTENT(IN),  DIMENSION(10)               :: GUESS
  REAL(DP), INTENT(OUT), DIMENSION(10)               :: RES
  REAL(DP), INTENT(OUT), DIMENSION(8)                :: ERR
  !---------------------------------------------------------------------
  REAL(DP),        DIMENSION(NBINS,4)   :: OBS, SCAT
  REAL(DP),        DIMENSION(4)         :: WEIGHTS
  REAL(DP),        DIMENSION(10)        :: BESTMODEL, MODELG, LASTGOODMODEL, DMODEL, SIGMA
  REAL(DP)                              :: BESTCHI2, LASTGOODCHI2, TOTPOL, NEWCHI2, LASTGOODLAMBDA
  REAL(DP),        DIMENSION(ITER)      :: LAMBDA, CHI2
  REAL(DP),        DIMENSION(10,ITER)   :: MODEL
  INTEGER                               :: I, K, M, FLAG, WRONG, LASTGOODITER, BESTITER
  LOGICAL                               :: DERIVATIVE
  !---------------------------------------------------------------------
  REAL(DP),     DIMENSION(NBINS,4)      :: SYN, LASTGOODSYN, BESTSYN
  REAL(DP),     DIMENSION(10,NBINS,4)   :: DSYN, LASTGOODDSYN
  !---------------------------------------------------------------------
  CALL SORT_OBS(OBS_LONG,OBS)
  CALL SORT_OBS(SCAT_LONG,SCAT)
  !
  IF (MAXVAL(OBS(:,1)).GT.1E-2) THEN
     !---------------------------------------------------------------------------------
     ! This IF statement checks wheter the intensity is large enough to invert the data
     ! Otherwsie the pixel is ignored. This is used to avoid pixels off-the=limb
     !---------------------------------------------------------------------------------
     CALL GET_TOTPOL(OBS,TOTPOL)
     CALL GET_WEIGHT(OBS,WEIGHTS)
     !-----------------------------
     ! Inversion initial values
     !-----------------------------
     BESTMODEL(:)=GUESS
     MODEL(:,1)=GUESS
     IF (FREE(4).EQ..FALSE.) THEN
        MODEL(4,1)=0.5D0
        BESTMODEL(4)=0.5D0
     ENDIF
     IF (FREE(10).EQ..FALSE.) THEN
        MODEL(10,1)=1D0
        BESTMODEL(10)=1D0
     ENDIF
     LASTGOODCHI2=1E24
     LAMBDA(1)=0.1D0
     BESTCHI2=1E24
     I=1
     FLAG=0
     WRONG=1
     !-----------------------------
     ! Start LOOP iteration
     !-----------------------------
     DO WHILE (I .LT. ITER .AND. BESTCHI2 .GT. GOODFIT)
        !----------------------------------------
        ! Too many consecutive bad iterations
        !------------------------------------
        IF (FLAG.GE.2.AND.BESTCHI2.GT.10D0.AND.I.GT.ITER/2) THEN
           MODELG=BESTMODEL
           CALL RANDOM_MODEL_JUMP(RANDOM_JUMP*WRONG,MODELG)
           MODEL(:,I)=MODELG
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=1E12*WRONG
           FLAG=0
           LAMBDA(I)=10D0*WRONG
           WRONG=WRONG+1
        ENDIF
        IF (FLAG.GT.3) THEN
           MODELG=BESTMODEL
           CALL RANDOM_MODEL_JUMP(RANDOM_JUMP*WRONG,MODELG)
           MODEL(:,I)=MODELG
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=1E12*WRONG
           FLAG=0
           LAMBDA(I)=10D0*WRONG
           WRONG=WRONG+1
        ENDIF
        !--------------------------------------
        ! Synthetic profiles but no derivatives
        !--------------------------------------
        DERIVATIVE=.TRUE. 
        CALL SYNTHESIS(MODEL(:,I),SCAT,DERIVATIVE,SYN,DSYN)
        !-------------
        ! Getting CHI2
        !-------------
        CALL GET_CHI2(SYN,OBS,WEIGHTS,NEWCHI2)
        ! Checking for NAN in NEWCHI2
        IF (NEWCHI2.EQ.NEWCHI2+1D0) THEN
           PRINT*,'NaN detected in Subroutine GETCHI2'
           STOP
        ENDIF
        CHI2(I)=NEWCHI2
        !-------------------
        ! Have we improved ?
        !-------------------
        ! YES ##############
        IF (NEWCHI2.LT.LASTGOODCHI2) THEN
           IF (LASTGOODCHI2-NEWCHI2.LT.0.1) FLAG=FLAG+1
           IF (LASTGOODCHI2-NEWCHI2.GE.0.1) FLAG=0
           ! Synthetic profiles and derivatives
           DERIVATIVE=.TRUE.
           CALL SYNTHESIS(MODEL(:,I),SCAT,DERIVATIVE,SYN,DSYN)
           ! Normalizing Derivatives
           CALL NORMALIZE_DSYN(DSYN)
           ! Make FREE=.FALSE. derivatives = 0
           CALL ZERO_DSYN(DSYN)
           !
           LASTGOODLAMBDA=LAMBDA(I)
           LASTGOODITER=I
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=NEWCHI2
           LASTGOODSYN=SYN
           LASTGOODDSYN=DSYN
           WRONG=1D0
           CALL GET_LAMBDA(LAMBDA(I),.TRUE.,LAMBDA(I+1))
           IF (NEWCHI2.LT.BESTCHI2) THEN
              BESTITER=I
              BESTMODEL=MODEL(:,I)
              BESTCHI2=NEWCHI2
              BESTSYN=SYN
           ENDIF
        ENDIF
        ! NO ##############
        IF (NEWCHI2.GT.LASTGOODCHI2) THEN
           CALL GET_LAMBDA(LAMBDA(I),.FALSE.,LAMBDA(I+1))
           SYN=LASTGOODSYN
           DSYN=LASTGOODDSYN
           MODEL(:,I)=LASTGOODMODEL
           FLAG=FLAG+1
        ENDIF
        !---------------------------
        ! Getting Divergence of CHI2
        !---------------------------
        CALL GET_DIVC(LASTGOODSYN,OBS,LASTGOODDSYN,WEIGHTS)
        !------------------------------------------------
        ! Getting Hessian matrix
        !------------------------------------------------ 
        CALL GET_HESS(LASTGOODDSYN,LAMBDA(I+1),WEIGHTS)
        !------------------------------------------------------------
        ! SVD of the Hessian matrix to get perturbations to old model
        !------------------------------------------------------------
        CALL SVDSOL(DMODEL)
        !---------------------------------------------------------
        ! Make sure non-inverted parameters are not affected
        !---------------------------------------------------------
        DO M=1,10
           IF (FREE(M).EQ..FALSE.) DMODEL(M)=0D0
        ENDDO
        !---------------------------------------------------
        ! Checking for NaN in DMODEL
        !---------------------------------------------------
        DO M=1,10
           IF (DMODEL(M).EQ.DMODEL(M)+1D0) DMODEL(M)=0D0
           IF (ABS(DMODEL(M)).GT.1E10) THEN 
              DMODEL(M)=0D0
              LAMBDA(I+1)=LAMBDA(I)*10D0
           ENDIF
        ENDDO
        !-----------------------------------------
        ! Normalizing, trimming and tunning DMODEL
        !-----------------------------------------
        CALL NORMALIZE_DMODEL(DMODEL)
        CALL CUT_DMODEL(DMODEL)
        MODEL(:,I+1)=LASTGOODMODEL+DMODEL
        CALL FINE_TUNE_MODEL(MODEL(:,I+1))
        ! New Iteration
        I=I+1
     ENDDO
     !---------------------------------------
     ! End LOOP iteration: Inversion Finished
     !---------------------------------------
  ENDIF
  !---------------
  ! Getting Errors
  !---------------
  RES=BESTMODEL
  CALL GET_HESS(LASTGOODDSYN,LAMBDA(LASTGOODITER),WEIGHTS)
  CALL GET_ERR(BESTCHI2,ERR,SIGMA)
  !
END SUBROUTINE VFISV
