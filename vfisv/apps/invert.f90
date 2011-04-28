SUBROUTINE INVERT (OBS_LONG,SCAT_LONG,GUESS,RES,ERR, FILTERS_LONG, CONVERGENCE_FLAG, WEIGHTS)
  !
  ! J M Borrero
  ! Dec 16, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  ! By RCE, Apr 21, 2010: Passing FILTERS variable to INVERT to get 
  ! Sebastien's filter profiles. Defining FILTERS variable inside INVERT.
  ! Passing FILTERS to SYNTHESIS module.
  ! Element [0][0] in C corresponds to (1,1) in Fortran. Element [1,25] in 
  ! C corresponds to (26,2) in Fortran; hence exchanging the indices in definition
  ! of variable: FILTERS (NUMW, NBINS)
  !
  ! By RCE, Apr 23, 2010: Changed azimuth reference by 90 degrees at the end (after the inversion loop).
  !
  ! By RCE, Apr 23, 2010: Normalize all filters to central filter area hard-coded value:
  ! The integrated areas of the six filters for the central pixel (pixel # 8386559) are:    
  ! 1.43858806880917  1.49139978559615   1.52324167542270  1.52811487224149  1.50984871281028  1.46553486521323
  ! the average area being:	1.49279. We will normalize ALL the filters by this value.
  !
  ! By RCE, Apr 23, 2010: Reverse the wavelength order of the observed Stokes profiles so
  ! that the convention of the filter profiles matches the convention for the observations.
  !
  ! By RCE, Nov 1, 2010: Changed the convergence criterium to CHI2 .LT. GOODFIT*(ICONT/6.5d4) to scale it to the intensity 
  ! of the pixel with respect to disk center. Disk center value is HARD CODED!!!! (6.5D4)
  ! 
  ! By RCE, Nove 8, 2010: Changed convergence criterium again. Now, when the NEWCHI2 is an improvement over the old one, 
  ! and the difference between the two is less than GOODFIT, then we stop iterating. This is, when CHI2 doesn't improve
  ! fast enough, we get out of the loop. The absolute check for CHI2 (CHI2 .LT. GOODFIT) has been eliminated because it
  ! wasn't a good measure of convergence.
  ! 
  ! By RCE Feb 2011: Added the weights for the Stokes profiles as an input parameter rather than having them hard-coded
  ! in VFISV. The call to GET_WEIGHTS now just normalizes the weights and multiplies by the noise.
  ! 
  ! By RCE, April 2011: Changes to the LM algorithm:
  !  - Reset to zero the FLAG that counts consecutive bad iterations when a good iteration is achieved (line
  !  - Change random perturbations to angles inside RANDOM_JUMP routine in INV_UTILS
  !  - Put lower bound on the LM lambda parameter so that it doesn't decrease beyond a value (in GET_LAMBDA in INV_UTILS) 
  !  - Test flag (TMP_FLAG) to count the number of resets of the model parameters. Saved in one of the error variables
  !
  ! By RCE: Re-normalizing filters to different value to account for +/-2A coverage
  !
  ! By RCE April 2011: reset FLAG that counts number of consecutive interations to zero when one successful iteration is achieved.
  ! By RCE April 2011: Included filter hack. FILTERS_LONG are the filters calculated in the full wavelength range.
  ! FILTERS are the filter profiles in the narrow wavelength range, where we do the forward modeling.
  ! We do this to avoid computing the forward model very far into the continuum. It takes too much time. So we do
  ! a quick hack to compute the forward model only in an inner wavelength region and then add the filter contribution 
  ! outside of this region, in the continuum around the spectral line.

  USE FILT_PARAM
  USE CONS_PARAM
  USE LINE_PARAM
  USE INV_PARAM
  USE SVD_PARAM
  USE INV_UTILS
  USE FORWARD
  USE RAN_MOD
  IMPLICIT NONE
  !---------------------------------------------------------------------
  REAL(DP), INTENT(IN),  DIMENSION(NBINS*4)          :: OBS_LONG
  REAL(DP), INTENT(IN),  DIMENSION(NBINS*4)          :: SCAT_LONG
  REAL(DP),              DIMENSION(NUMW_LONG, NBINS) :: FILTERS_LONG
  REAL(DP), INTENT(IN),  DIMENSION(10)               :: GUESS
  REAL(DP),              DIMENSION(NUMW,NBINS)       :: FILTERS
  REAL(DP), INTENT(OUT), DIMENSION(10)               :: RES
  REAL(DP), INTENT(OUT), DIMENSION(12)               :: ERR

  !---------------------------------------------------------------------
  REAL(DP),        DIMENSION(NBINS,4)   :: OBS, SCAT, OBS_REV
  REAL(DP),        DIMENSION(NBINS)     :: INTEG_FILTERS
  REAL(DP),        DIMENSION(4)         :: WEIGHTS
  REAL(DP),        DIMENSION(10)        :: BESTMODEL, MODELG, LASTGOODMODEL, DMODEL
  REAL(DP),        DIMENSION(16)        :: SIGMA
  REAL(DP)                              :: BESTCHI2, LASTGOODCHI2, TOTPOL, NEWCHI2, LASTGOODLAMBDA, ICONT
  REAL(DP),        DIMENSION(ITER)      :: LAMBDA, CHI2
  REAL(DP),        DIMENSION(10,ITER)   :: MODEL
  INTEGER                               :: I, K, M, FLAG, WRONG, LASTGOODITER, BESTITER, NPOINTS
  LOGICAL                               :: DERIVATIVE
  INTEGER                               :: CONV_FLAG, NAN_FLAG, CONVERGENCE_FLAG, EPSILON_FLAG, NRESET
  !---------------------------------------------------------------------
  REAL(DP),     DIMENSION(NBINS,4)      :: SYN, LASTGOODSYN, BESTSYN
  REAL(DP),     DIMENSION(10,NBINS,4)   :: DSYN, LASTGOODDSYN
  !---------------------------------------------------------------------

  CALL SORT_OBS(OBS_LONG,OBS)
  CALL SORT_OBS(SCAT_LONG,SCAT)


  ! By RCE April 2011: Select the region of the spectrum of the filters that we are going to use in the 
  ! synthesis and the region that we are going to integrate assuming that they're in the continuum of the
  ! spectral line. 
  
  ! Difference between NUMW wavelength points for synthesis and NUMW_LONG points in which filters were computed.
  ! We divide by 2 to get the number of wavelength points at each side of the "synthesis range"
  NPOINTS = (NUMW_LONG - NUMW)/2  
  
  ! Select region of filters corresponding to forward modeling wavelength region
  FILTERS(:,:) = FILTERS_LONG(NPOINTS+1:NUMW+NPOINTS,:) 
  ! Integrate the filters in the remaining regions at each side of the forward modeling region.
  DO I = 1, NBINS 
     INTEG_FILTERS(I) = SUM(FILTERS_LONG(1:NPOINTS,I)) + SUM(FILTERS_LONG(NUMW_LONG-NPOINTS:,I))
  ENDDO
  ! Making sure we add no filter integral if we're doing the forward modeling in the full wavelength range.
  IF (NPOINTS .EQ. 0) THEN 
     INTEG_FILTERS(:) = 0D0
  ENDIF
  


  !By RCE, Apr 23, 2010: Normalizing filters

!FILTERS(:,:) = FILTERS(:,:) / 1.49279D0 ! Normalization for +/- 0.648A coverage and 27 mA sampling
FILTERS(:,:) = FILTERS(:,:) /1.5799      ! Normalization for +/- 2A coverage and 27mA sampling

  !By RCE, Apr 23, 2010: Reversing the wavelength order of the observations 
  ! so that they are in increasing wavelength order

  DO k=1,NBINS
     DO m = 1,4
	OBS_REV(k,m) = OBS(NBINS-k+1,m)		
     ENDDO
  ENDDO

  OBS(:,:) = OBS_REV(:,:)

  ! By RCE: CONV_FLAG checks for for convergence of SVD and NAN_FLAG checks for NaNs 
  ! in chisqr. 
  ! EPSILON_FLAG checks the convergence rate (how fast CHI2 diminishes). If it's not 
  ! fast enough (NEWCHI2-OLDCHI2 LT epsilon), the iterative process stops.
  ! CONVERGENCE_FLAG is sent out to the wrapper with different values depending 
  ! on whether the algorithm converges or not and why.
  
  CONV_FLAG = 0
  NAN_FLAG = 0
  CONVERGENCE_FLAG = 0
  EPSILON_FLAG = 0
  
  
  ICONT=MAXVAL(OBS(:,1))
  IF (ICONT .GT. TREIC) THEN
     !---------------------------------------------------------------------------------
     ! This IF statement checks wheter the intensity is large enough to invert the data
     ! Otherwsie the pixel is ignored. This is used to avoid pixels off-the-limb
     !---------------------------------------------------------------------------------
     CALL GET_TOTPOL(OBS,TOTPOL)
     CALL GET_WEIGHT(OBS,WEIGHTS)
     !----------------------------------------------------
     ! Commented by RCE: check wfa_guess.f90 routine later
     !     Initial guess for Doppler velocity only
     !----------------------------------------------------

!     IF (TOTPOL.GT.1E-3) THEN
        CALL WFA_GUESS(OBS/ICONT,LANDA0,GUESS(7))
!     ENDIF

     !----------------------------------------------------
     ! Added by RCE: New initialization of Source Function
     ! used when profiles come in photon counts
     !----------------------------------------------------
     !-----------------------------
     ! Inversion initial values
     !-----------------------------
     MODEL(:,1)=GUESS
     MODEL(8,1)=0.4D0*0.98D0*ICONT
     MODEL(9,1)=0.6D0*0.98D0*ICONT
     BESTMODEL(:)=MODEL(:,1)

     IF (FREE(4).EQ..FALSE.) THEN
        MODEL(4,1)=0.5D0
        BESTMODEL(4)=0.5D0
     ENDIF
     IF (FREE(10).EQ..FALSE.) THEN
        MODEL(10,1)=1D0
        BESTMODEL(10)=1D0
     ENDIF

  ! By RCE: Default values for lambda and chi2 so that I can find out when they are not being updated.

     LAMBDA(:) = 33.3333
     CHI2(:) = 33.3333333
     LASTGOODCHI2=1E24
     LAMBDA(1)=0.1D0
     BESTCHI2=1E24
     LASTGOODITER = 1 
     I=1
     FLAG=0
     WRONG=1
     ! Counter for the number of random resets that the model parameters suffer.
     NRESET=0

     !-----------------------------
     ! Start LOOP iteration
     !-----------------------------
     DO WHILE (I .LT. ITER .AND. EPSILON_FLAG .EQ. 0 .AND. NAN_FLAG .EQ. 0 .AND. CONV_FLAG .EQ. 0)
        !----------------------------------------
        ! Too many consecutive bad iterations
        !------------------------------------
        ! Slightly modified by JM Borrero. Apr 15, 2010
        IF (FLAG.GT.5) THEN
           MODELG=BESTMODEL
           CALL RANDOM_MODEL_JUMP(RANDOM_JUMP*WRONG,MODELG,ICONT)
           MODEL(:,I)=MODELG
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=1E12*WRONG
           FLAG=0
           LAMBDA(I)=10D0*WRONG
           WRONG=WRONG+1
           NRESET = NRESET+1
        ENDIF

        !--------------------------------------
        ! Synthetic profiles but no derivatives
        !--------------------------------------
        DERIVATIVE=.FALSE. 
        CALL SYNTHESIS(MODEL(:,I),SCAT,DERIVATIVE,SYN,DSYN, FILTERS, INTEG_FILTERS)
        !-------------
        ! Getting CHI2
        !-------------
        CALL GET_CHI2(SYN,OBS,WEIGHTS,NEWCHI2)
        ! Checking for NAN in NEWCHI2
        IF (NEWCHI2.EQ.NEWCHI2+1D0) THEN
           PRINT*,'NaN detected in Subroutine GETCHI2'
           NAN_FLAG = 1
        ELSE 
        CHI2(I)=NEWCHI2
        !-------------------
        ! Have we improved ?
        !-------------------
        ! YES ##############
        IF (NEWCHI2.LE.LASTGOODCHI2) THEN
           ! Synthetic profiles and derivatives
           DERIVATIVE=.TRUE.
           CALL SYNTHESIS(MODEL(:,I),SCAT,DERIVATIVE,SYN,DSYN, FILTERS, INTEG_FILTERS)
           ! Normalizing Derivatives
           CALL NORMALIZE_DSYN(DSYN,ICONT)
           ! Set to Zero unneeded derivatives
           CALL ZERO_DSYN(DSYN)
           !
           !CHECKING THE RATE OF CONVERGENCE OF CHI2
           IF ((LASTGOODCHI2 - NEWCHI2) .LT. (GOODFIT*ICONT*GOODFIT*ICONT)) EPSILON_FLAG = 1
           !
           LASTGOODLAMBDA=LAMBDA(I)
           LASTGOODITER=I
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=NEWCHI2
           LASTGOODSYN=SYN
           LASTGOODDSYN=DSYN
           WRONG=1D0
	   ! By RCE: Resetting to zero FLAG that counts consecutive bad iterations
	   FLAG = 0

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
        CALL SVDSOL(DMODEL, CONV_FLAG)

        IF (CONV_FLAG .NE. 1) THEN
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
        CALL NORMALIZE_DMODEL(DMODEL,ICONT)
        CALL CUT_DMODEL(DMODEL,ICONT)
        MODEL(:,I+1)=LASTGOODMODEL+DMODEL
        CALL FINE_TUNE_MODEL(MODEL(:,I+1),ICONT)
        ! New Iteration
        I=I+1
     ENDIF ! Checking for CONV_FLAG for SVD convergence
     END IF ! NaN detected in chi2 routine

     ENDDO ! Main loop

     !---------------
     ! Getting Errors
     !---------------

     IF ((NAN_FLAG .EQ. 0) .AND. (CONV_FLAG .EQ. 0)) THEN
	! Compute non-modified Hessian (lambda = 0)
        CALL GET_HESS(LASTGOODDSYN,0D0,WEIGHTS)
        ! This Hessian is constructed with normalized derivatives
        CALL GET_ERR(BESTCHI2,ICONT,SIGMA)
        ERR(1)=SIGMA(6)  ! Bfield
        ERR(2)=SIGMA(2)  ! Gamma
        ERR(3)=SIGMA(3)  ! Phi
        ERR(4)=SIGMA(7)  ! VLOS
        ERR(5)=SIGMA(10) ! Alpha_mag
        ERR(6)=SIGMA(11) ! Bfield-Gamma
        ERR(7)=SIGMA(12) ! Bfield-Phi
        ERR(8)=SIGMA(13) ! Gamma-Phi
        ERR(9)=SIGMA(14) ! Bfield-Alpha
        ERR(10)=SIGMA(15)! Gamma-Alpha
        ERR(11)=SIGMA(16)! Phi-Alpha
        ERR(12)=BESTCHI2 ! Bestchi2
        !--------------
        ! Final Result
        !--------------
        RES=BESTMODEL
        
        ! By RCE, Apr 23, 2010: Juanma's correction to azimuth (which is rotated by 90 degrees with
        ! respect to the convention that people want.
        RES(3) = RES(3) + 90D0
        IF (RES(3) .GT. 180D0) RES(3) = RES(3) - 180D0
        
      ELSE
           IF (NAN_FLAG .EQ. 1) CONVERGENCE_FLAG = 4
           IF (CONV_FLAG .EQ. 1) CONVERGENCE_FLAG = 5
      ENDIF ! CONV_FLAG and NAN_FLAG are equal to 0.

   ELSE ! Intensity above threshold to invert
        CONVERGENCE_FLAG = 6
   ENDIF
     
     IF (I .EQ. ITER) THEN 
        CONVERGENCE_FLAG = 2
        IF (FLAG .GT. 5) CONVERGENCE_FLAG = 3
     ENDIF
     IF (BESTCHI2 .LT. GOODFIT) CONVERGENCE_FLAG = 0


! CONVERGENCE_FLAG values
!    1: (CHI2old-CHI2new) LT epsilon (converged!)
!    2: Maximum number of iterations reached
!    3: Too many non-consecutive successful iterations
!    4: NaN in CHI2 detected
!    5: NaN in SVD
!    6: Pixel not inverted due to ICONT LT threshold intensity.


END SUBROUTINE INVERT
