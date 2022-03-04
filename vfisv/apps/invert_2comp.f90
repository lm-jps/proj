SUBROUTINE INVERT (OBS_LONG,SCAT_LONG,GUESS,RES,ERR, FILTERS_LONG, CONVERGENCE_FLAG, WEIGHTS, POL_LOW, POL_UPP)  !ABGM
!!!SUBROUTINE INVERT (OBS_LONG,SCAT_LONG,GUESS,RES,ERR, FILTERS_LONG, CONVERGENCE_FLAG, WEIGHTS)
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
  ! By RCE, Apr 23, 2010: (OBSOLETE, see May 17 2011 comment) Normalize all filters to central filter area hard-coded value:
  ! The integrated areas of the six filters for the central pixel (pixel # 8386559) are:    
  ! 1.43858806880917  1.49139978559615   1.52324167542270  1.52811487224149  1.50984871281028  1.46553486521323
  ! the average area being: 1.49279. We will normalize ALL the filters by this value.
  !
  ! By RCE, Apr 23, 2010: Reverse the wavelength order of the observed Stokes profiles so
  ! that the convention of the filter profiles matches the convention for the observations.
  !
  ! By RCE, Nov 1, 2010: OBSOLETE. See March 2012 comment. 
  ! Changed the convergence criterium to CHI2 .LT. GOODFIT*(ICONT/6.5d4) to scale it to the intensity 
  ! of the pixel with respect to disk center. Disk center value is HARD CODED!!!! (6.5D4)
  ! 
  ! By RCE, Nove 8, 2010:  OBSOLETE: See March 2012 comment.
  ! Changed convergence criterium again. Now, when the NEWCHI2 is an improvement over the old one, 
  ! and the difference between the two is less than GOODFIT, then we stop iterating. This is, when CHI2 doesn't improve
  ! fast enough, we get out of the loop. The absolute check for CHI2 (CHI2 .LT. GOODFIT) has been eliminated because it
  ! wasn't a good measure of convergence.
  ! 
  ! By RCE Feb 2011: Added the weights for the Stokes profiles as an input parameter rather than having them hard-coded
  ! in VFISV. The call to GET_WEIGHTS now just normalizes the weights and multiplies by the noise.
  !
  ! By RCE: Re-normalizing filters to different value to account for +/-2A coverage (OBSOLETE, see May 17 2011 comment)
  !
  ! By RCE April 2011: Included filter hack. FILTERS_LONG are the filters calculated in the full wavelength range.
  ! FILTERS are the filter profiles in the narrow wavelength range, where we do the forward modeling.
  ! We do this to avoid computing the forward model very far into the continuum. It takes too much time. So we do
  ! a quick hack to compute the forward model only in an inner wavelength region and then add the filter contribution 
  ! outside of this region, in the continuum around the spectral line.
  !
  ! By RCE, May 17, 2011: Normalization of filters is now done in the wrapper. 
  ! This way, when we change parameters such as the wavelength range
  ! or the wavelength sampling, the normalization factor is computed for the correct parameters. 
  !
  ! By RCE May 24 2011: Defined the CHANGEVAR_FLAG, that decides whether the variable changes (defined in change_var.f90)
  ! for the inversion become in effect or not. There are IF statements in invert.f90 and the
  ! synthesis that check the flag before doing variable conversion.
  !
  ! By RCE, July 2011: I added a module to the code that performs variable changes/combinations of the model atmosphere.
  ! Instead of inverting for eta0 and DoppWidth, now we invert for eta0 and sqrt(eta0)*DoppWidth.
  ! The routines in the CHANGE_VAR module are called from invert at two different points. First after computing the derivatives. 
  ! The variable change is done for the derivatives. Then the gradient and Hessian are computed in order to obtain the perturbation
  ! to the model. This perturbation is given in the new variables, so we call a routine that undoes the change, so that we can
  ! add the perturbation to the model in order to move on to the next iteration.
  ! CHANGEVAR_FLAG is a flag that allows us to decide whether we want the inversion code to work with the changes or not.
  !
  ! By RCE Jan 2012:
  ! Introducing CHI2 REGULARIZATION: There is a regularization term in chi2. It penalizes chi2 for high values of eta0. We
  ! are using this to get rid of the double minima problem, by which, chi2 showed two minima at very different values of eta0
  ! (one low and one high), that was compensated by other parameters like the field strength and the doppler width.
  ! REGUL_FLAG switches the regularization on (1) and off (0).
  !
  ! By RCE March 2012:
  ! Jesper Schou has changed the convergence criterion and the requirements for random resets of the guess model for the 
  ! inversion of a given pixel. 
  ! Convergence criterion is now based on the LM lambda parameter. In general, if lambda reaches its maximum allowed value
  ! then the inversion will be considered converged. Pixels with high field strength (B>500G) will be forced to go to the maximum
  ! number of iterations, regardless of lambda. However, every time that lambda reaches lambda_max, the high field pixel will
  ! be given a new random initial guess.
  ! There are a number of cases, picked out empirically, that force a re-initialization of the inversion for a given pixel. 
  ! They are mostly related with low inclination magnetic fields (mostly vertical), especially in quiet Sun. 
  ! NOTE: The parametrization of the decision to re-initialize an inversion is very specific to the hmi.S_720s series. If the data 
  ! processing changes (especially if the noise changes), all of the parameters related to the resets will have to be tweaked.

  USE FILT_PARAM
  USE CONS_PARAM
  USE LINE_PARAM
  USE INV_PARAM
  USE INV_UTILS
  USE FORWARD
  USE RAN_MOD
  USE CHANGE_VAR
  IMPLICIT NONE
  !---------------------------------------------------------------------
  REAL(DP), INTENT(IN),  DIMENSION(NBINS*4)          :: OBS_LONG
  REAL(DP), INTENT(INOUT),  DIMENSION(NBINS*4)       :: SCAT_LONG   !ABGM
  REAL(DP),              DIMENSION(NUMW_LONG, NBINS) :: FILTERS_LONG
  REAL(DP), INTENT(INOUT),  DIMENSION(10)            :: GUESS
  REAL(DP),              DIMENSION(NUMW,NBINS)       :: FILTERS
  REAL(DP), INTENT(OUT), DIMENSION(10)               :: RES
  REAL(DP), INTENT(OUT), DIMENSION(12)               :: ERR

  !---------------------------------------------------------------------
  REAL(DP),        DIMENSION(NBINS,4)   :: OBS, SCAT, OBS_REV
  REAL(DP),        DIMENSION(NBINS)     :: INTEG_FILTERS
  REAL(DP),        DIMENSION(4)         :: WEIGHTS
  REAL(DP),        DIMENSION(10)        :: BESTMODEL, MODELG, LASTGOODMODEL, DMODEL, BESTMINMODEL
  REAL(DP),        DIMENSION(16)        :: SIGMA
  REAL(DP)                              :: BESTCHI2, LASTGOODCHI2, TOTPOL, NEWCHI2, LASTGOODLAMBDA, DELTACHI2, LAMBDAX, BESTMINCHI2
  REAL(DP),        DIMENSION(ITER)      :: LAMBDA, CHI2
  REAL(DP),        ALLOCATABLE          :: HESS(:,:), DIVC(:)
  REAL(DP),        DIMENSION(10,ITER)   :: MODEL
  REAL(DP)                              :: SSUM,SRATIO,BNOISE, B_LOW, B_HIGH, LAMBDA_START, LAMBDA_RESET
  INTEGER                               :: I, J, K, M, NPOINTS
  INTEGER                               :: CONV_FLAG, CONVERGENCE_FLAG, NRESET, CHANGEVAR_FLAG, REGUL_FLAG, RESET_FLAG, DONE, ITSTART, ABANDON, N_ABANDON
  !---------------------------------------------------------------------
  REAL(DP),     DIMENSION(NBINS,4)      :: SYN, LASTGOODSYN, BESTSYN
  REAL(DP),     DIMENSION(10,NBINS,4)   :: DSYN, LASTGOODDSYN
  !---------------------------------------------------------------------
  REAL(DP)                              :: POL         !!ABGM
  LOGICAL                               :: CALL_SYN2   !!ABGM
  LOGICAL                               :: FREE_INI_FF !!ABGM
  REAL(DP)                              :: POL_LOW     !!ABGM
  REAL(DP)                              :: POL_UPP     !!ABGM
  !---------------------------------------------------------------------
  CHARACTER(LEN=20), PARAMETER :: FMT = '("6f14.10")'
! Some variables for random number initialization
  integer                               :: dt_return(1:8), nseed,clock_count
  integer, dimension(:), allocatable    :: seed


  ! Allocate the Hessian and Divergence matrices
  ALLOCATE(HESS(NUMFREE_PARAM,NUMFREE_PARAM),DIVC(NUMFREE_PARAM))

  ! Re-dimension the observations and scattered ligth 
  CALL SORT_OBS(OBS_LONG,OBS)
  CALL SORT_OBS(SCAT_LONG,SCAT)

  FREE_INI_FF = FREE(10)  !ABGM

! ABGM
!!POL_LOW = 0.0025   !!polarization value used in Ana and Yang papers -> 0.25%
  IF (FREE(10).EQ..FALSE.) THEN
    CALL_SYN2 = .FALSE.
  ELSE
    POL = SUM(SQRT(SUM(OBS(:,2:4)**2,dim=2)),dim=1) / SUM(OBS(:,1),dim=1)
    IF ((POL.GT.POL_LOW).AND.(POL.LT.POL_UPP)) THEN
      CALL_SYN2 = .TRUE.
      GUESS(10) = 0.5
      FREE(10) = .TRUE.
    ELSE
      CALL_SYN2 = .FALSE.
      GUESS(10) = 1.0
      FREE(10) = .FALSE.
    ENDIF
  ENDIF
! ABGM


  ! --- BEGIN THE WAVELENGTH HACK
  ! By RCE April 2011: Select the region of the spectrum of the filters that we are 
  ! going to use in the synthesis and the region that we are going to integrate 
  ! assuming that they're in the continuum of the spectral line. 
  
  ! Difference between NUMW wavelength points for synthesis and NUMW_LONG points in 
  ! which filters were computed.
  ! We divide by 2 to get the number of wavelength points at each side of the
  ! "synthesis range"

  NPOINTS = (NUMW_LONG - NUMW)/2  

  ! Select region of filters corresponding to forward modeling wavelength region
  FILTERS(:,:) = FILTERS_LONG(NPOINTS+1:NUMW+NPOINTS,:) 
  ! Integrate the filters in the remaining regions at each side of the forward
  ! modeling region.
  DO I = 1, NBINS 
     INTEG_FILTERS(I) = SUM(FILTERS_LONG(1:NPOINTS,I)) +  &
          SUM(FILTERS_LONG(NUMW_LONG-NPOINTS+1:,I))
  ENDDO

  ! Making sure we add no filter integral if we're doing the forward modeling in the
  ! full wavelength range.
  IF (NPOINTS .EQ. 0) THEN 
     INTEG_FILTERS(:) = 0D0
  ENDIF
  ! --- END OF WAVELENGTH HACK

  !By RCE, Apr 23, 2010: Reversing the wavelength order of the observations 
  ! so that they are in increasing wavelength order

  DO k=1,NBINS
     DO m = 1,4
        OBS_REV(k,m) = OBS(NBINS-k+1,m)
     ENDDO
  ENDDO
  OBS(:,:) = OBS_REV(:,:)

  ! ----- FLAGS FOR THE CODE
  ! By RCE: 
  ! CONVERGENCE_FLAG is sent out to the wrapper with different values depending 
  !   on whether the algorithm converges or not and why.
  ! CHANGEVAR_FLAG decides whether variable change for inversion is implemented (1) or not (0).
  
  CONVERGENCE_FLAG = 0
  CHANGEVAR_FLAG = 1 ! Use new variables or not
  REGUL_FLAG = 1
  ! ---- END FLAGS

  !--- Iteration and re-initialization control
  ! The following parameters control the iterations, the convergence criterion and 
  ! if and when a pixel is re-setted (inversion re-initialized with different guess model)

  BNOISE=300.0D0 ! Fields above this are dominated by signal
  N_ABANDON=5 ! Abandon 1st reset if no better than this many iterations.
  DELTACHIMIN=1.0D-4 ! Min change of chi2 to decrease lambda for.
  LAMBDA_START=0.1D0 ! Initial value of lambda
  LAMBDA_RESET=0.1D0 ! Value of lambda after a reset
  LAMBDA_DOWN=5.0D0 ! How much to decrease lambda for good guess
  LAMBDA_UP=10.0D0 ! How much to increase lambda for bad guess
  LAMBDA_MIN=1.0D-4 ! Min lambda value
  LAMBDA_MAX=100.0D0 ! Max lambda value (exit criterion)
  

  !--- Get seed for random generator
  ! Both values incement with time, so not entirely independent
  call date_and_time(values=dt_return)
  call system_clock(count=clock_count,count_rate=i,count_max=j)
  
  call random_seed(size=nseed)
  allocate(seed(1:nseed))
  seed(:) = dt_return(8) ! Milliseconds of system clock
  seed(1)=clock_count ! Something to do with time spent


  !---------------------------------------------------------------------------------
  ! This IF statement checks whether the intensity is large enough to invert the data
  ! Otherwise the pixel is ignored. This is used to avoid pixels off-the-limb
  !---------------------------------------------------------------------------------
  IF (ICONT .GT. TREIC) THEN

     CALL GET_TOTPOL(OBS,TOTPOL)
     CALL GET_WEIGHT(OBS,WEIGHTS)


     !----------------------------------------------------
     ! Commented by RCE: 
     ! WFA_GUESS overwrites initial guess for Doppler velocity 
     ! eta0, doppler width, azimuth, source function and gradient.
     !----------------------------------------------------

     CALL WFA_GUESS(OBS/ICONT,LANDA0,GUESS(7), GUESS(3), GUESS(1), GUESS(5), GUESS(8), GUESS(9))

     ! Initialization of non-free parameters
     ! ABGM CHANGED THIS LINE FOR ff VARIABLE VERSION 
     IF (FREE(4).EQ..FALSE.) GUESS(4)=0.5D0
     !IF (FREE(10).EQ..FALSE.) GUESS(10)=1D0
     ! ABGM


     CALL FINE_TUNE_MODEL(GUESS)

     MODEL(:,1)=GUESS
     BESTMODEL(:)=GUESS

     ! By RCE: Weird default values for lambda and chi2 so that I can find out
     ! when they are not being updated.
     LAMBDA(:) = 33.3333
     CHI2(:) = 33.3333333
     LASTGOODCHI2=1E24
     LAMBDA(1)=LAMBDA_START
     BESTCHI2=1E24 ! Best chi2 ever seen
     BESTMINCHI2=1E24 ! Best value actually converged
     BESTMINMODEL=1E24 ! Best model actually converged
     ITSTART=1 ! Start of current iteration

     I=1

     ! Counter for the number of random resets that the model parameters suffer.
     ! RESET_FLAG to flag when the inversion of a pixel needs to be resetted.
     ! DONE flags the exit criterion for the main loop.
     NRESET=0
     RESET_FLAG=0
     DONE=0

     !-----------------------------
     ! Start LOOP iteration
     !-----------------------------
     DO WHILE ((I.LT.ITER).and.(DONE.EQ.0))
        ! Restart with new guess?
        IF (RESET_FLAG.ne.0) THEN
           BESTMINCHI2=BESTCHI2 ! Save current best chi2
           BESTMINMODEL=BESTMODEL ! Save current best chi2
           ITSTART=I ! Remember where we started
           MODELG=BESTMODEL
           CALL RANDOM_MODEL_JUMP(MODELG)
           if (mod(nreset,2).eq.0) then ! Do clever resets every other time.
              if (bestmodel(6).ge.BNOISE) then
                 modelg=bestmodel
                 if (bestmodel(6).gt.1000) then
                    modelg(1)=50.0d0
                    modelg(5)=12.0d0
                 else
                    if (mod(nreset,4).eq.2) then
                       modelg(1)=10.0d0
                       modelg(5)=7.0d0
                       sratio=0.2d0
                    else
                       modelg(1)=9.0d0
                       modelg(5)=30.0d0
                       sratio=1.0d0
                    endif
                    ssum=bestmodel(8)+bestmodel(9)
                    modelg(8)=ssum*sratio/(1.0d0+sratio)
                    modelg(9)=ssum/(1.0d0+sratio)
                 endif
              endif
              if (BESTMODEL(6).lt.BNOISE) then
                 modelg=bestmodel
                 modelg(1)=8.0d0
                 modelg(5)=0.75d0*bestmodel(5)
                 modelg(6)=85.0d0
! Predict S0/S1 from Doppler width. Maintain S0+S1
                 ssum=bestmodel(8)+bestmodel(9)
                 sratio=0.017d0*bestmodel(5)
                 modelg(8)=ssum*sratio/(1.0d0+sratio)
                 modelg(9)=ssum/(1.0d0+sratio)
              endif
           endif

           CALL FINE_TUNE_MODEL(MODELG)
           MODEL(:,I)=MODELG
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=1.0d24
           LAMBDA(I)=LAMBDA_RESET
           NRESET = NRESET+1
        ENDIF

        ! ABGM added this condition for ff variable version 
        IF (CALL_SYN2 .EQV. .TRUE.) THEN
          CALL SYNTHESIS2(MODEL(:,I),SCAT,.FALSE.,SYN,DSYN, FILTERS, INTEG_FILTERS)
        ELSE 
          CALL SYNTHESIS(MODEL(:,I),SCAT,.FALSE.,SYN,DSYN, FILTERS, INTEG_FILTERS)
        ENDIF
        ! ABGM

        !--- Get chi2 only, see if it is better and only calculate derivatives if it is.
        CALL GET_CHI2(MODEL(:,I),SYN,OBS,WEIGHTS,REGUL_FLAG,NEWCHI2)
        ! Checking for NAN in NEWCHI2
        IF (NEWCHI2.EQ.NEWCHI2+1.0D0) THEN
           PRINT*,'NaN detected in Subroutine GETCHI2'
           NEWCHI2=2*LASTGOODCHI2
        ENDIF
        CHI2(I)=NEWCHI2
        DELTACHI2=LASTGOODCHI2-NEWCHI2  ! change in chi2

        IF (DELTACHI2.GE.0.0d0) THEN ! ---- Things got better
           ! Get synthetic profiles and derivatives

           ! ABGM added this condition for ff variable version 
           IF (CALL_SYN2 .EQV. .TRUE.) THEN
             CALL SYNTHESIS2(MODEL(:,I),SCAT,.TRUE.,SYN,DSYN, FILTERS, INTEG_FILTERS)
           ELSE
             CALL SYNTHESIS(MODEL(:,I),SCAT,.TRUE.,SYN,DSYN, FILTERS, INTEG_FILTERS)
           ENDIF
           ! ABGM
           !--------------------------------------------
           ! Calling change of variables for derivatives
           !--------------------------------------------
           IF (CHANGEVAR_FLAG .EQ. 1) THEN
              CALL DO_CHANGE_DER(MODEL(:,I), DSYN)    
           ENDIF
           ! From now on, the derivatives correspond to the changed variable ones.
           ! Normalizing Derivatives
           CALL NORMALIZE_DSYN(DSYN)
           ! Set to Zero unneeded derivatives
           CALL ZERO_DSYN(DSYN)
           
           LASTGOODLAMBDA=LAMBDA(I)
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=NEWCHI2
           LASTGOODSYN=SYN
           LASTGOODDSYN=DSYN

           IF (NEWCHI2.LT.BESTCHI2) THEN
              BESTMODEL=MODEL(:,I)
              BESTCHI2=NEWCHI2
              BESTSYN=SYN
           ENDIF

        ELSE !---- Things did not get better
           SYN=LASTGOODSYN
           DSYN=LASTGOODDSYN
           MODEL(:,I)=LASTGOODMODEL
        ENDIF

        ! Levenberg-Marquardt lambda factor
        CALL GET_LAMBDA(DELTACHI2,LAMBDA(I),LAMBDA(I+1))

        ! Getting Divergence and Hessian.
        CALL GET_DIVC(LASTGOODSYN,OBS,LASTGOODDSYN,WEIGHTS,DIVC)
        CALL GET_HESS(LASTGOODDSYN,WEIGHTS,HESS)

        ! Get perturbations to old model
        CALL GET_DMODEL(LASTGOODMODEL, REGUL_FLAG, DIVC, HESS, LAMBDA(I+1), DMODEL, CONV_FLAG)
        ! Ignore CONV_FLAG and rely on DMODEL set to zero/LAMBDA increase.

        ! Normalizing, trimming and tuning DMODEL
        ! Undoing variable change to obtain the perturbations to
        ! the original model parameters. Overwriting DMODEL!
        CALL NORMALIZE_DMODEL(DMODEL)
        IF (CHANGEVAR_FLAG .EQ. 1) THEN
           CALL UNDO_CHANGE_DMODEL_FINITE(MODEL(:,I),DMODEL)   
        ENDIF
        CALL CUT_DMODEL(DMODEL,MODEL)
        MODEL(:,I+1)=LASTGOODMODEL+DMODEL
        CALL FINE_TUNE_MODEL(MODEL(:,I+1))

        !---- Decision making:
        !     Have we converged?
        !     Do we need to try a new initial guess for this pixel?
        !     Do we just carry on with the iterations?

        DONE=0
        RESET_FLAG=0
        ABANDON=0
        if ((nreset.eq.1).and.((i-itstart).ge.N_ABANDON).and.((lastgoodchi2-bestminchi2).gt.0).and.(bestmodel(6).lt.BNOISE)) ABANDON=1 ! Restart is not getting better.
        if ((LAMBDA(I+1).GE.LAMBDA_MAX).or.(ABANDON.ne.0)) then ! We are done with this round.

           ! Test if reset is desired.
           IF (NRESET .EQ. 0) THEN ! Only one reset allowed for the following cases
              if (BESTMODEL(6).lt.30.0D0) RESET_FLAG=1 ! Low field
              if (BESTMODEL(6).gt.BNOISE) RESET_FLAG=1 ! High field
              if (BESTCHI2.gt.20.00d0) RESET_FLAG=1
              if ((BESTCHI2.gt.2.75d0).and.(BESTMODEL(6).lt.70.0d0).and.(BESTMODEL(5).lt.15.0d0)) RESET_FLAG=1 ! A few really bad points
              if (abs(BESTMODEL(2)-90.0d0).gt.(10*log(BESTMODEL(6)/7))) RESET_FLAG=1 ! A stab at high incl points
              if ((BESTMODEL(1).lt.5.0d0).and.(BESTMODEL(5).gt.20.0d0).and.(BESTMODEL(6).lt.BNOISE)) RESET_FLAG=2
           ENDIF
              
           B_LOW=35.0D0-nreset
           B_HIGH=500.0D0
           if (abs(BESTMODEL(2)-90.0d0).gt.60.0d0) RESET_FLAG=1 ! Inclination spike
           if (BESTMODEL(6).lt.B_LOW) RESET_FLAG=1 ! Low field
           if (BESTMODEL(6).gt.B_HIGH) RESET_FLAG=1 ! High field
           if ((BESTMODEL(1).lt.3.5d0).and.(BESTMODEL(5).gt.30.0d0).and.(BESTMODEL(6).lt.45.0d0)) RESET_FLAG=2

           if (RESET_FLAG.EQ.0) DONE=1 ! No resets desired
        endif

        I=I+1
     ENDDO ! Main iteration loop

     ! Get errors
     IF (BESTCHI2.gt.1d20) THEN
        CONVERGENCE_FLAG = 4 ! Set flag and give up
     ELSE
       
        ! We compute the derivatives with the best model parameters (without
        ! doing the change of variable afterwards), and we compute the Hessian
        ! and the errors from there.
        
        ! Recalculate chi2 and derivatives. No regularization.
        ! JS: Not clear why NEWCHI2 is recalculated.
        CALL GET_CHI2(BESTMODEL,BESTSYN,OBS,WEIGHTS,0,NEWCHI2)

        ! ABGM added this condition for ff variable version 
        IF (CALL_SYN2 .EQV. .TRUE.) THEN
          CALL SYNTHESIS2(BESTMODEL,SCAT,.TRUE.,SYN,LASTGOODDSYN, FILTERS, INTEG_FILTERS)
        ELSE
          CALL SYNTHESIS(BESTMODEL,SCAT,.TRUE.,SYN,LASTGOODDSYN, FILTERS, INTEG_FILTERS)
        ENDIF
        ! ABGM

        CALL NORMALIZE_DSYN(LASTGOODDSYN)

        ! Compute Hessian
        Call GET_HESS(LASTGOODDSYN,WEIGHTS,HESS)
        ! This Hessian is constructed with normalized derivatives
        CALL GET_ERR(HESS, BESTCHI2,SIGMA,CONV_FLAG)

        if (CONV_FLAG.NE.0) CONVERGENCE_FLAG = 5 ! ERROR in SVD: Set flag, otherwise proceed.

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

        ! Final Result
        RES=BESTMODEL
        
        ! By RCE, Apr 23, 2010: Juanma's correction to azimuth (which is
        ! rotated by 90 degrees with respect to the convention people want).
        RES(3) = RES(3) + 90D0
        IF (RES(3) .GT. 180D0) RES(3) = RES(3) - 180D0

        IF ((I .EQ. ITER).and.(NRESET.eq.0)) THEN 
           CONVERGENCE_FLAG = 2 ! Not converged
        ELSE
           CONVERGENCE_FLAG = 0 ! Found at least one local minimum.
        ENDIF
     ENDIF ! Decent chisq

  ELSE ! Intensity too low to invert
     CONVERGENCE_FLAG = 1
  ENDIF

! CONVERGENCE_FLAG values
!    0: Converged!
!    1: Pixel not inverted due to ICONT LT threshold intensity.
!    2: Maximum number of iterations reached without convergence.
!    4: Never got decent chisq.
!    5: Could not get errors.

  DEALLOCATE(HESS, DIVC)

! ABGM
DO I=1,4 
  SCAT_LONG(NBINS*(I-1)+1:NBINS*I) = SCAT(:,I)
ENDDO

FREE(10) = FREE_INI_FF
! ABGM

END SUBROUTINE INVERT


! subroutines added by K.H.
! consolidate (some of) allocate arrays 
SUBROUTINE VFISVALLOC(NUM_LAMBDA_FILTER, NUM_LAMBDA_LONG, NUM_LAMBDA)
  USE FILT_PARAM
  USE LINE_PARAM
  USE CONS_PARAM
  USE INV_PARAM
  IMPLICIT NONE
  INTEGER,  INTENT(IN)           :: NUM_LAMBDA_FILTER, NUM_LAMBDA_LONG, NUM_LAMBDA
  
  NUMW_LONG = NUM_LAMBDA_LONG
  NBINS = NUM_LAMBDA_FILTER

  ALLOCATE (FILTER(NUMW,NBINS),TUNEPOS(NBINS))

  NUMW=NUM_LAMBDA
  ALLOCATE (WAVE(NUMW))
     
END SUBROUTINE VFISVALLOC

!CVSVERSIONINFO "$Id: invert_2comp.f90,v 1.1 2022/03/04 17:31:38 arta Exp $"
