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
  ! By RCE, Apr 23, 2010: (OBSOLETE, see May 17 2011 comment) Normalize all filters to central filter area hard-coded value:
  ! The integrated areas of the six filters for the central pixel (pixel # 8386559) are:    
  ! 1.43858806880917  1.49139978559615   1.52324167542270  1.52811487224149  1.50984871281028  1.46553486521323
  ! the average area being: 1.49279. We will normalize ALL the filters by this value.
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
  !  - Put lower bound on the LM lambda parameter so that it doesn't decrease beyond a value (in GET_LAMBDA in INV_UTILS) 
  !  - Test flag (TMP_FLAG) to count the number of resets of the model parameters. Saved in one of the error variables
  !
  ! By RCE: Re-normalizing filters to different value to account for +/-2A coverage (OBSOLETE, see May 17 2011 comment)
  !
  ! By RCE April 2011: Included filter hack. FILTERS_LONG are the filters calculated in the full wavelength range.
  ! FILTERS are the filter profiles in the narrow wavelength range, where we do the forward modeling.
  ! We do this to avoid computing the forward model very far into the continuum. It takes too much time. So we do
  ! a quick hack to compute the forward model only in an inner wavelength region and then add the filter contribution 
  ! outside of this region, in the continuum around the spectral line.
  !
  ! By RCE, May 17, 2011: Normalization of filters is now done in the wrapper. This way, when we change parameters such as the wavelength range
  ! or the wavelength sampling, the normalization factor is computed for the correct parameters. 
  ! By RCE May 24 2011: Defined the CHANGEVAR_FLAG, that decides whether the variable changes
  ! for the inversion become in effect or not. There are IF statements in invert.f90 and the
  ! synthesis that check the flag before doing variable conversion.
  !
  ! By RCE, July 2011: I added a module to the code that performs variable changes/combinations of the model atmosphere.
  ! Instead of inverting for eta0 and DoppWidth, now we invert for eta0 and sqrt(eta0)*DoppWidth.
  ! Instead of inverting B0 and B1 (source function and gradient) we invert B0 and B0+B1.
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
  REAL(DP),        DIMENSION(10)        :: BESTMODEL, MODELG, LASTGOODMODEL, DMODEL,test1,test2,GUESS_TWEAK,BESTMINMODEL
  REAL(DP),        DIMENSION(16)        :: SIGMA
  REAL(DP)                              :: BESTCHI2, LASTGOODCHI2, TOTPOL, NEWCHI2, LASTGOODLAMBDA, ICONT, DELTACHI2, DUMMY, LAMBDAX, BESTMINCHI2
  REAL(DP),        DIMENSION(ITER)      :: LAMBDA, CHI2
  REAL(DP),        ALLOCATABLE          :: HESS(:,:), DIVC(:)
  REAL(DP),        DIMENSION(10,ITER)   :: MODEL
  REAL(DP)                              :: SSUM,SRATIO,BNOISE, divc_norm, B_LOW, B_HIGH, LAMBDA_START, LAMBDA_RESET
  INTEGER                               :: I, J, K, M, NPOINTS
  INTEGER                               :: CONV_FLAG, CONVERGENCE_FLAG, NRESET, CHANGEVAR_FLAG, REGUL_FLAG, RESET_FLAG, DONE, ITSTART, ABANDON, N_ABANDON
  REAL(DP),        DIMENSION(NBINS)     :: pweights
  !---------------------------------------------------------------------
  REAL(DP),     DIMENSION(NBINS,4)      :: SYN, LASTGOODSYN, BESTSYN
  REAL(DP),     DIMENSION(10,NBINS,4)   :: DSYN, LASTGOODDSYN
  !---------------------------------------------------------------------
  CHARACTER(LEN=20), PARAMETER :: FMT = '("6f14.10")'
! Some variables for random number initialization
  integer                               :: dt_return(1:8), nseed,clock_count
  integer, dimension(:), allocatable :: seed


  CALL SORT_OBS(OBS_LONG,OBS)
  CALL SORT_OBS(SCAT_LONG,SCAT)

  ALLOCATE(HESS(NUMFREE_PARAM,NUMFREE_PARAM),DIVC(NUMFREE_PARAM))


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

! Some constants for later for control of iterations etc.
  BNOISE=300D0 ! Fields above this are dominated by signal
  N_ABANDON=5 ! Abandon 1st reset if no better than this many iterations.
  DELTACHIMIN=1D-4 ! Min change of chi2 to decrease lambda for.
  LAMBDA_START=0.1D0 ! Initial value of lambda
  LAMBDA_RESET=0.1D0 ! Value of lambda after a reset
  LAMBDA_DOWN=5D0 ! How much to decrease lambda for good guess
  LAMBDA_UP=10D0 ! How much to increase lambda for bad guess
  LAMBDA_MIN=1D-4 ! Min lambda value
  LAMBDA_MAX=100D0 ! Max lambda value (exit criterion)
  
! Get seed for random generator
! Both values incement with time, so not entirely independent
  call date_and_time(values=dt_return)
  call system_clock(count=clock_count,count_rate=i,count_max=j)
  
  call random_seed(size=nseed)
  allocate(seed(1:nseed))
  seed(:) = dt_return(8) ! Milliseconds of system clock
  seed(1)=clock_count ! Something to do with time spent

  ICONT=MAXVAL(OBS(:,1))
  IF (ICONT .GT. TREIC) THEN
     CALL INV1_INIT (ICONT,REGUL_FLAG) ! Initialize limits and norms
     !---------------------------------------------------------------------------------
     ! This IF statement checks whether the intensity is large enough to invert the data
     ! Otherwise the pixel is ignored. This is used to avoid pixels off-the-limb
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
     GUESS_TWEAK=GUESS
     GUESS_TWEAK(1)=5
     GUESS_TWEAK(3)=atan2(sum(obs(:,3)),sum(obs(:,2)))*90/DPI
! The following is slightly better
!    pweights=[0.4,1.0,1.0,1.0,1.0,0.4]
!    GUESS_TWEAK(3)=atan2(sum(pweights*obs(:,3)),sum(pweights*obs(:,2)))*90/DPI
     GUESS_TWEAK(5)=20 ! More reasonable Doppler width
     GUESS_TWEAK(8)=0.15D0*0.98D0*ICONT
     GUESS_TWEAK(9)=0.85D0*0.98D0*ICONT
!     call random_seed(put=seed)
!     GUESS_TWEAK(1)=exp(5*ran1()) ! eta0
!     GUESS_TWEAK(2)=180*ran1() ! incl
!     GUESS_TWEAK(3)=180*ran1() ! az
!     GUESS_TWEAK(5)=5+50*ran1() ! Dop
!     GUESS_TWEAK(6)=exp(2.5+5*ran1()) ! B
!     GUESS_TWEAK(7)=GUESS_TWEAK(7)+(ran1()-0.5)*100000 ! vlos +/- 500m/s
!     GUESS_TWEAK(8)=(0.15+0.90*ran1())*ICONT
!     GUESS_TWEAK(9)=ICONT-GUESS_TWEAK(8)

     IF (FREE(4).EQ..FALSE.) GUESS_TWEAK(4)=0.5D0
     IF (FREE(10).EQ..FALSE.) GUESS_TWEAK(10)=1D0

     CALL FINE_TUNE_MODEL(GUESS_TWEAK)

     MODEL(:,1)=GUESS_TWEAK
     BESTMODEL(:)=GUESS_TWEAK

     ! By RCE: Default values for lambda and chi2 so that I can find out
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
     NRESET=0
     DUMMY=0
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
           CALL RANDOM_MODEL_JUMP(RANDOM_JUMP,MODELG)
           if (mod(nreset,2).eq.0) then ! Do clever resets every other time.
              if (bestmodel(6).ge.BNOISE) then
                 modelg=bestmodel
                 if (bestmodel(6).gt.1000) then
                    modelg(1)=50
                    modelg(5)=12
                 else
                    if (mod(nreset,4).eq.2) then
                       modelg(1)=10
                       modelg(5)=7
                       sratio=.2
                    else
                       modelg(1)=9
                       modelg(5)=30
                       sratio=1
                    endif
                    ssum=bestmodel(8)+bestmodel(9)
                    modelg(8)=ssum*sratio/(1+sratio)
                    modelg(9)=ssum/(1+sratio)
                 endif
              endif
              if (BESTMODEL(6).lt.BNOISE) then
                 modelg=bestmodel
                 modelg(1)=8
                 modelg(5)=0.75*bestmodel(5)
                 modelg(6)=85
! Predict S0/S1 from Doppler width. Maintain S0+S1
                 ssum=bestmodel(8)+bestmodel(9)
                 sratio=0.017*bestmodel(5)
                 modelg(8)=ssum*sratio/(1+sratio)
                 modelg(9)=ssum/(1+sratio)
              endif
           endif
           CALL FINE_TUNE_MODEL(MODELG)
           MODEL(:,I)=MODELG
           LASTGOODMODEL=MODEL(:,I)
           LASTGOODCHI2=1E24
           LAMBDA(I)=LAMBDA_RESET
           NRESET = NRESET+1
        ENDIF

        ! Get chi2 only, see if it is better and only calculate derivatives if it is.
        CALL SYNTHESIS(MODEL(:,I),SCAT,.FALSE.,SYN,DSYN, FILTERS, INTEG_FILTERS)
        CALL GET_CHI2(MODEL(:,I),SYN,OBS,WEIGHTS,REGUL_FLAG,NEWCHI2)
        ! Checking for NAN in NEWCHI2
        IF (NEWCHI2.EQ.NEWCHI2+1D0) THEN
           PRINT*,'NaN detected in Subroutine GETCHI2'
           NEWCHI2=2*LASTGOODCHI2
        ENDIF
        CHI2(I)=NEWCHI2
        DELTACHI2=LASTGOODCHI2-NEWCHI2
        IF (DELTACHI2.GE.0.0) THEN ! Things got better
           ! Get synthetic profiles and derivatives
           CALL SYNTHESIS(MODEL(:,I),SCAT,.TRUE.,SYN,DSYN, FILTERS, INTEG_FILTERS)
           !--------------------------------------------
           ! Calling change of variables for derivatives
           !--------------------------------------------
           IF (CHANGEVAR_FLAG .EQ. 1) THEN
              CALL DO_CHANGE_DER(MODEL(:,I), DSYN)     !!NEW!!
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
        ELSE ! Things did not get better
           SYN=LASTGOODSYN
           DSYN=LASTGOODDSYN
           MODEL(:,I)=LASTGOODMODEL
        ENDIF
        CALL GET_LAMBDA(DELTACHI2,LAMBDA(I),LAMBDA(I+1))

        ! Getting Divergence and Hessian.
        CALL GET_DIVC(LASTGOODSYN,OBS,LASTGOODDSYN,WEIGHTS,DIVC)
        divc_norm=sqrt(sum(divc*divc))
        CALL GET_HESS(LASTGOODDSYN,WEIGHTS,HESS)

        ! Get perturbations to old model
        CALL GET_DMODEL(LASTGOODMODEL, REGUL_FLAG, DIVC, HESS, LAMBDA(I+1), DMODEL, CONV_FLAG)
        ! Ignore CONV_FLAG and rely on DMODEL set to zero/LAMBDA increase.

        ! Normalizing, trimming and tuning DMODEL
        ! Undoing variable change to obtain the perturbations to
        ! the original model parameters. Overwriting DMODEL!
        CALL NORMALIZE_DMODEL(DMODEL)
        IF (CHANGEVAR_FLAG .EQ. 1) THEN
           CALL UNDO_CHANGE_DMODEL_FINITE(MODEL(:,I),DMODEL)    !!NEW!!
        ENDIF
        CALL CUT_DMODEL(DMODEL,MODEL)
        MODEL(:,I+1)=LASTGOODMODEL+DMODEL
        CALL FINE_TUNE_MODEL(MODEL(:,I+1))

        DONE=0
        RESET_FLAG=0
        ABANDON=0
        if ((nreset.eq.1).and.((i-itstart).ge.N_ABANDON).and.((lastgoodchi2-bestminchi2).gt.0).and.(bestmodel(6).lt.BNOISE)) ABANDON=1 ! Restart is not getting better.
        if ((LAMBDA(I+1).GE.LAMBDA_MAX).or.(ABANDON.ne.0)) then ! We are done with this round.

           ! Test if reset is desired.
           if (abs(BESTMODEL(2)-90).gt.45) RESET_FLAG=1 ! Inclination spike
           if (BESTMODEL(6).lt.30) RESET_FLAG=1 ! Low field
           if (BESTMODEL(6).gt.BNOISE) RESET_FLAG=1 ! High field
           if (BESTCHI2.gt.20.00) RESET_FLAG=1
           if ((BESTCHI2.gt.2.75).and.(BESTMODEL(6).lt.70).and.(BESTMODEL(5).lt.15)) RESET_FLAG=1 ! A few really bad points
           if (abs(BESTMODEL(2)-90).gt.(10*log(BESTMODEL(6)/7))) RESET_FLAG=1 ! A stab at high incl points
           if ((BESTMODEL(1).lt.5).and.(BESTMODEL(5).gt.20).and.(BESTMODEL(6).gt.000).and.(BESTMODEL(6).lt.BNOISE)) RESET_FLAG=2

           if (NRESET.ne.0) RESET_FLAG=0 ! Force a max of one reset for items above.

           B_LOW=35D0-nreset
           B_HIGH=500D0
           if (abs(BESTMODEL(2)-90).gt.60) RESET_FLAG=1 ! Inclination spike
           if (BESTMODEL(6).lt.B_LOW) RESET_FLAG=1 ! Low field
           if (BESTMODEL(6).gt.B_HIGH) RESET_FLAG=1 ! High field
           if ((BESTMODEL(1).lt.3.5).and.(BESTMODEL(5).gt.30).and.(BESTMODEL(6).lt.45)) RESET_FLAG=2

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
        CALL SYNTHESIS(BESTMODEL,SCAT,.TRUE.,SYN,LASTGOODDSYN, FILTERS, INTEG_FILTERS)

        CALL NORMALIZE_DSYN(LASTGOODDSYN)

        ! Compute Hessian
        Call GET_HESS(LASTGOODDSYN,WEIGHTS,HESS)
        ! This Hessian is constructed with normalized derivatives
        CALL GET_ERR(HESS, BESTCHI2,SIGMA,CONV_FLAG)
        if (CONV_FLAG.NE.0) CONVERGENCE_FLAG = 5 ! Set flag, otherwise proceed.

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
        ! Replace some errors with assorted other parameters
        !ERR(9)=Nreset ! Save number of resets
        !ERR(10)=I! Save number of iterations
        !ERR(11)=ICONT

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

END SUBROUTINE INVERT


! subroutines added by K.H.
! consolidate (some of) allocate arrays 
SUBROUTINE VFISVALLOC(NUM_LAMBDA_FILTER, NUM_TUNNING, NUM_LAMBDA_LONG, NUM_LAMBDA)
  USE FILT_PARAM
  USE LINE_PARAM
  USE CONS_PARAM
  USE INV_PARAM
  IMPLICIT NONE
  INTEGER,  INTENT(IN)           :: NUM_LAMBDA_FILTER, NUM_TUNNING, NUM_LAMBDA_LONG, NUM_LAMBDA
  
  NUMW_LONG = NUM_LAMBDA_LONG
  NBINS = NUM_LAMBDA_FILTER
  NTUNE = NUM_TUNNING
  ALLOCATE (FILTER(NUMW,NBINS),TUNEPOS(NBINS))

  NUMW=NUM_LAMBDA
  ALLOCATE (WAVE(NUMW))
     
END SUBROUTINE VFISVALLOC

!CVSVERSIONINFO "$Id: invert.f90,v 1.7 2012/04/09 22:21:17 keiji Exp $"
