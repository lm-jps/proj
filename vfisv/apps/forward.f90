MODULE FORWARD
  !
  ! J M Borrero
  ! Jan 10, 2007
  ! HAO-NCAR for HMI-Stanford
  !
  ! By RCE: Added FILTERS as a parameter passed to SYNTHESIS
  ! Defined FILTERS in SYNTHESIS
  ! Changed all the FILTER (Juanma's filter profiles) instances 
  ! by FILTERS (Sebastien's filter profiles). Filters is now defined
  ! as a matrix of [NUMW, NBINS] elements. We've exchanged the 
  ! indices with respect to the C-wrapper definition because C and Fortran
  ! read the elements in different order.
  ! By RCE, Feb 2011: Implemented changes in the derivatives that J.M Borrero
  ! suggested after changing the voigt.f function by a factor of 2. It affects
  ! the derivatives of PHI and PSI with respect to the Damping and FREC(B,R,P)
  !
  ! By RCE, April 2011: Adding the integral of the filter profiles times the continuum
  ! (S0+S1) to Stokes I, for the outer wavelength range where the forward modeling is not done.

CONTAINS
  !!
  !! SUBROUTINE SYNTHESIS
  !!
  SUBROUTINE SYNTHESIS(MODEL,SCAT,DERIVATIVE,SYN,DSYN, FILTERS, INTEG_FILTERS)
    USE FILT_PARAM
    USE LINE_PARAM
    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP), INTENT(IN),  DIMENSION(10)          :: MODEL
    REAL(DP), INTENT(IN),  DIMENSION(NBINS,4)     :: SCAT
    REAL(DP), INTENT(IN),  DIMENSION(NUMW, NBINS)  :: FILTERS
    REAL(DP), INTENT(IN),  DIMENSION(NBINS)       :: INTEG_FILTERS
    LOGICAL,  INTENT(IN)                          :: DERIVATIVE
    REAL(DP), INTENT(OUT),  DIMENSION(NBINS,4)    :: SYN
    REAL(DP), INTENT(OUT),  DIMENSION(10,NBINS,4) :: DSYN
    REAL(DP),               DIMENSION(NBINS,4)    :: SYN_MAG
    REAL(DP),               DIMENSION(9,NBINS,4)  :: DSYN_MAG
    !------------------------------------------------------
    REAL(DP),    DIMENSION(NUMW)       :: ETAI, ETAQ, ETAU, ETAV, RHOQ, RHOU, RHOV
    REAL(DP),    DIMENSION(7,NUMW)     :: DerETAI, DerETAQ, DerETAU, DerETAV
    REAL(DP),    DIMENSION(7,NUMW)     :: DerRHOQ, DerRHOU, DerRHOV
    REAL(DP),    DIMENSION(NUMW)       :: EXTRA, DET_MAT
    REAL(DP),    DIMENSION(9,NUMW,4)   :: DSTOKES_MAG
    REAL(DP),    DIMENSION(NUMW,4)     :: STOKES_MAG
    REAL(DP),    DIMENSION(7,NUMW)     :: DEXTRA, DDMAT
    REAL(DP)                           :: S0, S1, ALPHAM
    REAL(DP),    DIMENSION(NUMW)       :: A1, A2, A3, A4, A5, A6, A7
    REAL(DP),    DIMENSION(NUMW)       :: B1, B2, B3, B4, B5, B6, B7
    REAL(DP),    DIMENSION(NUMW)       :: C1, C2, C3, C4, C5, C6
    REAL(DP),    DIMENSION(NUMW)       :: D1, D2, D3, D4, D5, D6
    REAL(DP),    DIMENSION(NUMW)       :: PART1, PART2
    INTEGER                            :: I, J, K, M
    !------------------------------------------------------
    S0=MODEL(8)
    S1=MODEL(9)
    ALPHAM=MODEL(10)
    !
    SYN(:,:)=0D0
    DSYN(:,:,:)=0D0
    SYN_MAG(:,:)=0D0
    DSYN_MAG(:,:,:)=0D0
    ETAI(:)=0D0
    ETAQ(:)=0D0
    ETAU(:)=0D0
    ETAV(:)=0D0
    RHOQ(:)=0D0
    RHOU(:)=0D0
    RHOV(:)=0D0
    DerETAI(:,:)=0D0
    DerETAQ(:,:)=0D0
    DerETAU(:,:)=0D0
    DerETAV(:,:)=0D0
    DerRHOQ(:,:)=0D0
    DerRHOU(:,:)=0D0
    DerRHOV(:,:)=0D0
    STOKES_MAG(:,:)=0D0
    EXTRA(:)=0D0
    DET_MAT(:)=0D0
    DSTOKES_MAG(:,:,:)=0D0
    DEXTRA(:,:)=0D0
    DDMAT(:,:)=0D0

    CALL ABSMAT(MODEL, DERIVATIVE, ETAI, ETAQ, ETAU, ETAV, RHOQ, RHOU, RHOV, &
         DerETAI, DerETAQ, DerETAU, DerETAV, DerRHOQ, DerRHOU, DerRHOV)
    ! Common parts
    EXTRA = ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV
    DET_MAT = (ETAI**2D0)*(ETAI**2D0-ETAQ**2D0-ETAU**2D0-ETAV**2D0+ &
         RHOQ**2D0+RHOU**2D0+RHOV**2D0)-EXTRA**2D0
    !----------------------------------------------------------------------------
    ! Solution to the Unno-Rachkovski equations
    ! This is the Stokes vector coming from the magnetic atmosphere
    !----------------------------------------------------------------------------
    STOKES_MAG(:,1) = S0+(1D0/DET_MAT)*ETAI*(ETAI**2D0+RHOQ**2D0+RHOU**2D0+RHOV**2D0)*S1
    STOKES_MAG(:,2) = -(1D0/DET_MAT)*(ETAI**2D0*ETAQ+ETAI*(ETAV*RHOU-ETAU*RHOV)+RHOQ*EXTRA)*S1
    STOKES_MAG(:,3) = -(1D0/DET_MAT)*(ETAI**2D0*ETAU+ETAI*(ETAQ*RHOV-ETAV*RHOQ)+RHOU*EXTRA)*S1
    STOKES_MAG(:,4) = -(1D0/DET_MAT)*(ETAI**2D0*ETAV+ETAI*(ETAU*RHOQ-ETAQ*RHOU)+RHOV*EXTRA)*S1
    !-----------------------------------------------------------------------------
    ! Now we apply HMI Filters
    !-----------------------------------------------------------------------------
    DO K=1,4
       DO J=1,NBINS
          SYN_MAG(J,K)=SUM(FILTERS(:,J)*STOKES_MAG(:,K))     
       ENDDO
    ENDDO

    ! By RCE, April 2011: Adding integral of filters outside wavelength range for Stokes I
    ! We're assuming that this outer wavelength range corresponds to continuum, hence we 
    ! multiply the integral of the filters by the continuum for Stokes I and assume it is 0
    ! for Stokes Q, U, and V
    DO J =1, NBINS
      SYN_MAG(J, 1) = SYN_MAG(J,1)+INTEG_FILTERS(J)*(S0+S1)
    ENDDO

    !---------------------------------------------------------
    ! Total Stokes vector including the non-magnetic component
    !---------------------------------------------------------
    SYN=(1D0-ALPHAM)*SCAT+ALPHAM*SYN_MAG
    !---------------------------------------------------------
    ! Derivatives
    !---------------------------------------------------------
    IF (DERIVATIVE.EQ..TRUE.) THEN
       ! Derivatives of the Stokes Parameters (emering from Magnetic component)
       ! with respect to the 7 free parameters: eta0, gam, phi, dam, dldop, B, Vlos
       ! plus two new dependences: S0, S1. Total 9 free parameters.
       ! First derivatives of EXTRA and DET_MAT
       ! These ones do not depend on S0, S1
       DO I=1,7
          A1 = ETAI
          A2 = ETAQ
          A3 = ETAU
          A4 = ETAV
          A5 = RHOQ
          A6 = RHOU
          A7 = RHOV
          B1 = DerETAI(I,:)
          B2 = DerETAQ(I,:)
          B3 = DerETAU(I,:)
          B4 = DerETAV(I,:)
          B5 = DerRHOQ(I,:)
          B6 = DerRHOU(I,:)
          B7 = DerRHOV(I,:)
          DEXTRA(I,:)=A2*B5+A5*B2+A3*B6+A6*B3+A4*B7+A7*B4
          DDMAT(I,:)=2D0*A1*B1*(A1**2D0-A2**2D0-A3**2D0-A4**2D0+A5**2D0+ &
               A6**2D0+A7**2D0)+2D0*A1**2D0*(A1*B1-A2*B2-A3*B3-A4*B4+A5*B5+ &
               A6*B6+A7*B7)-2D0*EXTRA*DEXTRA(I,:)
       ENDDO
       ! Now derivatives of Stokes I with respect to
       ! 7 regular free parameters
       DO I=1,7
          A1 = ETAI
          B1 = DerETAI(I,:)
          A2 = RHOQ
          B2 = DerRHOQ(I,:)
          A3 = RHOU
          B3 = DerRHOU(I,:)
          A4 = RHOV
          B4 = DerRHOV(I,:)
          DSTOKES_MAG(I,:,1)=S1*((1D0/DET_MAT)*(B1*(A1**2D0+A2**2D0+A3**2D0+&
               A4**2D0)+2D0*A1*(A1*B1+A2*B2+A3*B3+A4*B4))- &
               (DDMAT(I,:)/DET_MAT**2D0)*A1*(A1**2D0+A2**2D0+A3**2D0+A4**2D0))
          !print*,minval(dstokesm(i,:,1)),maxval(dstokesm(i,:,1))
       ENDDO
       ! Derivatives of Stokes Q, U, V with respect to
       ! 7 regular free parameters.
       DO K=2,4
          DO I=1,7
             A1 = ETAI
             B1 = DerETAI(I,:)
             A2 = EXTRA
             B2 = DEXTRA(I,:)
             SELECT CASE (K)
             CASE(2)
                C1=ETAQ
                D1=DerETAQ(I,:)
                C2=ETAV
                D2=DerETAV(I,:)
                C3=RHOU
                D3=DerRHOU(I,:)
                C4=ETAU
                D4=DerETAU(I,:)
                C5=RHOV
                D5=DerRHOV(I,:)
                C6=RHOQ
                D6=DerRHOQ(I,:)
             CASE(3)
                C1=ETAU
                D1=DerETAU(I,:)
                C2=ETAQ
                D2=DerETAQ(I,:)
                C3=RHOV
                D3=DerRHOV(I,:)
                C4=ETAV
                D4=DerETAV(I,:)
                C5=RHOQ
                D5=DerRHOQ(I,:)
                C6=RHOU
                D6=DerRHOU(I,:)
             CASE(4)
                C1=ETAV
                D1=DerETAV(I,:)
                C2=ETAU
                D2=DerETAU(I,:)
                C3=RHOQ
                D3=DerRHOQ(I,:)
                C4=ETAQ
                D4=DerETAQ(I,:)
                C5=RHOU
                D5=DerRHOU(I,:)
                C6=RHOV
                D6=DerRHOV(I,:)
             END SELECT
             PART1 = (1D0/DET_MAT)*(2D0*A1*B1*C1+D1*A1**2D0+&
                  B1*(C2*C3-C4*C5)+A1*(C2*D3+D2*C3-D4*C5-C4*D5)+ &
                  D6*A2+C6*B2)
             PART2 = (DDMAT(I,:)/DET_MAT**2D0)*(A1**2D0*C1+A1* &
                  (C2*C3-C4*C5)+C6*EXTRA)
             DSTOKES_MAG(I,:,K) = S1*(PART2-PART1)
          ENDDO
       ENDDO
       ! Derivatives of I, Q, U, V with respect to S0 and S1
       DSTOKES_MAG(8,:,1) = 1D0
       DSTOKES_MAG(8,:,2:4) = 0D0
       DSTOKES_MAG(9,:,1) = (STOKES_MAG(:,1)-S0)/S1
       DSTOKES_MAG(9,:,2:4) = STOKES_MAG(:,2:4)/S1
       !-----------------------------------------------------------------------------
       ! Now we apply HMI Filters
       !-----------------------------------------------------------------------------
       DO M=1,9
          DO K=1,4
             DO J=1,NBINS
                DSYN_MAG(M,J,K)=SUM(FILTERS(:,J)*DSTOKES_MAG(M,:,K))     
             ENDDO
          ENDDO
       ENDDO

       ! By RCE: The derivatives of the filtered Stokes parameters with respect to S0 and S1 have to include the 
       ! derivative of the hack for the wavelength coverage. We basically added (S0+S1)*C to the filtered Stokes I,
       ! where C is the integral of the filters in the outer range of the wavelength vector. So I have to add C
       ! to the derivative of Stokes I with respect to S0 and S1.
       DO J = 1, NBINS
          DSYN(8,J,1) = DSYN(8,J,1) + INTEG_FILTERS(J)
 	  DSYN(9,J,1) = DSYN(9,J,1) + INTEG_FILTERS(J)
       ENDDO
 
       DSYN(1:9,:,:)=ALPHAM*DSYN_MAG
       DSYN(10,:,:)=SYN_MAG-SCAT
    ENDIF
  END SUBROUTINE SYNTHESIS
  !!
  !! SOUBROUTINE ABSMAT
  !!
  SUBROUTINE ABSMAT(MODEL, DERIVATIVE, ETAI, ETAQ, ETAU, ETAV, RHOQ, RHOU, RHOV, &
       DerETAI, DerETAQ, DerETAU, DerETAV, DerRHOQ, DerRHOU, DerRHOV)
    !
    ! J M Borrero
    ! Jan 7, 2007
    ! HAO-NCAR for HMI-Stanford
    !
    USE CONS_PARAM
    USE LINE_PARAM
    USE INV_PARAM
    USE VOIGT_DATA
    IMPLICIT NONE
    !----------------------------------------------------------
    REAL(DP), INTENT(IN), DIMENSION(10)     :: MODEL
    LOGICAL, INTENT(IN)                     :: DERIVATIVE
    !----------------------------------------------------------
    REAL(DP),    DIMENSION(NUMW)   :: FRECR, FRECP, FRECB
    REAL(DP),    DIMENSION(NUMW)   :: PHIR, PHIP, PHIB, PSIR, PSIP, PSIB
    REAL(DP),    DIMENSION(NUMW)   :: DerPHIR_DerDAM, DerPHIR_DerFRECR
    REAL(DP),    DIMENSION(NUMW)   :: DerPSIR_DerDAM, DerPSIR_DerFRECR
    REAL(DP),    DIMENSION(NUMW)   :: DerPHIB_DerDAM, DerPHIB_DerFRECB
    REAL(DP),    DIMENSION(NUMW)   :: DerPSIB_DerDAM, DerPSIB_DerFRECB
    REAL(DP),    DIMENSION(NUMW)   :: DerPHIP_DerDAM, DerPHIP_DerFRECP
    REAL(DP),    DIMENSION(NUMW)   :: DerPSIP_DerDAM, DerPSIP_DerFRECP
    REAL(DP),    DIMENSION(NUMW)   :: DerFRECR_DerVLOS, DerFRECR_DerDLDOP, DerFRECR_DerB
    REAL(DP),    DIMENSION(NUMW)   :: DerFRECB_DerVLOS, DerFRECB_DerDLDOP, DerFRECB_DerB
    REAL(DP),    DIMENSION(NUMW)   :: DerFRECP_DerVLOS, DerFRECP_DerDLDOP, DerFRECP_DerB
    !
    REAL(DP),    DIMENSION(NUMW)   :: ETAI, ETAQ, ETAU, ETAV, RHOQ, RHOU, RHOV
    REAL(DP),    DIMENSION(7,NUMW) :: DerETAI, DerETAQ, DerETAU, DerETAV
    REAL(DP),    DIMENSION(7,NUMW) :: DerRHOQ, DerRHOU, DerRHOV
    !
    REAL(DP)                     :: VLOS, DLDOP, BFIELD, DAM, ETA0, GAM, PHI
    REAL(DP)                     :: SIN2INC, COS2AZI, SIN2AZI, SININC, COSINC, SINCOSINC
    REAL(DP),    DIMENSION(NUMW) :: ABSOR1, ABSOR3, DISPE1, DISPE3
    INTEGER                      :: I
    !
    FRECR(:)=0D0
    FRECP(:)=0D0
    FRECB(:)=0D0
    PHIR(:)=0D0
    PHIP(:)=0D0
    PHIB(:)=0D0
    PSIR(:)=0D0
    PSIP(:)=0D0
    PSIB(:)=0D0
    ETAI(:)=0D0
    ETAQ(:)=0D0
    ETAU(:)=0D0
    ETAV(:)=0D0
    RHOQ(:)=0D0
    RHOU(:)=0D0
    RHOV(:)=0D0
    DerETAI(:,:)=0D0
    DerETAQ(:,:)=0D0
    DerETAU(:,:)=0D0
    DerETAV(:,:)=0D0
    DerRHOQ(:,:)=0D0
    DerRHOU(:,:)=0D0
    DerRHOV(:,:)=0D0
    DerPHIR_DerDAM(:)=0D0
    DerPHIR_DerFRECR(:)=0D0
    DerPSIR_DerDAM(:)=0D0 
    DerPSIR_DerFRECR(:)=0D0
    DerPHIB_DerDAM(:)=0D0
    DerPHIB_DerFRECB(:)=0D0
    DerPSIB_DerDAM(:)=0D0
    DerPSIB_DerFRECB(:)=0D0
    DerPHIP_DerDAM(:)=0D0
    DerPHIP_DerFRECP(:)=0D0
    DerPSIP_DerDAM(:)=0D0
    DerPSIP_DerFRECP(:)=0D0
    DerFRECR_DerVLOS(:)=0D0
    DerFRECR_DerDLDOP(:)=0D0
    DerFRECR_DerB(:)=0D0
    DerFRECB_DerVLOS(:)=0D0 
    DerFRECB_DerDLDOP(:)=0D0 
    DerFRECB_DerB(:)=0D0
    DerFRECP_DerVLOS(:)=0D0
    DerFRECP_DerDLDOP(:)=0D0
    DerFRECP_DerB(:)=0D0
    VLOS=0D0
    DLDOP=0D0
    BFIELD=0D0
    DAM=0D0
    ETA0=0D0
    GAM=0D0
    PHI=0D0
    ABSOR1(:)=0D0
    ABSOR3(:)=0D0
    DISPE1(:)=0D0
    DISPE3(:)=0D0
    
    ! Model parameters: magnetic component
    ETA0=MODEL(1)
    GAM=MODEL(2)*D2R
    PHI=MODEL(3)*D2R
    DAM=MODEL(4)
    DLDOP=MODEL(5)
    BFIELD=MODEL(6)
    VLOS=MODEL(7)
    ! Frecuency arrays
    DO I=1,NUMW
       FRECR(I)=(WAVE(I)-1.D3*VLOS*LANDA0/LIGHT+BFIELD*SHIFT)/DLDOP
       FRECB(I)=(WAVE(I)-1.D3*VLOS*LANDA0/LIGHT-BFIELD*SHIFT)/DLDOP
       FRECP(I)=(WAVE(I)-1.D3*VLOS*LANDA0/LIGHT)/DLDOP
    ENDDO
    ! Absortion-dispersion profiles: slow calculation
    IF (FREE(4).EQ..TRUE.) THEN
       CALL VOIGT(NUMW,DAM,FRECR,PHIR,PSIR)
       CALL VOIGT(NUMW,DAM,FRECB,PHIB,PSIB)
       CALL VOIGT(NUMW,DAM,FRECP,PHIP,PSIP)
    ENDIF
    IF (FREE(4).EQ..FALSE.) THEN
       CALL VOIGT_TAYLOR(DAM,FRECR,PHIR,PSIR)
       CALL VOIGT_TAYLOR(DAM,FRECB,PHIB,PSIB)
       CALL VOIGT_TAYLOR(DAM,FRECP,PHIP,PSIP)
    ENDIF
    ! Common parts
    SIN2INC=DSIN(GAM)**2D0
    COS2AZI=DCOS(2D0*PHI)
    SIN2AZI=DSIN(2D0*PHI)
    SINCOSINC=DSIN(GAM)*DCOS(GAM)
    SININC=DSIN(GAM)
    COSINC=DCOS(GAM)
    ABSOR1=PHIP-0.5D0*(PHIB+PHIR)
    DISPE1=PSIP-0.5D0*(PSIB+PSIR)
    ABSOR3=PHIR-PHIB
    DISPE3=PSIR-PSIB
    !
    ! ETAI:
    !
    ETAI = 1D0+(ETA0/2D0)*(PHIP*SIN2INC+0.5D0*(PHIR+PHIB)*(2D0-SIN2INC))
    !
    ! ETAQ:
    !
    ETAQ = 0.5D0*ETA0*ABSOR1*SIN2INC*COS2AZI
    !
    ! ETAU:
    !
    ETAU = 0.5D0*ETA0*ABSOR1*SIN2INC*SIN2AZI
    !
    ! ETAV:
    !
    COSINC=DCOS(GAM)
    ETAV = -0.5D0*ETA0*ABSOR3*COSINC
    !
    ! RHOQ:
    !
    RHOQ = 0.5D0*ETA0*DISPE1*SIN2INC*COS2AZI
    !
    ! RHOU:
    !
    RHOU = 0.5D0*ETA0*DISPE1*SIN2INC*SIN2AZI
    !
    ! RHOV:
    !
    RHOV = -0.5D0*ETA0*DISPE3*COSINC
    
    IF (DERIVATIVE.EQ..TRUE.) THEN
       
       ! Derivatives of absortion-dispersion profiles
       ! Sigma-Red component
       DerPHIR_DerDAM = -2D0/DSQRT(DPI)+2D0*(DAM*PHIR+FRECR*PSIR)
       DerPHIR_DerFRECR = 2D0*DAM*PSIR-2D0*FRECR*PHIR
       DerPSIR_DerDAM = DerPHIR_DerFRECR
       DerPSIR_DerFRECR = -DerPHIR_DerDAM
       ! Sigma-Blue component
       DerPHIB_DerDAM = -2D0/DSQRT(DPI)+2D0*(DAM*PHIB+FRECB*PSIB)
       DerPHIB_DerFRECB = 2D0*DAM*PSIB-2D0*FRECB*PHIB
       DerPSIB_DerDAM = DerPHIB_DerFRECB
       DerPSIB_DerFRECB = -DerPHIB_DerDAM
       ! Simga-Pi component
       DerPHIP_DerDAM = -2D0/DSQRT(DPI)+2D0*(DAM*PHIP+FRECP*PSIP)
       DerPHIP_DerFRECP = 2D0*DAM*PSIP-2D0*FRECP*PHIP
       DerPSIP_DerDAM = DerPHIP_DerFRECP
       DerPSIP_DerFRECP = -DerPHIP_DerDAM 


       ! Derivatives of the frecuency with respect to
       ! the field strength, LOS velocity, Doppler width.
       DerFRECR_DerVLOS = -1D3*LANDA0/(LIGHT*DLDOP)
       DerFRECB_DerVLOS = DerFRECR_DerVLOS
       DerFRECP_DerVLOS = DerFRECR_DerVLOS
       !
       DerFRECR_DerDLDOP = - FRECR/DLDOP
       DerFRECB_DerDLDOP = - FRECB/DLDOP
       DerFRECP_DerDLDOP = - FRECP/DLDOP
       !
       DerFRECR_DerB = SHIFT/DLDOP
       DerFRECB_DerB = - DerFRECR_DerB
       DerFRECP_DerB(:) = 0D0 
       ! DerFRECP_DerB won't be used to speed up and reduce noise.
       !
       ! ETAI derivatives:
       !
       DerETAI(1,:) = (ETAI(:)-1D0)/ETA0
       DerETAI(2,:) = ETA0*SINCOSINC*ABSOR1*D2R
       DerETAI(3,:) = 0D0
       DerETAI(4,:) = 0.5D0*ETA0*(DerPHIP_DerDAM*SIN2INC+0.5D0*(DerPHIB_DerDAM+ &
            DerPHIR_DerDAM)*(2D0-SIN2INC))
       DerETAI(5,:) = 0.5D0*ETA0*(DerPHIP_DerFRECP*DerFRECP_DerDLDOP*SIN2INC+ &
            0.5D0*(DerPHIB_DerFRECB*DerFRECB_DerDLDOP+DerPHIR_DerFRECR* &
            DerFRECR_DerDLDOP)*(2D0-SIN2INC))
       DerETAI(6,:) = 0.25D0*ETA0*(DerPHIB_DerFRECB*DerFRECB_DerB+DerPHIR_DerFRECR* &
            DerFRECR_DerB)*(2D0-SIN2INC)
       DerETAI(7,:) = 0.5D0*ETA0*(DerPHIP_DerFRECP*DerFRECP_DerVLOS*SIN2INC+ &
            0.5D0*(DerPHIB_DerFRECB*DerFRECB_DerVLOS+DerPHIR_DerFRECR* &
            DerFRECR_DerVLOS)*(2D0-SIN2INC))
       !
       ! ETA Q derivatives
       !
       DerETAQ(1,:) = ETAQ/ETA0
       DerETAQ(2,:) = SINCOSINC*COS2AZI*ETA0*ABSOR1*D2R
       DerETAQ(3,:) = -ETA0*ABSOR1*SIN2INC*SIN2AZI*D2R
       DerETAQ(4,:) = 0.5D0*ETA0*(DerPHIP_DerDAM-0.5D0*(DerPHIB_DerDAM+ &
            DerPHIR_DerDAM))*SIN2INC*COS2AZI
       DerETAQ(5,:) = 0.5D0*ETA0*(DerPHIP_DerFRECP*DerFRECP_DerDLDOP- &
            0.5D0*(DerPHIB_DerFRECB*DerFRECB_DerDLDOP+DerPHIR_DerFRECR* &
            DerFRECR_DerDLDOP))*SIN2INC*COS2AZI
       DerETAQ(6,:) = -0.25D0*ETA0*(DerPHIB_DerFRECB*DerFRECB_DerB+ &
            DerPHIR_DerFRECR*DerFRECR_DerB)*SIN2INC*COS2AZI
       DerETAQ(7,:) = 0.5D0*ETA0*(DerPHIP_DerFRECP*DerFRECP_DerVLOS- &
            0.5D0*(DerPHIB_DerFRECB*DerFRECB_DerVLOS+DerPHIR_DerFRECR* &
            DerFRECR_DerVLOS))*SIN2INC*COS2AZI
       !
       ! ETAU derivatives
       !
       DerETAU(1,:) = ETAU/ETA0
       DerETAU(2,:) = SINCOSINC*SIN2AZI*ETA0*ABSOR1*D2R
       DerETAU(3,:) = ETA0*ABSOR1*SIN2INC*COS2AZI*D2R
       DerETAU(4,:) = 0.5D0*ETA0*(DerPHIP_DerDAM-0.5D0*(DerPHIB_DerDAM+ &
            DerPHIR_DerDAM))*SIN2INC*SIN2AZI
       DerETAU(5,:) = 0.5D0*ETA0*(DerPHIP_DerFRECP*DerFRECP_DerDLDOP- &
            0.5D0*(DerPHIB_DerFRECB*DerFRECB_DerDLDOP+DerPHIR_DerFRECR* &
            DerFRECR_DerDLDOP))*SIN2INC*SIN2AZI
       DerETAU(6,:) = -0.25D0*ETA0*(DerPHIB_DerFRECB*DerFRECB_DerB+ &
            DerPHIR_DerFRECR*DerFRECR_DerB)*SIN2INC*SIN2AZI
       DerETAU(7,:) = 0.5D0*ETA0*(DerPHIP_DerFRECP*DerFRECP_DerVLOS- &
            0.5D0*(DerPHIB_DerFRECB*DerFRECB_DerVLOS+DerPHIR_DerFRECR* &
            DerFRECR_DerVLOS))*SIN2INC*SIN2AZI
       !
       ! ETAV derivatives
       !
       DerETAV(1,:) = ETAV/ETA0
       DerETAV(2,:) = 0.5D0*ETA0*ABSOR3*SININC*D2R
       DerETAV(3,:) = 0D0
       DerETAV(4,:) = -0.5D0*ETA0*(DerPHIR_DerDAM-DerPHIB_DerDAM)*COSINC
       DerETAV(5,:) = -0.5D0*ETA0*(DerPHIR_DerFRECR*DerFRECR_DerDLDOP- &
            DerPHIB_DerFRECB*DerFRECB_DerDLDOP)*COSINC
       DerETAV(6,:) = -0.5D0*ETA0*(DerPHIR_DerFRECR*DerFRECR_DerB- &
            DerPHIB_DerFRECB*DerFRECB_DerB)*COSINC
       DerETAV(7,:) = -0.5D0*ETA0*(DerPHIR_DerFRECR*DerFRECR_DerVLOS- &
            DerPHIB_DerFRECB*DerFRECB_DerVLOS)*COSINC
       
       !
       ! RHOQ derivatives:
       !
       DerRHOQ(1,:) = RHOQ/ETA0
       DerRHOQ(2,:) = SINCOSINC*COS2AZI*ETA0*DISPE1*D2R
       DerRHOQ(3,:) = -ETA0*DISPE1*SIN2INC*SIN2AZI*D2R
       DerRHOQ(4,:) = 0.5D0*ETA0*(DerPSIP_DerDAM-0.5D0*(DerPSIB_DerDAM+ &
            DerPHIR_DerDAM))* SIN2INC*COS2AZI
       DerRHOQ(5,:) = 0.5D0*ETA0*(DerPSIP_DerFRECP*DerFRECP_DerDLDOP- &
            0.5D0*(DerPSIB_DerFRECB*DerFRECB_DerDLDOP+DerPSIR_DerFRECR* &
            DerFRECR_DerDLDOP))*SIN2INC*COS2AZI
       DerRHOQ(6,:) = -0.25D0*ETA0*(DerPSIB_DerFRECB*DerFRECB_DerB+ &
            DerPSIR_DerFRECR*DerFRECR_DerB)*SIN2INC*COS2AZI
       DerRHOQ(7,:) = 0.5D0*ETA0*(DerPSIP_DerFRECP*DerFRECP_DerVLOS- &
            0.5D0*(DerPSIB_DerFRECB*DerFRECB_DerVLOS+DerPSIR_DerFRECR* &
            DerFRECR_DerVLOS))*SIN2INC*COS2AZI
       !
       ! RHOU derivatives
       !
       
       DerRHOU(1,:) = RHOU/ETA0
       DerRHOU(2,:) = SINCOSINC*SIN2AZI*ETA0*DISPE1*D2R
       DerRHOU(3,:) = ETA0*DISPE1*SIN2INC*COS2AZI*D2R
       DerRHOU(4,:) = 0.5D0*ETA0*(DerPSIP_DerDAM-0.5D0*(DerPSIB_DerDAM+ &
            DerPSIR_DerDAM))*SIN2INC*SIN2AZI
       DerRHOU(5,:) = 0.5D0*ETA0*(DerPSIP_DerFRECP*DerFRECP_DerDLDOP- &
            0.5D0*(DerPSIB_DerFRECB*DerFRECB_DerDLDOP+DerPSIR_DerFRECR* &
            DerFRECR_DerDLDOP))*SIN2INC*SIN2AZI
       DerRHOU(6,:) = -0.25D0*ETA0*(DerPSIB_DerFRECB*DerFRECB_DerB+ &
            DerPSIR_DerFRECR*DerFRECR_DerB)*SIN2INC*SIN2AZI
       DerRHOU(7,:) = 0.5D0*ETA0*(DerPSIP_DerFRECP*DerFRECP_DerVLOS- &
            0.5D0*(DerPSIB_DerFRECB*DerFRECB_DerVLOS+DerPSIR_DerFRECR* &
            DerFRECR_DerVLOS))*SIN2INC*SIN2AZI
       !
       ! RHOV derivatives
       !
       RHOV = -0.5D0*ETA0*DISPE3*COSINC
       DerRHOV(1,:) = RHOV/ETA0
       DerRHOV(2,:) = 0.5D0*ETA0*DISPE3*SININC*D2R
       DerRHOV(3,:) = 0D0
       DerRHOV(4,:) = -0.5D0*ETA0*(DerPSIR_DerDAM-DerPSIB_DerDAM)*COSINC
       DerRHOV(5,:) = -0.5D0*ETA0*(DerPSIR_DerFRECR*DerFRECR_DerDLDOP- &
            DerPSIB_DerFRECB*DerFRECB_DerDLDOP)*COSINC
       DerRHOV(6,:) = -0.5D0*ETA0*(DerPSIR_DerFRECR*DerFRECR_DerB- &
            DerPSIB_DerFRECB*DerFRECB_DerB)*COSINC
       DerRHOV(7,:) = -0.5D0*ETA0*(DerPSIR_DerFRECR*DerFRECR_DerVLOS- &
            DerPSIB_DerFRECB*DerFRECB_DerVLOS)*COSINC
    ENDIF
    
  ENDSUBROUTINE ABSMAT
END MODULE FORWARD

