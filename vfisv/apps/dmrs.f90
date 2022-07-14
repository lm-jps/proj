PROGRAM DMRS
  !
  ! J M Borrero
  ! Dec 14, 2009
  ! HAO-NCAR for HMI-Stanford
  !
  IMPLICIT NONE
  INTEGER, PARAMETER             :: NPROF = 10000  ! NPROF is the variable name,interger is the type specifier and here NPROF is the name for integer value 10000 
  REAL(8), DIMENSION(24,NPROF)   :: OBS  !OBS is a real variable of dimension (24,10000) not yet assigned any values
  REAL(8), DIMENSION(24,NPROF)   :: SCAT 
  REAL(8), DIMENSION(10,NPROF)   :: GUESS, RES 
  REAL(8), DIMENSION(11,NPROF)    :: ERR 
  INTEGER                        :: I, M  ! I,M are integers
  INTEGER, DIMENSION(10)         :: LIST_FREE_PARAMS
  !--------------------------------------------------------------------------
  REAL(8)                        :: LAMBDA_0, LAMBDA_B, NOISE_LEVEL
  INTEGER                        :: NUM_LAMBDA 
  REAL(8)                        :: LAMBDA_MIN, DELTA_LAMBDA
  INTEGER                        :: NUM_LAMBDA_FILTER, NUM_TUNNING, CONTINUUM
  REAL(8)                        :: LYOTFWHM, WNARROW, WSPACING
  INTEGER                        :: NUM_ITERATIONS
  REAL(8)                        :: SVD_TOLERANCE, CHI2_STOP, POLARIZATION_THRESHOLD
  REAL(8)                        :: INTENSITY_THRESHOLD, PERCENTAGE_JUMP
  !--------------------------------------------------------------------------
  ! Some defaults
  !--------------------------------------------------------------------------
  NUM_ITERATIONS = 30
  SVD_TOLERANCE = 1D-32  ! Exponentials declared with D rather than E for the exponent are immediately of kind double-precision real by default
  CHI2_STOP = 1D-6
  POLARIZATION_THRESHOLD = 1E-2
  INTENSITY_THRESHOLD = 0.8D0
  PERCENTAGE_JUMP = 10D0
  LAMBDA_MIN = -432D0
  !---------------------------------------------------------------------------
  ! Defaults for Fe I 6302:
  !---------------------------------------------------------------------------
  NUM_LAMBDA = 33
  LAMBDA_0 = 6302.4936D0
  LAMBDA_B = 0.046243577D0
  DELTA_LAMBDA = 27.0D0
  NOISE_LEVEL = 3D-3
  NUM_LAMBDA_FILTER = 6
  NUM_TUNNING = 6
  CONTINUUM = 0
  LYOTFWHM = 433.0D0
  WNARROW = 176.0D0
  WSPACING = 70.4D0
  !----------------------------------------------------------------------------
  ! Defaults for Fe I 6173:
  !----------------------------------------------------------------------------
  !NUM_LAMBDA = 33
  !LAMBDA_0 = 6173.334D0
  !LAMBDA_B = 0.046243577D0
  !DELTA_LAMBDA = 27.0D0
  !NOISE_LEVEL = 3E-3
  !NUM_LAMBDA_FILTER = 6
  !NUM_TUNNINNG = 6
  !CONTINUUM = 0
  !LYOTFWHM = 424.0D0
  !WNARROW = 172.0D0
  !WSPACING = 69.0D0
  !----------------------------------------------------------------------------
  ! Free Parameters. We set defaults here,but they are also allowed to be 
  ! changed from the command line. All it is needed is to pass LIST_FREE_PARAMS = [0,0,1,1,0, etc]
  ! to free init.f90
  !----------------------------------------------------------------------------
  LIST_FREE_PARAMS(:)=(/1,1,1,0,1,1,1,1,1,1/)
  ! eta0, gamma, phi, damping, dopplerw, bfield, vlos, s0, s1, alpha_mag
  !----------------------------------------------------------------------------
  ! This is the initialization part:
  ! When writting dmrs.f90 into dmrs.c calls should be
  ! subtituted from (e.g.) CALL WAVE_INIT to WAVE_INIT_
  !----------------------------------------------------------------------------
  CALL LINE_INIT (LAMBDA_0, LAMBDA_B, NOISE_LEVEL)
  CALL WAVE_INIT (LAMBDA_MIN, DELTA_LAMBDA, NUM_LAMBDA)
  CALL FILT_INIT (NUM_LAMBDA_FILTER, NUM_TUNNING, CONTINUUM, LYOTFWHM, WNARROW, WSPACING)
  CALL INV_INIT(NUM_ITERATIONS, SVD_TOLERANCE, CHI2_STOP, POLARIZATION_THRESHOLD &
       , INTENSITY_THRESHOLD, PERCENTAGE_JUMP)
  CALL FREE_INIT(LIST_FREE_PARAMS)
  CALL SVD_INIT
  IF (LIST_FREE_PARAMS(4).EQ.0) CALL VOIGT_INIT
  !------------------------------------------------------------------
  ! Somewhere here the DMRS should read:
  ! OBS: observations for each pixel
  ! SCAT: scattered light profile for each pixel
  ! GUESS: initial guess atmosphere (read from previous inversion?)
  !------------------------------------------------------------------
  ! Example of observations. Real ones have to be given by the DMRS
  OBS(:,:)=0D0
  OBS(1:6,:)=1D0
  SCAT(:,:)=0D0
  GUESS(:,:)=0D0
  !------------------------------------------------------------------
  ! This is the Loop that needs to be paralellized with MPICH-2
  !------------------------------------------------------------------
  DO I=1,NPROF
     ! Example of an intial guess.
     ! This is fact has to be provided by the DMRS
     GUESS(1,I)=15D0
     GUESS(2,I)=90D0
     GUESS(3,I)=45D0
     GUESS(4,I)=0.5D0
     GUESS(5,I)=40D0
     GUESS(6,I)=500D0
     GUESS(7,I)=0D0
     GUESS(8,I)=0.4D0
     GUESS(9,I)=0.6D0
     GUESS(10,I)=1D0
     ! When written in C it should read VFISV_
     CALL VFISV(OBS(:,I),SCAT(:,I),GUESS(:,I),RES(:,I),ERR(:,I))
  ENDDO
  ! RES are in the same order as LIST_FREE_PARAMS
  ! ERR are: sigma_B^2, sigma_gamma^2, sigma_phi^2, sigma_vlos^2, sigma_alpham^2, sigma_(b,gamma), sigma)_(b,phi), sigma_(gamma,phi)
  PRINT*,'Finished'
  CALL FREE_MEMORY

END PROGRAM DMRS
