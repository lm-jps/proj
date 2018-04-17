subroutine doppcal_estimate(naxis1, naxis2, naxis3i, naxis3o, & ! N_x, N_y; naxis3 = 7 unneeded?
     data_in, data_out, & ! input and output arrays
     crpix, obs_v, rsun_obs, rsun_ref, crlt_obs, crota2, & ! input HMI keywords
     doppcal_bias, & ! output estimate of bias velocity
     max_ang_in, thresh_blos_in, pix_in, &  ! optional input parameters
     thresh_bmag_in, rad_los_sep_in, diff_a0_in, diff_a2_in, diff_a4_in)  ! more optional input parameters

  
  !===================================================================
  ! This subroutine uses the median of Doppler shifts on LoS PILs
  ! to estimate zero-point offset in Doppler shift.
  !===================================================================
  !
  ! "naxis1","naxis2","naxis3" should be a integer(8) array dimensions 
  !
  ! "data_in" should be a real(8) (NAXIS1, NAXIS2, NAXIS3) array of 
  !     naxis3 = 7 maps, each a (4096x4096) array, corresponding to: 
  !       1. azimuth
  !       2. conf_disambig
  !       3. confid
  !       4. disambig
  !       5. field
  !       6. inclination
  !       7. vlos [uncorrected for S/C, diff. rot., or convective blueshift]
  !
  ! "crpix" is a real*8, two-element array of (x,y) pix. address 
  !     of disk center 
  !
  ! "obs_v" is a real*8, three-element array of observer's (HMI's)
  !     velocity (in m/s) in radial, west, and north directions
  !
  ! "rsun_obs" is a real*8, equal to sun's radius in arcsec
  !     default pix. size = HMI's (0.504276''), others poss. (see below)
  !
  ! "rsun_ref" is a real*8, equal to sun's radius in meters,
  !     used for removing rotation
  !
  ! "crlt_obs" is a real*8, equal to Heliogrph Latitude of 
  !     observer's disk center, to remove differential rotation.
  !
  !  Uses Snodgrass '82/'84 coeffs from Rick Bogart's mtrack.c -- quote:
  !      a0      float   -0.02893	Coefficients in sin^2 (latitude)
  !      a2      float   -0.3441	expansion of rotation rate minus
  !      a4      float   -0.5037	Carrington rotation (urad/sec)
  !
  ! "doppcal_bias" is OUTPUT, real*8; is D.C. bias of Doppler vel,
  !    subtract for true
  !
  ! "data_out" is OUTPUT, real*4 (single precision, to keep file size
  !    small), (naxis1,naxis2,3) array of spacecraft and
  !    solar-differential-rotation-corrected Dopplergram; map of LOS
  !    PILs; map of radial-field PILs
  !
  ! The following are optional, real*8 input parameters ("keywords"): 
  !------------------------------------------------------------------
  !   "max_ang_in" = greatest angular distance (in degrees) from disk 
  !        center for which PILs will be identified; default = 60.0
  !
  !   "pix_in" = user-set pixel size (arcsec), default = HMI = 0.504276
  !
  !   "thresh_blos_in" = threshold LOS magnetic flux density above which  
  !        fields are treated as significantly (statistically) differ-
  !        end from zero (i.e., truly + or -); default: 25 Mx cm^2
  !
  !   "thresh_bmag_in" = threshold LOS magnetic flux density above which  
  !        fields are treated as significantly (statistically) differ-
  !        end from zero (i.e., truly + or -); default: 25 Mx cm^2
  !
  !   "rad_los_sep_in" = maximum separation allowed between LoS PIL
  !        pixels and nearest radial-field PIL
  !
  !   "diff_a0, diff_a2, diff_a4" are a0/a2/a4 above. added X Sun Apr 16 2018
  !
  !==========================================================
  !
  ! METHOD:
  !   (1) Compute pixels' heliographic latitude & longitudes
  !   (2) W/observer's lat/lon, remove spacecraft vel. from V_los 
  !   (3) W/Sun's lat/lon, remove differ'l rotation from V_los
  !   (4) Get pix. on PILs of B_LOS field within "max_ang" of disk
  !       center. To do so:
  !        - Create binary 2D masks of +/- LoS pixels above "blos_thr"
  !        - Use dilate.f90 to shift (+) mask in 4 directions, same w/(-), 
  !          the where product of masks is nonzero = LoS PIL  
  !   (5) Get pix. on PILs of B_radial within "max_ang" of disk
  !       center. To do so:
  !        - Create binary 2D masks of +/- pixels w/|B| above "bmag_thr"
  !          and |B_r| above "blos_thr"
  !        - Use dilate.f90 to shift (+) mask in 4 directions, same w/(-), 
  !          the where product of masks is nonzero = radial-field PIL  
  !   (6) Find median of vlos in all pixels on LoS PILs within
  !       "rad_los_sep" pixels on radial-field PILs  
  !   (7) Find median of vlos in all pixels on LoS+radial-field PILs  
  !
  ! HISTORY:
  !  2016/06/29, BT Welsch: Started adaption from various IDL codes.
  !  2017/07/31, BTW: Started wholesale import from doppcal_test.f90
  !  2017/08/08, BTW: Added array sizes in calling arguments
  !  2017/09/19, BTW: Finishing up, I hope!
  !   
  !=====================================================================
  
  implicit none
  ! Input
  integer*4, intent(in) :: naxis1, naxis2, naxis3i, naxis3o !, naxis3 (unneeded?)
  real*8, intent(in) :: data_in(naxis1,naxis2,naxis3i)
  real*4, intent(inout) :: data_out(naxis1,naxis2,naxis3o)
  real*4, intent(inout) :: doppcal_bias ! bias of Doppler vel, subtract for true

  ! Define optional inputs
  real*8, intent(in), optional :: max_ang_in, thresh_blos_in, thresh_bmag_in
  real*8, intent(in), optional :: diff_a0_in, diff_a2_in, diff_a4_in
  real*8, intent(in), optional :: pix_in
  integer, intent(in), optional :: rad_los_sep_in 
  !
  real*8, intent(in) :: crpix(2) ! 2-element pix. address of disk center
  real*8, intent(in) :: obs_v(3) ! 3-element spacecraft vel. (R,N,W)
  real*8, intent(in) :: rsun_obs ! radius of the sun in arcsec, to HMI
  real*8, intent(in) :: rsun_ref ! radius of the sun in m, 696,000,000
  real*8, intent(in) :: crlt_obs ! Heliogrph Latit. of observer's disk ctr
  real*8, intent(in) :: crota2   ! HMI uses (-) P-angle, so...
  !                             ...(+)=solar N west of observer N
    
  real*8, allocatable :: blos(:,:), azim(:,:), brad(:,:)  
  real*8, allocatable :: vlos(:,:)
  real*8, allocatable :: lonlat(:,:,:), lonlat_noBP(:,:,:) ! 
  integer, allocatable :: pilmap_los(:,:), pilmap_rad(:,:) !
  integer, allocatable :: pilmap_rad_dil(:,:), pilmap_both(:,:)
  real*8, allocatable :: vlos_pils(:)
    
  real*8 pix ! pixel size in arc sec
  real*8 rsun_pix ! radius of the sun in pixels
  real*8 max_ang ! only use data < max_ang degrees from disk ctr; default=60
  real*8 dtor ! degrees-to-radians conversion
  real*8 obsv_x, obsv_y ! projections of SDO velocity onto CCD x & y axes

  integer dilation_param ! max # pix between +/- still called a PIL (use 1)
  integer*8 pil_count
  real*8 thresh_blos ! use flux densities above this (use |B_los|.ge.60)
  real*8 thresh_bmag ! use pix w/avg |B| above this (use |B|.ge.250)
  integer rad_los_sep
  
  real*8 vlos_pils_median ! median of vlos on PILs

  real*8 diff_a0, diff_a2, diff_a4

  real*8, parameter :: pi = 3.1415927 ! 3.1415926535897932

  integer alloc_count, dealloc_count

  dtor = pi/180

  alloc_count = 0 
  dealloc_count = 0 

  if (present(max_ang_in)) then ! max ang dist [deg.] from disk ctr
     max_ang = max_ang_in 
  else 
     max_ang = 60.0
  end if

  if (present(thresh_blos_in)) then ! B_LOS above this assumed not noise
     thresh_blos = thresh_blos_in 
  else 
     thresh_blos = 60.0
  end if

  if (present(thresh_bmag_in)) then ! |B| above this assumed not noise
     thresh_bmag = thresh_bmag_in 
  else 
     thresh_bmag = 250.0
  end if

  if (present(rad_los_sep_in)) then ! max sep [pix] btwn B_r & B_L PILs
     rad_los_sep = rad_los_sep_in
  else 
     rad_los_sep = 2
  end if

  if (present(pix_in)) then ! pixel size in arcsec
     pix = pix_in
  else 
     pix = 0.504276     
  end if

  if (present(diff_a0_in)) then ! pixel size in arcsec
     diff_a0 = diff_a0_in
  else
     diff_a0 = 2.71390
  end if

  if (present(diff_a2_in)) then ! pixel size in arcsec
     diff_a2 = diff_a2_in
  else
     diff_a2 = - 0.40500
  end if

  if (present(diff_a4_in)) then ! pixel size in arcsec
     diff_a4 = diff_a4_in
  else
     diff_a4 = - 0.42200
  end if

!  print *, diff_a0, diff_a2, diff_a4

  ! Alloc. intermediate variables 
  allocate( lonlat( naxis1, naxis2, 2)) 
  allocate( lonlat_noBP( naxis1, naxis2, 2)) 

  allocate( pilmap_los( naxis1, naxis2)) 
  allocate( pilmap_rad( naxis1, naxis2)) 
  allocate( pilmap_rad_dil( naxis1, naxis2))
  allocate( pilmap_both( naxis1, naxis2))

  allocate( blos( naxis1, naxis2)) 
  allocate( azim( naxis1, naxis2)) 
  allocate( brad( naxis1, naxis2)) 

  allocate( vlos( naxis1, naxis2)) 

  alloc_count = alloc_count + 10

  ! Get observer's & Sun's longit & latit
  !========================================
  rsun_pix = rsun_obs/pix ! default HMI pix size (arcsec)

  call pix2helio(naxis1, naxis2, lonlat_noBP, crpix, rsun_pix, 0.d0, 0.d0, max_ang)

  call pix2helio(naxis1, naxis2, lonlat, crpix, rsun_pix, crlt_obs, crota2, max_ang)

  ! Get PILs
  !========================================
  azim = data_in(:,:,1)

  ! need conf_disambig of 60 or 90, and within max_ang of disk center
  where ((data_in(:,:,2).ge.60).and.(lonlat_noBP(:,:,1).ne.0.)) 

     blos = data_in(:,:,5)*cos(data_in(:,:,6)*dtor)

     where (data_in(:,:,4).eq.7) 
        
        azim = azim - 180.  ! disambig=7 ==> add 180
        
     end where
     
     ! This assumes azim=0 is straight down, increasing CCW
     ! note (-) for sin longitude in "upside-down" noBP frame (line 3)
     ! note (+) for sin latitude in "upside-down" noBP frame (line 4)
     brad = blos*cos(lonlat_noBP(:,:,2))*cos(lonlat_noBP(:,:,1)) &
      + data_in(:,:,5)*sin(data_in(:,:,6)*dtor)*( &
      - sin(azim*dtor)*cos(lonlat_noBP(:,:,2))*sin(lonlat_noBP(:,:,1)) &
      + cos(azim*dtor)*sin(lonlat_noBP(:,:,2)) )

  elsewhere

     blos = 0.
     brad = 0.

  end where

  dilation_param = 1

  call get_pils_los(naxis1, naxis2, blos, pilmap_los, thresh_blos, dilation_param)

  call get_pils_rad(naxis1, naxis2, brad, data_in(:,:,5), pilmap_rad, &
       thresh_blos, thresh_bmag, dilation_param)

  ! Remove spacecraft vel., using lonlat_noBP
  !===========================================

  ! To account for HMI being *almost* upside-down, apply rotation
  ! matrix (thru P-angle) to obs_v's west (W=right) & north (N) 
  ! components, which are heliocentric.
  !
  ! NOTE HMI P SIGN: (+)ive when camera is rotated CCW w/r.t. fixed Sun, 
  !                  or Sun's pole is rotated clockwise w/r.t. the frame
  !
  ! Result is v_x & v_y in HMI's CCD coords (x=right, y=up)

  ! obs_v(3) is (R,N,W)
  obsv_x =  obs_v(3)*cos(crota2*dtor) + obs_v(2)*sin(crota2*dtor)
  obsv_y = -obs_v(3)*sin(crota2*dtor) + obs_v(2)*cos(crota2*dtor)

  vlos = data_in(:,:,7)
 
  where (lonlat_noBP(:,:,1).ne.0.) ! only inside "max_ang" from disk ctr

     ! remove projection of spacecraft's 3-component obs_v onto LoS
     ! r-comp (+ = away from Sun), drop 2d-order dependence on (x,y) 
     ! y-comp (note 1 AU = 215 R_sun)
     ! x-comp
     vlos = vlos & ! obs_v in m/s, so 100 to convert from m/s --> cm/s 
     - 100.*obs_v(1) & ! 
     + 100.*obsv_y*sin(lonlat_noBP(:,:,2))/215. &
     + 100.*obsv_x*sin(lonlat_noBP(:,:,1))*cos(lonlat_noBP(:,:,2))/215.

     ! Following mtrack.c, remove diff. rot'n (use lonlat with B & P)
     !=================================================================
     !  a0   float   -0.02893	Coefficients in sin^2 (latitude)
     !  a2   float   -0.3441	expansion of rotation rate minus
     !  a4   float   -0.5037	Carrington rotation (urad/sec)
     !
     !  define CARR_RATE       (2.86532908457)
     !  vrot = (a0 + CARR_RATE) *  RSunMm * coslat * uonlos;
     !  vrot += a2 *  RSunMm * coslat * uonlos * sin2lat;
     !  vrot += a4 *  RSunMm * coslat * uonlos * sin2lat * sin2lat;

     ! 1e+2 for rsun_ref in m --> cm & 1e-6 for microradians 
     ! UNCOMMENT THE FOLLOWING FOUR LINES TO USE MTRACK APPROACH
     !     vlos = vlos -rsun_ref *1e-4 *sin(lonlat(:,:,1))*cos(crlt_obs*dtor) &
     !     *( (-0.02893 + 2.86532908457) &    ! a0 term
     !     - 0.3441 * sin(lonlat(:,:,2))**2 & ! a2 term
     !     - 0.5037 * sin(lonlat(:,:,2))**4 ) ! a4 term

  
     ! From /cvs/JSOC/proj/lev1.5_hmi/libs/lev15/rotcoef_file.txt 
     !=================================================================
     !  a0   float    2.71390	Coefficients in sin^2 (latitude)
     !  a2   float   -0.405000	expansion of rotation (urad/sec) 
     !  a4   float   -0.422000	
     !
     ! 1e+2 for rsun_ref in m --> to cm; and 1e-6 for microradians 
     vlos = vlos - rsun_ref *1e-4 *sin(lonlat(:,:,1)) *cos(crlt_obs*dtor) &
          *( diff_a0 &    ! a0 term
          + diff_a2 * sin(lonlat(:,:,2))**2 & ! a2 term
          + diff_a4 * sin(lonlat(:,:,2))**4 ) ! a4 term

  elsewhere

     vlos = 0.

  end where

  ! Dilate radial-field PILs
  call dilate(naxis1, naxis2, pilmap_rad, pilmap_rad_dil, rad_los_sep)

  pilmap_both = pilmap_los*pilmap_rad_dil ! mask of B_LOS PILs near B_r PILs

  ! Remove pixels where vlos has been masked to zero
  where (vlos .eq. 0)
    pilmap_both = 0
  end where

  pil_count = sum(pilmap_both) ! # of identified PILs

  allocate( vlos_pils(pil_count))

  vlos_pils = pack( vlos*pilmap_both, vlos*pilmap_both /= 0)

  call median_sub(vlos_pils, pil_count, vlos_pils_median)

  deallocate(vlos_pils)
  
  doppcal_bias = real(vlos_pils_median, 4) ! in cm/s

!  print *,'Estimated Doppler bias, from PIL median (cm/s):',doppcal_bias

  data_out(:,:,1) = real(vlos(:,:), 4)
  data_out(:,:,2) = real(pilmap_los(:,:), 4)
  data_out(:,:,3) = real(pilmap_rad(:,:), 4)

  deallocate( lonlat) 
  deallocate( lonlat_noBP) 

  deallocate( pilmap_los)
  deallocate( pilmap_rad)
  deallocate( pilmap_rad_dil)
  deallocate( pilmap_both)

  deallocate( blos)
  deallocate( azim)
  deallocate( brad)

  deallocate( vlos)

  dealloc_count = dealloc_count + 10
  
  ! Check for matched allocations and deallocations.
  check_alloc: if (alloc_count /= dealloc_count) then
     print *,'Found allocate/ deallocate mismatch.'
     print *,'  alloc_count:',alloc_count
     print *,'dealloc_count:',dealloc_count
  endif check_alloc

end subroutine doppcal_estimate
