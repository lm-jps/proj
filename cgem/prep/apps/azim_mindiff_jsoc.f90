subroutine azim_mindiff_jsoc(azim_in, azim_out, nx, ny, nstep, flip_thr)

  !=================================================================
  ! This subroutine looks for back-and-forth flips in azimuth in the 
  ! "central" image that are large, perhaps due to ambiguity resolution
  ! error, and reverses identified flips *IF* doing so minimizes 
  ! the mean of unsigned changes in azimuth over (2*nstep_in+1) steps. 
  !=================================================================
  !
  ! "azim_in" should be a real (NAXIS1, NAXIS2, 2*NSTEP+1) array of 
  ! azimuths, corresponding to a set of (2*NSTEP+1) frames of 
  ! (N_x, N_y) images.
  !
  ! "azim_out" is one (N_x, N_y) image, corresponding to time 
  ! NSTEP
  !
  !==========================================================
  !
  ! METHOD:
  !   (1) Compute frame-to-frame difference arrays.  
  !   (2) Identify pix w/aziuth changes \Delta phi that both:
  !        - increase mean, unsigned azimuth change over +/- nstep
  !        - are larger than a threhsold; default, > 120 deg
  !
  ! HISTORY: 2013/07/25, BT Welsch: Hacked from IDL azim_thresh_mindiff.pro
  !          2016/01/14, BT Welsch: Now acts only on 2*NSTEP_IN+1 frames.
  !			 2016/10/11, X Sun: modified input
  ! 
  !=====================================================================

  implicit none
  ! Input
  integer, intent(in) :: nx, ny, nstep, flip_thr
  real, intent(in) :: azim_in(nx,ny,2*nstep+1)
  real, intent(inout) :: azim_out(nx,ny)
  ! Define variables to be used.
  integer naxis1, naxis2, naxis3
  real, allocatable :: cur_azim(:,:), hyp_azim(:,:)      ! intermediate vars
  real, allocatable :: mean_act_diff(:,:), mean_hyp_diff(:,:)   
  real, allocatable :: nbr_azim(:,:,:) 
  real, allocatable :: dazim_hyp(:,:,:), dazim_act(:,:,:) 
  integer j
  real pi, dtor

  pi=3.14159 ! accurate enough
  dtor = pi/180

  naxis1 = size(azim_in, 1)
  naxis2 = size(azim_in, 2)
  naxis3 = size(azim_in, 3)

  ! Alloc. intermediate variables 
  allocate( cur_azim( naxis1, naxis2)) 
  allocate( hyp_azim( naxis1, naxis2)) 
  allocate( mean_act_diff( naxis1, naxis2)) 
  allocate( mean_hyp_diff( naxis1, naxis2)) 
  allocate( nbr_azim( naxis1, naxis2, 2*nstep)) 
  allocate( dazim_hyp( naxis1, naxis2, 2*nstep)) 
  allocate( dazim_act( naxis1, naxis2, 2*nstep)) 

  cur_azim = azim_in(:,:,nstep+1)

  ! azimuths nsteps before & after current -- uses *updated* azim
  nbr_azim(:,:,1:nstep) = azim_in(:,:,1:nstep)
  nbr_azim(:,:,nstep+1:2*nstep) = azim_in(:,:,nstep+2:2*nstep+1)

  ! hypothetical, flipped azimuths  
  hyp_azim = mod(cur_azim + 180., 360.) ! azims w/flip

  ! Compute azimuth changes for actual & hypothetical azimuths
  ! (beware of branch cut at 0=360) 
  do j = 1,nstep 

     ! Diff. btwn. nbr & cur azimuths, 
     dazim_act(:,:,j     ) = & ! before nstep
          180/pi* &
          atan2(cos(cur_azim*dtor)*sin(nbr_azim(:,:,j     )*dtor) - &
          sin(cur_azim*dtor)*cos(nbr_azim(:,:,j     )*dtor), & 
          cos(cur_azim*dtor)*cos(nbr_azim(:,:,j     )*dtor) + &
          sin(cur_azim*dtor)*sin(nbr_azim(:,:,j     )*dtor) )
     
     dazim_act(:,:,j+nstep) = & ! after nstep
          180/pi* &
          atan2(cos(cur_azim*dtor)*sin(nbr_azim(:,:,j+nstep)*dtor) - &
          sin(cur_azim*dtor)*cos(nbr_azim(:,:,j+nstep)*dtor), &
          cos(cur_azim*dtor)*cos(nbr_azim(:,:,j+nstep)*dtor) + &
          sin(cur_azim*dtor)*sin(nbr_azim(:,:,j+nstep)*dtor) )
     
     ! Diff. btwn. nbr & hyp azimuths, 
     dazim_hyp(:,:,j     ) = & ! before nstep
          180/pi* &
          atan2(cos(hyp_azim*dtor)*sin(nbr_azim(:,:,j     )*dtor) - &
          sin(hyp_azim*dtor)*cos(nbr_azim(:,:,j     )*dtor), & 
          cos(hyp_azim*dtor)*cos(nbr_azim(:,:,j     )*dtor) + &
          sin(hyp_azim*dtor)*sin(nbr_azim(:,:,j     )*dtor) )

     dazim_hyp(:,:,j+nstep) = & ! after nstep
          180/pi* &
          atan2(cos(hyp_azim*dtor)*sin(nbr_azim(:,:,j+nstep)*dtor) - &
          sin(hyp_azim*dtor)*cos(nbr_azim(:,:,j+nstep)*dtor), &
          cos(hyp_azim*dtor)*cos(nbr_azim(:,:,j+nstep)*dtor) + &
          sin(hyp_azim*dtor)*sin(nbr_azim(:,:,j+nstep)*dtor) )
     
  end do

  ! Compare differences w/flips & w/o - Don't need to normalize
  mean_act_diff = sum(abs(dazim_act),3)/(2*nstep) 
  mean_hyp_diff = sum(abs(dazim_hyp),3)/(2*nstep)
  
  ! If flipping would lower mean diff, and azim change > thresh, then flip!
  where ((mean_hyp_diff .lt. mean_act_diff) .and. &
       (abs(dazim_act(:,:,nstep)) .gt. flip_thr)) cur_azim = hyp_azim 

  azim_out = cur_azim   ! store flipped data
   
  deallocate( cur_azim)
  deallocate( hyp_azim)
  deallocate( mean_act_diff)
  deallocate( mean_hyp_diff)
  deallocate( nbr_azim)
  deallocate( dazim_act)
  deallocate( dazim_hyp)


end subroutine azim_mindiff_jsoc
