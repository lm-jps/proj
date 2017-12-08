subroutine pix2helio(naxis1, naxis2, lonlat, crpix, rsun_pix, crlt_obs, crota2, max_ang)
!
! PURPOSE: Given input "lonlat" array (Nx,Ny) corresponding to CCD pixels, along w/
!          pointing & solar disk params, this returns (longitude, latitude) for all
!          pixels closer to disk center than "max_ang" 
!
! HISTORY: 2016/10/18, BT Welsch: finished coding, sensible output for (latit, longit)
!
  integer*4, intent(in) :: naxis1, naxis2
  real*8, intent(inout) :: lonlat(naxis1, naxis2, 2) ! input gives size, output: longit, latit [radian]
  
  real*8 crpix(2) ! 2-element pix. address (x,y) of disk center
  real*8 rsun_pix ! radius of the sun in pixels
  real*8 crlt_obs ! Heliogrph Latit. of observ.'s disk ctr (B-angle)
  real*8 crota2   ! HMI uses (-) P-angle; so (+)=solar N west of Earth N [degrees]
  real*8 max_ang  ! input max. angular distance [degrees] from disk ctr (less than 90!)
  real*8 dpi, dtor ! double-precision pi and degrees-to-radians conversion

  real*8, allocatable :: x1d(:),y1d(:) ! 1D
  real*8, allocatable :: r2d(:,:),disktheta(:,:) ! 2D
  real*8, allocatable :: r3d(:,:,:) ! 3D
  integer i ! counting index for initializing x1d, y1d

!  integer naxis1,naxis2

  real*8 xhcs(3),yhcs(3),zhcs(3) ! 3-element unit vectors in heliographic coord. system

  integer alloc_count, dealloc_count ! bookkeeping of allocated arrays

  dpi = 4.D0*atan(1.D0)
  dtor = dpi/180.D0

  alloc_count = 0 
  dealloc_count = 0 

!  naxis1 = size(lonlat, 1)
!  naxis2 = size(lonlat, 2)

  allocate( x1d( naxis1))
  allocate( y1d( naxis2))
  allocate( r2d( naxis1, naxis2))
  allocate( disktheta( naxis1, naxis2))
  allocate( r3d( naxis1, naxis2, 3))   ! r3d = (x2d, y2d, z2d)

  alloc_count = alloc_count + 5

  x1d = (/ (i, i = 1, naxis1) /) 
  r3d(1:naxis1,1:naxis2,1) = spread(x1d, 2, naxis2) ! x2d

  y1d = (/ (i, i = 1, naxis2) /) 
  r3d(1:naxis1,1:naxis2,2)  = spread(y1d, 1, naxis1) ! y2d

  ! center coords on sun's center, scale to r_sun = 1.
  r3d(1:naxis1,1:naxis2,1) = (r3d(1:naxis1,1:naxis2,1)  - crpix(1))/ rsun_pix
  r3d(1:naxis1,1:naxis2,2) = (r3d(1:naxis1,1:naxis2,2)  - crpix(2))/ rsun_pix

  r2d = sqrt(r3d(1:naxis1,1:naxis2,1)**2 + r3d(1:naxis1,1:naxis2,2)**2 )

  ! For all (x,y) get "disk theta" (Cartesian --> polar), with P-angle removed
  where (r2d.gt.0) 
     disktheta = atan2(r3d(1:naxis1,1:naxis2,2), r3d(1:naxis1,1:naxis2,1)) + crota2*dpi/180.d0
  elsewhere
     disktheta = 0.
  end where

  ! Only get long./lat. for pix inside max_ang deg. from disk ctr.
  where (r2d.lt.sin(max_ang*dtor)) 
     r3d(1:naxis1,1:naxis2,1) = r2d*cos(disktheta)
     r3d(1:naxis1,1:naxis2,2) = r2d*sin(disktheta)
     r3d(1:naxis1,1:naxis2,3) = sqrt(1.d0 - r2d**2)
  elsewhere
     r3d(1:naxis1,1:naxis2,1) = 0.d0
     r3d(1:naxis1,1:naxis2,2) = 0.d0
     r3d(1:naxis1,1:naxis2,3) = 0.d0
  end where

  ! Get vectors for rotation to remove B-angle, in heliosph. coord syst.
  ! xhcs to right
  ! yhcs is up
  ! zhcs is toward observer
  xhcs(1:3) = (/ 1.d0,                0.d0,               0.d0  /) 
  yhcs(1:3) = (/ 0.d0,  cos(crlt_obs*dtor),  sin(crlt_obs*dtor) /) 
  zhcs(1:3) = (/ 0.d0, -sin(crlt_obs*dtor),  cos(crlt_obs*dtor) /) 

  !new x = (old x, y, z) # xhcs ! +x = toward West
  !new y = (old x, y, z) # yhcs ! +y = toward N. pole
  !new z = (old x, y, z) # zhcs ! +z = toward observer

  !longitude = atan( new x, new z)
  lonlat(1:naxis1,1:naxis2,1) = atan2( &
       r3d(1:naxis1,1:naxis2,1)*xhcs(1) + &  ! first arg is new x
       r3d(1:naxis1,1:naxis2,2)*xhcs(2) + &
       r3d(1:naxis1,1:naxis2,3)*xhcs(3) , &
       r3d(1:naxis1,1:naxis2,1)*zhcs(1) + &  ! second arg is new z
       r3d(1:naxis1,1:naxis2,2)*zhcs(2) + &
       r3d(1:naxis1,1:naxis2,3)*zhcs(3) )

  ! Check output:
  !print *, minval(lonlat(:,:,1))/dtor, maxval(lonlat(:,:,1))/dtor

  !latitude = 0.5*dpi - acos(new y) 
  lonlat(1:naxis1,1:naxis2,2) = 0.5*dpi - acos( &
       r3d(1:naxis1,1:naxis2,1)*yhcs(1) + &
       r3d(1:naxis1,1:naxis2,2)*yhcs(2) + &
       r3d(1:naxis1,1:naxis2,3)*yhcs(3) ) 

  ! Check output:
  !print *, minval(lonlat(:,:,2))/dtor, maxval(lonlat(:,:,2))/dtor

  !  Cleanup remaining allocated variables prior to exiting
  deallocate( x1d) 
  deallocate( y1d) 
  deallocate( r2d) 
  deallocate( disktheta) 
  deallocate( r3d) 
  dealloc_count = dealloc_count + 5

  ! Check for matching allocations and deallocations.
  check_alloc: if (alloc_count /= dealloc_count) then
     print *,'Found allocate/ deallocate mismatch.'
  endif check_alloc

end subroutine pix2helio
