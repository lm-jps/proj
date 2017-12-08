subroutine get_pils_rad(naxis1, naxis2, Brad, Bmag, pilmap, &
     thresh_brad, thresh_bmag, dilation_param)

!
! PURPOSE: Given input "Brad" array (Nx,Ny) corresponding to magnetogram 
!          pixels, this returns "pilmap" an (Nx, Ny) bitmap (of integer 
!          type, with default precision) for all pixels closer to 
!          opposite polarity than dilation_param (integer) 
!
! HISTORY: 2016/10/25, BT Welsch: get_pils.f90 started
!          2017/08/01, BTW: copied from get_pils, added bmag input
!
  integer*4, intent(in) :: naxis1, naxis2
  real*8, intent(in) :: brad(naxis1, naxis2) ! input magnetogram, B_los or B_r
  real*8, intent(in) :: bmag(naxis1, naxis2) ! input mag. vector field strength, |B|
  integer, intent(inout) :: pilmap(naxis1, naxis2) ! output map of PIL pixels (integer)
  integer, intent(in) :: dilation_param ! max. dist. btwn. +/- to be a PIL
  
  real*8 thresh_brad ! threshold in |B_r| below which pixel values are ignored
  real*8 thresh_bmag ! threshold in |B| below which pixel values are ignored

  integer, allocatable :: posmap(:,:),negmap(:,:) ! 2D
  integer, allocatable :: dilpos(:,:),dilneg(:,:) ! 2D

!  integer naxis1,naxis2
  integer alloc_count, dealloc_count ! bookkeeping of allocated arrays

  alloc_count = 0 
  dealloc_count = 0 

!  naxis1 = size(brad, 1)
!  naxis2 = size(brad, 2)

  allocate( posmap( naxis1, naxis2))
  allocate( negmap( naxis1, naxis2))
  allocate( dilpos( naxis1, naxis2))
  allocate( dilneg( naxis1, naxis2))
  alloc_count = alloc_count + 4

  where ((brad.ge.thresh_brad).and.(bmag.ge.thresh_bmag))
     posmap = 1
  elsewhere
     posmap = 0
  end where

  where ((brad.le.-thresh_brad).and.(bmag.ge.thresh_bmag))
     negmap = 1
  elsewhere
     negmap = 0
  end where

  call dilate(naxis1, naxis2, posmap, dilpos, dilation_param)
  call dilate(naxis1, naxis2, negmap, dilneg, dilation_param)

  !print *,dilpos ! for diagnostics of dilation subroutine

  pilmap = dilpos*dilneg

  !  Cleanup remaining allocated variables prior to exiting
  deallocate( posmap)
  deallocate( negmap)
  deallocate( dilpos)
  deallocate( dilneg)
  dealloc_count = dealloc_count + 4

  ! Check for matching allocations and deallocations.
  check_alloc: if (alloc_count /= dealloc_count) then
     print *,'Found allocate/ deallocate mismatch.'
  endif check_alloc

end subroutine get_pils_rad
