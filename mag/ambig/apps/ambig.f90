!! This subroutine is the Fortran wrapper for the disambiguation code
!! by Graham Barnes. Takes necessary input parameters for the global
!! variables and then does the computation.
!!
!! Written by:	 Xudong Sun
!!
!! Description
!!		v1.0
!!		Modified from original ambig.f90, cut away the F wrapper
!!		Included functionality of rd_field.f90 and write_field.f90
!!		Dropped various flags for I/O
!!		Made geometry a variable instead of constant
!!
!! Version:
!!		v1.0	Jan 25 2010


subroutine ambig(&
		geometry_p, &			! set_geometry
		weak_p, &			! set_weak_field_treatment
		nx_p, ny_p, &			! sizes
		npad_p, &			! pad
		xcen_p, ycen_p, &		! disk_center
		verb_p, &			! verbose
		lambda_p, &			! WeightingFactor
		neq_p, tfactr_p, tfac0_p, &	! anneal
		seed_p, &			! ranseed
		ntx, nty, &			! ntiles
		Bx_p, By_p, Bz_p, &		! mgram_data
		probBa_p, &			! uncertainty
		bitmap_p, &     		! maskvec
		radius_p, &	        	! ephemeris
		nap)

!==============================================================!

! from original ambig.f90
   use set_geometry
   use sizes
   use maskvec
   use mgram_data
   use pad
   use point
   use pot_field
   use verbose
   use WeightingFactor
   use anneal
   use ranseed
! from rd_field.f90
   use bobs
   use bounds
   use constant
   use disk_center
   use ephemeris
   use pix_size
   use spherical_position

!==============================================================!

   implicit none

! set_geometry
   integer, intent(in) :: geometry_p
! set_weak_field_treatment
   integer, intent(in) :: weak_p 
! sizes
   integer, intent(in) :: nx_p, ny_p
! pad
   integer, intent(in) :: npad_p
! disk_center
   real, intent(in) :: xcen_p, ycen_p
! verbose
   integer, intent(in) :: verb_p
! WeightingFactor
   real, intent(in) :: lambda_p
! anneal
   integer, intent(in) :: neq_p
   real, intent(in) :: tfactr_p, tfac0_p
! ranseed
   integer, intent(in) :: seed_p
! ntiles
   integer, intent(in) :: ntx, nty
! mgram_data /maskvec
   real, dimension(nx_p,ny_p), intent(inout), target :: Bx_p, By_p, Bz_p
! maskvec
   integer, dimension(nx_p,ny_p), intent(inout), target :: bitmap_p
! ephemeris
   real, intent(in) :: radius_p
!
   integer, intent(in) :: nap
   real, dimension(nx_p,ny_p), intent(out), target :: probBa_p

! local temp variables
   integer :: i, j, nxp, nyp, ipflag
   real :: xx, yy, y2
   real,dimension(:,:),pointer :: probBa   ! Pointer to input array probBa_p
   external :: setup_OCBP_PF_dzh_4p, CalcE_OCBP_PF_dzh_4p, CalcDE_reconfig_OCBP_PF_dzh_4p
   external :: setup_spherical_PF_4p, CalcE_spherical_PF_4p, CalcDE_reconfig_spherical_PF_4p

!==============================================================!

! Assign pointers to input arrays
   Bx => Bx_p ; By => By_p ; Bz => Bz_p  ! mgram_data
   mask => bitmap_p                      ! maskvec
   probBa => probBa_p                    ! uncertainty

! sizes
   nx = nx_p
   ny = ny_p
! set_geometry
   geometry = geometry_p
! set_weak_field_treatment
   ipflag = weak_p
! npad
   npad = npad_p
! disk_center
! --> Compute position of center of FOV relative to disk center in same
! --> units as the radius/pixel size.
   xcen = xcen_p
   ycen = ycen_p
! verbose
   iverb = verb_p
! WeightingFactor
   lambda = lambda_p
! anneal
   neq = neq_p
   tfactr = tfactr_p
   tfac0 = tfac0_p
! ranseed
   seed = seed_p
! ephemeris
   radius = radius_p

! bobs
   allocate(bt(nx, ny))
   do j=1,ny
      do i=1,nx
         bt(i,j)=sqrt(Bx(i,j)**2+By(i,j)**2)
      enddo
   enddo

! allocate memory for vertical derivative
   allocate(dBzdz(nx,ny))
 
! --> Determine box enclosing all pixels contained in the mask.
 
   call boxit
   if(iverb.eq.1) write(*,*) 'boxit done'

! --> Ensure pixels along the edges of the box are excluded from the mask.

   do i=nxa,nxb
      mask(i,nya)=0
      mask(i,nyb)=0
   enddo
   do j=nya,nyb
      mask(nxa,j)=0
      mask(nxb,j)=0
   enddo

!======================================================================!

   select case (geometry)

!
! geometry=1: the divergence and vertical current density are
! calculated with the helioplanar approximation
!
      case(1)

! --> colat/lon of center of patch
         phi=asin(ycen/radius)
         theta=atan(xcen/sqrt(radius*radius-xcen*xcen-ycen*ycen))
write(*,*) 'xcen,ycen,radius',xcen,ycen,radius
write(*,*) 'phi,theta',phi,theta

! --> Construct coordinate transformation matrices.
         call transform(theta,phi)
         if(iverb.eq.1) write(*,*) 'transform calculated'
!
! --> Applying a Tukey windowing to smooth the transition to zero field
! --> at the edges. This reduces the impact of the periodic boundary
! --> conditions.
!
         call tukey(Bz,nx,ny,nap)
         if(iverb.eq.1) write(*,*) 'windowing done'
!
! --> Calculate transverse components of potential field, and vertical
! --> derivative of vertical component of potential field using FFTs.
!
         allocate(dBpzdz(nx,ny),Bpix(nx,ny),Bpiy(nx,ny))
         call potential(nx,ny)
         dBzdz=dBpzdz
         deallocate(dBpzdz,Bpix,Bpiy)
         if(iverb.eq.1) write(*,*) 'potential field done'
!
! --> Global minimization of "energy" to resolve ambiguity.
!
         call global(setup_OCBP_PF_dzh_4p, CalcE_OCBP_PF_dzh_4p, CalcDE_reconfig_OCBP_PF_dzh_4p)
         if(iverb.eq.1) write(*,*) 'optimization done'
!
! --> Perform acute angle to nearest neighbor ambiguity resolution on
! --> pixels below the noise threshold. 
!
         call nacute6
         if(iverb.eq.1) write(*,*) 'smoothing done'
!
! --> Temporary rough estimate of uncertainty in ambiguity resolution:
! --> 0.0 for pixels which were not disambiguated (outside mask)
! --> 0.6 for pixels which had acute angle method applied
! --> 0.9 for pixels which were above threshold and annealed
!
         do i=1,nx
            do j=1,ny
               if(mask(i,j).eq.0) then
                  probBa(i,j)=0.5
               else if(mask(i,j).eq.1) then
                  probBa(i,j)=0.6
               else if(mask(i,j).eq.2) then
                  probBa(i,j)=0.9
               endif
            enddo
         enddo
!
! geometry=2: the divergence and vertical current density are
! calculated in spherical coordinates
!
      case(2)
!
! --> Generate arrays of colatitude and longitude. 
         call colatlon

! --> Calculate potential field (and derivatives) using tiling.
         call tile(ntx,nty,nap,ipflag)
         if(iverb.eq.1) write(*,*) 'tiling done'

! --> Update mask to exclude pixels for which no potential field was calculated.
         do i=1,nx
            do j=1,ny
               mask(i,j)=mask(i,j)*tmask(i,j)
            enddo
         enddo

! --> Global minimization of "energy" to resolve ambiguity.
         call global(setup_spherical_PF_4p, CalcE_spherical_PF_4p, CalcDE_reconfig_spherical_PF_4p)
         if(iverb.eq.1) write(*,*) 'optimization done'
!
! --> Perform acute angle to nearest neighbor ambiguity resolution on
! --> pixels below the noise threshold. 
!
         call nacute6
         if(iverb.eq.1) write(*,*) 'smoothing done'
!
! --> Temporary rough estimate of uncertainty in ambiguity resolution:
! --> 0.0 for pixels which were not disambiguated (outside tmask)
! --> 0.5 for pixels which had weak field method applied
! --> 0.6 for pixels which had acute angle method applied
! --> 0.9 for pixels which were above threshold and annealed
!
         do i=1,nx
            do j=1,ny
               if(tmask(i,j).eq.0) then
                  probBa(i,j)=0.0
               elseif(mask(i,j).eq.0) then
                  probBa(i,j)=0.5
               else if(mask(i,j).eq.1) then
                  probBa(i,j)=0.6
               else if(mask(i,j).eq.2) then
                  probBa(i,j)=0.9
               endif
            enddo
         enddo

! --> Deallocate memory for spherical coordinates.

         deallocate(t,p,sint,cost,sinp,cosp)
         deallocate(tmask)
!
! geometry=other: error
!
      case default

         if(iverb.eq.1) write(*,*) 'wrong geometry'
         return

   end select

   deallocate(bt,dBzdz)

end subroutine ambig
