!***********************************************************************************************************************************
subroutine minimise_energy(CalcE,CalcDE_reconfig)
!===================================================================================================================================
! This subroutine attempts to find the configuration of azimuthal angles that corresponds to the minimum of the "energy"
! summed over the region of interest. Simulated annealing is used to perform the search with a "temperature" that can vary from
! pixel to pixel. There are 3 main input parameters:
! tfac0 scales the initial temperature
! neq controls the number of reconfigurations that are attempted at each temperature setting
! tfactr is the cooling rate
!
! 2009,2010, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use anneal
   use bounds
   use energy_arrays
   use ran_pix
   use ranseed
   use spherical_deriv_coefficients
   use verbose

   implicit none

!
! Parameters tol_conv and nconv_min: If the relative changes in the energy is less than tol_conv for nconv_min consecutative
! temperature increments then the search has converged.
!
   real,parameter :: tol_conv=1e-5
   integer,parameter :: nconv_min=10
!
! Parameter tstop_par when the temperature is less than tstop_par times the initial temperature stop the search
!
   real,parameter :: tstop_par=1e-7

   integer :: i,j,idone,itemp,nconv,nsucc,nxnynz,neq_init,ii
   real :: E,E_prev,de,ran3,max_de,tstop,t
   real,dimension(:,:),allocatable :: tvar
   external :: CalcE,CalcDE_reconfig
!
! Initialisations
!
   allocate(tvar(nx,ny),DivB(nx,ny),Jz(nx,ny))
   allocate(ivec(((nx/jump)+1)*((ny/jump)+1)),jvec(((nx/jump)+1)*((ny/jump)+1)))

   call CalcE(E)

   nxny=dnx*dny
   nxnynz=nxny
   neq=neq*nxnynz
   neq_init=100*nxnynz

   nxjump=dnx/jump
   if (nxjump.lt.1) nxjump=1
   nyjump=dny/jump
   if (nyjump.lt.1) nyjump=1
   neq_init=neq_init/(nxjump*nyjump)
   neq=neq/(nxjump*nyjump)
   ia_prev=-1
   ja_prev=-1
!
! Calculate the initial temperature by sampling a very large number of reconfigurations (all accepted)
!
   do i=nxa,nxb
      do j=nya,nyb
         tvar(i,j)=0.
      enddo
   enddo
   max_de=0.
   do itemp=1,neq_init
!
! Loop over the sequence of pixels with random start point and separation: jump
!
      call get_ran_pix
      do ii=1,ng
         i=ivec(ii)
         j=jvec(ii)
!
! Calculate the change in the energy due the reconfiguration at the selected pixel
!
         call CalcDE(i,j,de,CalcDE_reconfig)
         if (abs(de).gt.max_de) max_de=abs(de)
         if (abs(de).gt.tvar(i,j)) tvar(i,j)=abs(de)
!
! Update the total energy (for the initial temperature accept all reconfigurations)
!
         E=E+de
!
! Implement the reconfiguration
!
         call reconfig(i,j)
      enddo
   enddo
!
! Scale the initial temperature according to input parameter tfac0
!
   do i=nxa,nxb
      do j=nya,nyb
         tvar(i,j)=tfac0*tvar(i,j)
      enddo
   enddo
   t=tfac0*max_de
!
! More initialisations
!
   tstop=tstop_par*t
   call CalcE(E)
   E_prev=E
   nconv=0
   idone=0
   if (iverb.eq.2) write(6,'(i9,3x,e15.8,3x,e15.8,3x,i3)')0,t*tstop_par/tstop,E,nconv
   do while (idone.eq.0)
      nsucc=0
      do itemp=1,neq
!
! Loop over the sequence of pixels with random start point and separation: jump
!
         call get_ran_pix
         do ii=1,ng
            i=ivec(ii)
            j=jvec(ii)
!
! Calculate the change in the energy due the reconfiguration at the selected pixel
!
            call CalcDE(i,j,de,CalcDE_reconfig)
!
! Decide whether or not to accept the reconfiguration
!
            if ((de.lt.0.).or.(ran3(seed).lt.exp(-de/tvar(i,j)))) then !metropolis criteria (simulated annealing, spatially dependent)
               nsucc=nsucc+1
!
! Update the energy
!
                E=E+de
!
! Implement the reconfiguration
!
               call reconfig(i,j)
            endif
         enddo
      enddo
!
! Reduce the temperature according to input parameter tfactr
!
      do i=nxa,nxb
         do j=nya,nyb
            tvar(i,j)=tfactr*tvar(i,j)
         enddo
      enddo
      t=t*tfactr
!
! Convergence test
!
      if (abs(E-E_prev)/(abs(E)+abs(E_prev)).lt.tol_conv) then
         nconv=nconv+1
      else
         nconv=0
      endif
      E_prev=E
      if (iverb.eq.2) write(6,'(i9,3x,e15.8,3x,e15.8,3x,i3)')nsucc,t*tstop_par/tstop,E,nconv
!
! Stopping conditions:
! 1. no reconfigurations were accepted during this temperature increment
! 2. temperature is "small"
! 3. the change in energy is "small" over several consecutive temperature steps (converged)
!
      if ((nsucc.eq.0).or.(t.lt.tstop).or.(nconv.ge.nconv_min)) idone=1
   enddo

   deallocate(tvar,DivB,Jz,ivec,jvec)
   if(allocated(ddt1)) deallocate(ddt1)
   if(allocated(ddt3)) deallocate(ddt3)
   if(allocated(ddp)) deallocate(ddp)
   if(allocated(ddr)) deallocate(ddr)

end subroutine minimise_energy
!***********************************************************************************************************************************
