* -------------------------------------------------------------------------------------------
*
* o Solve time evolution of MHD (trans-Alfvenic) flow, for JSOC-pipeline and/or CGEM project
*
* o (Universal/solar) constants : cgs unit system
*    pi                          = 3.1415926535897932385E+00
*    Solar radius                = 6.96E+10[cm]
*    Gravitational constant      = 6.67259E-08
*    mass of proton              = 1.6726231E-24[g]
*    Boltzmann constant          = 1.380658E-16
*    solar mass                  = 1.989E+33[g]
*    Period of solar rotation (sidereal; viewed in the inertia)
*                                = 25.38 [day] (* Note this is approx.)
*    Carrington solar rot. period by definition
*                                = 27.275 [day] (*)
*
* o Frequently used numbers among arbitrary parameters.
*    Radius of the bottom of the domain
*                                = 1.010 solar radius
*    Specific heat ratio         = 1.05
*    Temperature at the critical points of Parker solution
*                                = 1 M[K]
*    Mass density at the inner boundary sphere
*                                = 1.67E-16 [g/cm^3] or 10^8[/cm^3]
*
* o Typical normalizing factors (with gamma = 1.0E+00)
*    Unit of velosity  v0  = sqrt(2kT/mp) = 1.285E+07[cm/s]
*    Unit of length    r0  = GM/2(v0^2)   = 4.019E+11[cm](=5.78R_sun)
*    Unit of time      t0  = r0/v0        = 3.128E+04[s]=8.72[h]=0.36[d]
*    Unit of density   n0  = arbitrary   (= 1.500E+10[count/cm^3])
*    Unit of pressure  p0  = (v0^2)*ro0   = 0.08255[erg/cm^3]
*    Unit of mag_field B0  = sqrt(4*pi*p0)= 1.019[Gauss]
*
** <History>
*
*  1) Feb. --- May, 2009
*     Modified to be called by C/Fortran Wrapper.
*     Base code file : q34a1a.nocommon_nonfix_wrapped4.for
*     A lot of tests, done in "daily-MHD" sim. for solar corona/wind.
*  2) April, 2009 -- June 2013, sporadically
*     Clean-up unused lines, and comments etc.
*
*                                                                                   K.Hayashi
*
* -------------------------------------------------------------------------------------------
*
       program mainlayer
       implicit none
       integer ii, jj, kk
       integer rte, rts, imagdata
!       write(*,*) 'Input II JJ KK'
!       read(*,*) ii, jj, kk
       ii = 72
       jj = 64
       kk = 128
       rts=1000
       rte=1000
       imagdata = 5 ! poly. term. Usually 5 or 7.
       call m3dq(ii,jj,kk,rte,rts,imagdata)
       stop
       endprogram
*
** -------------------------------------------
*
*       subroutine m3dq(ii,jj,kk)
       subroutine m3dq(ii,jj,kk,rte,rts,imagdata)
*       subroutine m3dq(ii,jj,kk,rte,rts,imagdata,pota,potb,
*     &                 ro7,pg7,tm7,ur7,ut7,up7,br7,bt7,bp7)
       implicit none
* number of grid for each three direction.
       intent(in) :: rte,rts,imagdata
       intent(in) :: ii,jj,kk
!       intent(in) :: pota,potb
!       intent(out):: ro7,pg7,tm7,ur7,ut7,up7,br7,bt7,bp7
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       real    pota(0:10,0:10), potb(0:10,0:10)
* variables for exporting. Currently, not used, but declared, hahaha.
!       real    ro7(-1:kk+2,-2:jj+2,-2:ii+3),pg7(-1:kk+2,-2:jj+2,-2:ii+3)
!       real    ur7(-1:kk+2,-2:jj+2,-2:ii+3),ut7(-1:kk+2,-2:jj+2,-2:ii+3)
!       real    up7(-1:kk+2,-2:jj+2,-2:ii+3),br7(-1:kk+2,-2:jj+2,-2:ii+3)
!       real    bt7(-1:kk+2,-2:jj+2,-2:ii+3),bp7(-1:kk+2,-2:jj+2,-2:ii+3)
!       real    tm7(-1:kk+2,-2:jj+2,-2:ii+3)
*
       integer rt               ! CR number
       integer imagdata         ! data type of obs magetic field
       integer iflow            ! hydrostatic or Parker
       integer irot             ! solar rotation
       integer iplasmap         ! coronal distribution
       real    mgfct            ! magnetic field strength at the North Pole
       real    dtfact           ! CFL num.
* choice of boundary treatment
       integer ichrsurf
       parameter(ichrsurf = 0) ! this may appear at other parts.
* counting etc.
       real    dt, dt0, mhdtime
       integer ncal, nref, nt, ntb
       integer iprc, iprc2
* save data
       integer nn, nums
       integer numnsave
       parameter(numnsave = 1000)
       integer nsave(numnsave)
* counting Carring num. etc.
       integer rts, rte, rt2
       integer ncomp, ne, n
       integer ncase
       parameter (ncase = 1000)
       integer rtstock0(ncase)  ! storage of CR Num.
* flags
       integer ihead, isvld, idummy
       logical lstop, lshow, lsymok
       logical lfind
* integer time factor
       real    tf2i ! time float to integer
* method
       integer irkstep ! def. Runge-Kutta step
       parameter(irkstep = 1)
* dB/dt by CT
!       logical lct
!       parameter(lct = .false.)
* magnetosonic wave
       integer ivmax, jvmax, kvmax, ivtrn, jvtrn, kvtrn
       integer ivalf, jvalf, kvalf, ivacs, jvacs, kvacs
       integer itmin, jtmin, ktmin
       real    vmax, valf, vtrn, vacs
* constant value(s) or normalizing factor(s)
       real    v0, r0, p0, b0, t0, rhoc, n0
       real    omega,  omega0
       real    gammax, gamma0
       real    vcrtrn, mcrtrn
* perturbation
       logical lpert, lpert0
       integer nstepert(2)
       integer cijkpert(3)
* number of cell averaged along longitude around axis
       integer kstep(0:jj)
* solar surface plasma map.
       real    denf0(0:jj,0:kk+1)
       real    tmpf0(0:jj,0:kk+1)
* inner boundary maps : suppose moving map in rest frame
!       logical lrestsys
!       parameter(lrestsys = .false.)
       real    shiftlng
       real    romap(0:jj,0:kk+1),pgmap(0:jj,0:kk+1)
       real    urmap(0:jj,0:kk+1),utmap(0:jj,0:kk+1)
       real    upmap(0:jj,0:kk+1),brmap(0:jj,0:kk+1)
       real    btmap(0:jj,0:kk+1),bpmap(0:jj,0:kk+1)
* working
       integer i, j, k
       character*48 flnsave
       integer kshift
       real    weir, weil
*
* constants
       real    pi
       parameter(pi = 3.1415926535897932385E+00)
*
* large arrays used somewhere in this program
*
* variables at previous step
       real    ro1(-2:jj+2,-1:kk+2,-2:ii+3),pg1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur1(-2:jj+2,-1:kk+2,-2:ii+3),ut1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up1(-2:jj+2,-1:kk+2,-2:ii+3),br1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt1(-2:jj+2,-1:kk+2,-2:ii+3),bp1(-2:jj+2,-1:kk+2,-2:ii+3)
* variables at the end of each step
       real    ro9(-2:jj+2,-1:kk+2,-2:ii+3),pg9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur9(-2:jj+2,-1:kk+2,-2:ii+3),ut9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up9(-2:jj+2,-1:kk+2,-2:ii+3),br9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt9(-2:jj+2,-1:kk+2,-2:ii+3),bp9(-2:jj+2,-1:kk+2,-2:ii+3)
* conservatives.
       real    ed1(-2:jj+2,-1:kk+2,-2:ii+3),mr1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    mt1(-2:jj+2,-1:kk+2,-2:ii+3),mp1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ed2(-2:jj+2,-1:kk+2,-2:ii+3),mr2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    mt2(-2:jj+2,-1:kk+2,-2:ii+3),mp2(-2:jj+2,-1:kk+2,-2:ii+3)
* variables at new step at hyperbolic stage
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
!       real    tm2(-2:jj+2,-1:kk+2,-2:ii+3)
* variables at new step at viscous (parabolic) stage.
       real    ro3(-2:jj+2,-1:kk+2,-2:ii+3),pg3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur3(-2:jj+2,-1:kk+2,-2:ii+3),ut3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up3(-2:jj+2,-1:kk+2,-2:ii+3),br3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt3(-2:jj+2,-1:kk+2,-2:ii+3),bp3(-2:jj+2,-1:kk+2,-2:ii+3)
* initial and boundary variables on the inner surface ! ... actually, only ??0(i,j,0) may be used.
       integer ipho
       parameter(ipho = 5)
       real    ro0(0:jj,1:kk,-1:ipho),pg0(0:jj,1:kk,-1:ipho)
       real    ur0(0:jj,1:kk,-1:ipho),ut0(0:jj,1:kk,-1:ipho)
       real    up0(0:jj,1:kk,-1:ipho),br0(0:jj,1:kk,-1:ipho)
       real    bt0(0:jj,1:kk,-1:ipho),bp0(0:jj,1:kk,-1:ipho)
* impose gradual increase of Br at/near surface
!       logical ldmagon
!       parameter(ldmagon = .false.)
       real    dbr0(0:jj,0:kk+1,-1:ipho)
       real    dbt0(0:jj,0:kk+1,-1:ipho)
       real    dbp0(0:jj,0:kk+1,-1:ipho)
* impose density etc.
!       real    dromap(0:jj,0:kk+1),dpgmap(0:jj,0:kk+1) ! differneces
!       real    durmap(0:jj,0:kk+1)
* dummy array....
       real    dumarray(-2:jj+2,-1:kk+2,-2:ii+3)
* coordinate system
       real    rr(-2:ii+3)
       real    theta(0:jj)
       real    sain(0:jj), kosa(0:jj)
       real    sainb(0:jj-1), kosab(0:jj-1)
       real    phi(-1:kk+2)
       real    sinphb(0:kk), cosphb(0:kk)
       real    sinph(0:kk+1), cosph(0:kk+1)
* Parker
       real    pgpark0(-2:ii+3), ropark0(-2:ii+3)
       real    pgparkb(-2:ii+2), roparkb(-2:ii+2)
* scalar for divB cleaner (Poisson solver)
       real    divbpsi1(-1:jj,0:kk+1,-1:ii+1)
       real    divbpsi2(-1:jj,0:kk+1,-1:ii+1)
       real    divmag(-1:jj,0:kk+1,-1:ii+1)
* scalar for Laplace Eq. solver
       real    lapsi1(0:jj,0:kk+1,-1:ii+1)
       real    lapsi2(0:jj,0:kk+1,-1:ii+1)
* Legendre
       real    lgndr(1:10,0:10,0:jj), dlgndr(1:10,0:10,0:jj) ! 10 is fixed.....

* ask indice
*       call askindice(rte,rts,imagdata,iflow,irot,iplasmap,lpert0)
       iflow = 2
       irot = 1
       iplasmap = 2
       lpert = .false.
!       if (lpert0) then
!         write(*,*) ' Input start and end step of perturbation in Nt'
!         read(*,*) (nstepert(n),n=1,2)
!         write(*,*) ' Perturbation Nt :',(nstepert(n),n=1,2)
!         write(*,*) ' Input central address i j k'
!         read(*,*) (cijkpert(n),n=1,3)
!       else
         nstepert(1) = 0
         nstepert(2) = 1
         cijkpert(1) = 0
         cijkpert(2) = 0
         cijkpert(3) = 0
!       endif

* find "POTrrrrX.dat"  or "STNrrrr.dat"
       do n = 1, ncase
         rtstock0(n) = -1
       enddo
       n = 0
       do rt = rts, rte
         lfind = .true.
         if (lfind) then
           n = n + 1
           rtstock0(n) = rt
           write(*,*) n, rtstock0(n)
         endif
       enddo
       ne = n
       write(*,'('' Num of Rot = '',i5)') ne

       do ncomp = 1, ncase ! ---------------------- start CR-loop

!         lpert = lpert0

         rt = rtstock0(ncomp)
         rt2 = rt
         if (rt .GE. 0) then
           write(flnsave,'(''stts'',i4.4,''.log'')') rt
           open(unit = 19, file = flnsave, status = 'unknown')

           write(*, '('' II JJ KK = '',3i5)') ii, jj, kk
           write(19,'('' II JJ KK = '',3i5)') ii, jj, kk

           write(*, '('' TimeStep = '',i5)') irkstep
           write(19,'('' TimeStep = '',i5)') irkstep

* read "NN.lst"
           do n = 1, numnsave
             nsave(n) = -1
           enddo
           flnsave = 'nn.lst'
           call setnn(nn, nref, nsave, nums, flnsave)
           write(*, '('' Nt_ref = '',i6)') nref
           write(*, '('' Nt_end = '',i6)') nn
           write(*, '('' Num. of Nt_save = '',i6)') nums
           write(19,'('' Nt_ref = '',i6)') nref
           write(19,'('' Nt_end = '',i6)') nn
           write(19,'('' Num. of Nt_save = '',i6)') nums

* set some important variables and show the values
*           dtfact = 0.25E+00 ! CFL
           dtfact = 0.40E+00 ! CFL, somewhat robust choice
*           dtfact = 0.60E+00 ! CFL, 2/3 is max for MUSCL

             gamma0 = 1.05E+00
             gammax = 1.05E+00
             mgfct  = 1.00E+00 ! in HMI synopotic map, all factors are corrected
!             mgfct  = 1.80E+00

           write(*, '('' DTfact= '',f11.6)') dtfact
           write(*, '('' GammaX= '',f11.6)') gammax
           write(*, '('' Gamma0= '',f11.6)') gamma0
           write(*, '('' Mgfct = '',f11.6)') mgfct
           write(19,'('' DTfact= '',f11.6)') dtfact
           write(19,'('' GammaX= '',f11.6)') gammax
           write(19,'('' Gamma0= '',f11.6)') gamma0
           write(19,'('' Mgfct = '',f11.6)') mgfct
           write(*,'(A)') '---'

*
* OpenMP parallel start........... any openMP directive should not apprear above.!!!
*
* First of all, all dimension arrays which will appear as "private"
*                 should be initalized here for full subscription ranges.

!$omp parallel do private(i,k,j)
             do i = -2, ii + 3
               if ((i .GE. -1) .AND. (i .LE. ipho)) then
                 do k = 1, kk
                 do j = 0, jj
                   ro0(j,k,i) = 0.0E+00
                   pg0(j,k,i) = 0.0E+00
                   ur0(j,k,i) = 0.0E+00
                   ut0(j,k,i) = 0.0E+00
                   up0(j,k,i) = 0.0E+00
                   br0(j,k,i) = 0.0E+00
                   bt0(j,k,i) = 0.0E+00
                   bp0(j,k,i) = 0.0E+00
                   dbr0(j,k,i) = 0.0E+00
                   dbt0(j,k,i) = 0.0E+00
                   dbp0(j,k,i) = 0.0E+00
                 enddo
                 enddo
               endif
             enddo
!$omp end parallel do
!$omp parallel do private(i,k,j)
             do i = -2, ii + 3
               do k = -1, kk + 2
               do j = -2, jj + 2
                 dumarray(j,k,i) = 0.0E+00
                 ro1(j,k,i) = 0.0E+00
                 pg1(j,k,i) = 0.0E+00
                 ur1(j,k,i) = 0.0E+00
                 ut1(j,k,i) = 0.0E+00
                 up1(j,k,i) = 0.0E+00
                 br1(j,k,i) = 0.0E+00
                 bt1(j,k,i) = 0.0E+00
                 bp1(j,k,i) = 0.0E+00
                 ro2(j,k,i) = 0.0E+00
                 pg2(j,k,i) = 0.0E+00
                 ur2(j,k,i) = 0.0E+00
                 ut2(j,k,i) = 0.0E+00
                 up2(j,k,i) = 0.0E+00
                 br2(j,k,i) = 0.0E+00
                 bt2(j,k,i) = 0.0E+00
                 bp2(j,k,i) = 0.0E+00
                 ro3(j,k,i) = 0.0E+00
                 pg3(j,k,i) = 0.0E+00
                 ur3(j,k,i) = 0.0E+00
                 ut3(j,k,i) = 0.0E+00
                 up3(j,k,i) = 0.0E+00
                 br3(j,k,i) = 0.0E+00
                 bt3(j,k,i) = 0.0E+00
                 bp3(j,k,i) = 0.0E+00
                 ro9(j,k,i) = 0.0E+00
                 pg9(j,k,i) = 0.0E+00
                 ur9(j,k,i) = 0.0E+00
                 ut9(j,k,i) = 0.0E+00
                 up9(j,k,i) = 0.0E+00
                 br9(j,k,i) = 0.0E+00
                 bt9(j,k,i) = 0.0E+00
                 bp9(j,k,i) = 0.0E+00
                 ed1(j,k,i) = 0.0E+00
                 mr1(j,k,i) = 0.0E+00
                 mt1(j,k,i) = 0.0E+00
                 mp1(j,k,i) = 0.0E+00
                 ed2(j,k,i) = 0.0E+00
                 mr2(j,k,i) = 0.0E+00
                 mt2(j,k,i) = 0.0E+00
                 mp2(j,k,i) = 0.0E+00
               enddo
               enddo
             enddo
!$omp end parallel do
!$omp parallel do private(i)
             do i = -2, ii + 3
               ropark0(i) = 0.0E+00
               pgpark0(i) = 0.0E+00
             enddo
!$omp end parallel do
!$omp parallel do private(i)
             do i = -2, ii + 3
               if ((i .GE. -2) .AND. (i .LE. ii + 2)) then
                 roparkb(i) = 0.0E+00
                 pgparkb(i) = 0.0E+00
               endif
             enddo
!$omp end parallel do
!$omp parallel do private(i,k,j)
             do i = -2, ii + 3
               if ((i .GE. -1) .AND. (i .LE. ii + 1)) then
                 do k = 0, kk + 1
                 do j = -1, jj
                   divbpsi1(j,k,i) = 0.0E+00
                   divbpsi2(j,k,i) = 0.0E+00
                   divmag(j,k,i) = 0.0E+00
                 enddo
                 enddo
               endif
             enddo
!$omp end parallel do
!$omp parallel do private(i,k,j)
             do i = -2, ii + 3
               if ((i .GE. -1) .AND. (i .LE. ii + 1)) then
                 do k = 0, kk + 1
                 do j = 0, jj
                   lapsi1(j,k,i) = 0.0E+00
                   lapsi2(j,k,i) = 0.0E+00
                 enddo
                 enddo
               endif
             enddo
!$omp end parallel do

!$omp parallel do private(i)
             do i = -2, ii + 3
               rr(i) = 0.0E+00
             enddo
!$omp end parallel do
             do j = 0, jj
               theta(j) = 0.0E+00
               sain(j)  = 0.0E+00
               kosa(j)  = 0.0E+00
             enddo
             do j = 0, jj - 1
               sainb(j) = 0.0E+00
               kosab(j) = 0.0E+00
             enddo
             do k =-1, kk + 2
               phi(k) = 0.0E+00
             enddo
             do k = 0, kk
               sinphb(k) = 0.0E+00
               cosphb(k) = 0.0E+00
             enddo
             do k = 0, kk + 1
               sinph(k) = 0.0E+00
               cosph(k) = 0.0E+00
             enddo

!$omp parallel do private(i,k,j)
             do i = -2, ii + 3
               if (i .EQ. 0) then
                 do  k = 0, kk + 1
                 do  j = 0, jj
                   romap(j,k) = 0.0E+00
                   pgmap(j,k) = 0.0E+00
                   urmap(j,k) = 0.0E+00
                   utmap(j,k) = 0.0E+00
                   upmap(j,k) = 0.0E+00
                   brmap(j,k) = 0.0E+00
                   btmap(j,k) = 0.0E+00
                   bpmap(j,k) = 0.0E+00
!                   dromap(j,k) = 0.0E+00
!                   dpgmap(j,k) = 0.0E+00
!                   durmap(j,k) = 0.0E+00
                 enddo
                 enddo
               endif
             enddo
!$omp end parallel do

* set initial value and coordinates...
           call setinit(ii,jj,kk,
     &                  pota, potb,
     &                  mgfct,rt2,imagdata,iflow,irot,iplasmap,
     &                  denf0,tmpf0,
     &                  omega, omega0, gamma0, gammax, vcrtrn, mcrtrn,
     &                  v0, r0, p0, b0, t0, rhoc, n0,
     &                  ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,dumarray,
     &                  rr,theta,phi,
     &                  sain,kosa,sainb,kosab,
     &                  sinph,cosph,sinphb,cosphb,
     &                  ropark0,pgpark0,roparkb,pgparkb,
     &                  lapsi1,lapsi2,lgndr,dlgndr)

           tf2i = t0/3.6E+01  ! 100 step = 1 hour

* set kstep
           call setkstep(ii,jj,kk,kstep,rr,sain)

           call axibnd(ii,jj,kk,
     &                 gammax,kstep,omega,
     &                 ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                 rr,sain,kosa,sinph,cosph)
           call lngbnd(ii,jj,kk,
     &                 ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)

!$omp parallel do private(i,k,j)
           do i = -1, ii + 2
             if (i .LE. ipho) then
               do  k = 1, kk
               do  j = 0, jj
                 ro0(j,k,i) = ro2(j,k,i)
                 pg0(j,k,i) = pg2(j,k,i)
                 ur0(j,k,i) = ur2(j,k,i)
                 ut0(j,k,i) = ut2(j,k,i)
                 up0(j,k,i) = up2(j,k,i)
                 br0(j,k,i) = br2(j,k,i)
                 bt0(j,k,i) = bt2(j,k,i)
                 bp0(j,k,i) = bp2(j,k,i)
               enddo
               enddo
             endif
           enddo
!$omp end parallel do
!$omp parallel do private(i,k,j)
           do i = -1, ii + 2
             if (i .EQ. 0) then
               do  k = 0, kk + 1
               do  j = 0, jj
                 romap(j,k) = ro2(j,k,i)
                 pgmap(j,k) = pg2(j,k,i)
                 urmap(j,k) = ur2(j,k,i)
                 utmap(j,k) = ut2(j,k,i)
                 upmap(j,k) = up2(j,k,i)
                 brmap(j,k) = br2(j,k,i)
                 btmap(j,k) = bt2(j,k,i)
                 bpmap(j,k) = bp2(j,k,i)
               enddo
               enddo
             endif
           enddo
!$omp end parallel do

* read plasma values if needed
           if (iflow .EQ. 0) then
             idummy = 0 ! read Plasma value, only. Mag_init will be left untouchted.
             call rdintvar(idummy,rt2,gamma0,
     &                     ii,jj,kk,
     &                     dumarray,ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
             call axibnd(ii,jj,kk,
     &                   gammax,kstep,omega,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                   rr,sain,kosa,sinph,cosph)
             call lngbnd(ii,jj,kk,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
           endif

           shiftlng = 0.0E+00
* spread kk range etc, and fix some geometrical requirement. this must be called after phi & theta a
           nt = 0
           ncal = 0
           iprc2 = 1
           call surcnd(ii,jj,kk,
     &                 lpert,nstepert,cijkpert,nt,
     &                 omega0,omega,shiftlng,
     &                 iprc2,gammax,gamma0,v0,n0,p0,r0,vcrtrn, ! MIND xx0(*,*,*) is referred.
     &                 ro0,pg0,ur0,ut0,up0,br0,bt0,bp0,
     &                 ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                 ed2,mr2,mt2,mp2,
     &                 rr,phi,sain,kosa,sainb,kosab,
     &                 sinph,cosph,ropark0,pgpark0)

           call axibnd(ii,jj,kk,
     &                 gammax,kstep,omega,
     &                 ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                 rr,sain,kosa,sinph,cosph)
           call lngbnd(ii,jj,kk,
     &                 ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
*
           ihead = 1
           isvld = 0 ! switch ... save var.
           call svldvar2(nt,rt2,gamma0,isvld,ihead,
     &                   ii,jj,kk,
     &                   dumarray,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)

           nt   = nref
           ntb  = nref
           if (nref .NE. 0) then
             isvld = 1 ! switch ... load var.
             ihead = 2 ! when idummy = 1, means nothing.
             call svldvar2(nref,rt2,gamma0,isvld,ihead,
     &                     ii,jj,kk,
     &                     dumarray,
     &                     ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
             call axibnd(ii,jj,kk,
     &                   gammax,kstep,omega,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                   rr,sain,kosa,sinph,cosph)

             call lngbnd(ii,jj,kk,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
             call surcnd(ii,jj,kk,
     &                   lpert,nstepert,cijkpert,nt,
     &                   omega0,omega,shiftlng,
     &                   iprc2,gammax,gamma0,v0,n0,p0,r0,vcrtrn,
     &                   ro0,pg0,ur0,ut0,up0,br0,bt0,bp0,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                   ed2,mr2,mt2,mp2,
     &                   rr,phi,sain,kosa,sainb,kosab,
     &                   sinph,cosph,ropark0,pgpark0)
             call axibnd(ii,jj,kk,
     &                   gammax,kstep,omega,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                   rr,sain,kosa,sinph,cosph)
             call lngbnd(ii,jj,kk,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
             nt = 99999
             isvld = 0 ! switch ... save var.
             ihead = 1
             call svldvar2(nt,rt2,gamma0,isvld,ihead,
     &                     ii,jj,kk,
     &                     dumarray,
     &                     ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
           endif

           lshow  = .false.
           lsymok = .true.

** do iteration
           nt   = nref
           ntb  = nref
           dt0 = 0.0E+00
           mhdtime = 0.0E+00
           lstop = .false.
           do while ((nt .LE. nn) .AND. (.NOT. lstop)) ! --- start of  Nstep-loop

             ncal = ncal + 1
!             if (lchkstp)  then
!               nt = ncal + nref
!             else
               nt = nref + int(mhdtime * tf2i + 1.0E-10)
!             endif

             call copyvmax(ii,jj,kk,
     &                     ncal,dt0,nt,gammax,
     &                     lstop,kstep,
     &                     ivmax,jvmax,kvmax,ivtrn,jvtrn,kvtrn,
     &                     ivalf,jvalf,kvalf,ivacs,jvacs,kvacs,
     &                     itmin,jtmin,ktmin,
     &                     vmax,valf,vtrn,vacs,
     &                     rr,sain,
     &                     ro1,pg1,ur1,ut1,up1,br1,bt1,bp1,
     &                     ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)

             if (lstop) then
               if (nt .EQ. nref) nt = nt + 1
               isvld = 0 ! switch ... save var.
               ihead = 1
               call svldvar2(nt,rt2,gamma0,isvld,ihead,
     &                       ii,jj,kk,
     &                       dumarray,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
               write(*,*) 'Some abnormal : Ncal = ', ncal
             else

               dt0 = dtfact * dt0           ! dt0 = 360.0E+00 / t0 ! fixed at 360 sec...

               do i = 1, nums ! i = 1, numnsave
                if ((nt  .GE. nsave(i)) .AND.
     &              (ntb .LT. nsave(i))) lshow = .true.
               enddo
               if ((lshow) .OR. (ncal .EQ. 1))
     &           call showstat(ncal,nn,nt,dt0,mhdtime,t0,v0,
     &                   ivmax, jvmax, kvmax, ivtrn, jvtrn, kvtrn,
     &                   ivalf, jvalf, kvalf, ivacs, jvacs, kvacs,
     &                   itmin, jtmin, ktmin,
     &                   vmax, valf, vtrn, vacs)

                 shiftlng = 0.0E+00 ! this may be used at posedmag() and posedpls()
                 weil = 1.0E+00
                 weir = 0.0E+00
                 kshift = 0

* R-K loop
               do iprc = 1, irkstep
                 iprc2 = iprc

                 dt = dt0 / float(irkstep - iprc + 1)

                 call prc1st(ii,jj,kk,
     &                       dt,iprc2,omega,gammax,kstep,vcrtrn,mcrtrn,
     &                       denf0,tmpf0,gamma0,
     &                       v0,r0,p0,b0,t0,rhoc,n0,ncal,lstop,
     &                       ro0,pg0,            br0,bt0,bp0,
     &                       ro1,pg1,ur1,ut1,up1,br1,bt1,bp1,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                           pg3,            br3,bt3,bp3,
     &                       ro9,pg9,ur9,ut9,up9,br9,bt9,bp9,
     &                       ed1,mr1,mt1,mp1,ed2,mr2,mt2,mp2,
     &                       rr,sain,kosa,sainb,kosab,
     &                       ropark0,pgpark0,roparkb,pgparkb,
     &                       sinph,cosph,sinphb,cosphb,dbr0,dbt0,dbp0)
                 if (lstop) then
                   if (nt .EQ. nref) nt = nt + 1
                   isvld = 0 ! switch ... save var.
                   ihead = 1
                   call svldvar2(nt,rt2,gamma0,isvld,ihead,
     &                           ii,jj,kk,
     &                           dumarray,
     &                           ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
                   write(*,*) 'Some abnormal : Ncal = ', ncal
                   stop
                 endif

                 call axibnd(ii,jj,kk,
     &                       gammax,kstep,omega,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       rr,sain,kosa,sinph,cosph)
                 call lngbnd(ii,jj,kk,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
                 call surcnd(ii,jj,kk,
     &                       lpert,nstepert,cijkpert,nt,
     &                       omega0,omega,shiftlng,
     &                       iprc2,gammax,gamma0,v0,n0,p0,r0,vcrtrn,
     &                       ro0,pg0,ur0,ut0,up0,br0,bt0,bp0,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       ed2,mr2,mt2,mp2,
     &                       rr,phi,sain,kosa,sainb,kosab,
     &                       sinph,cosph,ropark0,pgpark0)
                 call axibnd(ii,jj,kk,
     &                       gammax,kstep,omega,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       rr,sain,kosa,sinph,cosph)
                 call lngbnd(ii,jj,kk,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)

                 call divbcln(ii,jj,kk,
     &                        gammax,nt,ncal,dt,iprc2,kstep, ! divB-cleaner
     &                        pg2,ur2,ut2,up2,br2,bt2,bp2,br3,bt3,bp3,
     &                        rr,theta,phi,sain,kosa,sainb,
     &                        divbpsi1,divbpsi2,divmag)
                 call axibnd(ii,jj,kk,
     &                       gammax,kstep,omega,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       rr,sain,kosa,sinph,cosph)
                 call lngbnd(ii,jj,kk,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
                 call surcnd(ii,jj,kk,
     &                       lpert,nstepert,cijkpert,nt,
     &                       omega0,omega,shiftlng,
     &                       iprc2,gammax,gamma0,v0,n0,p0,r0,vcrtrn,
     &                       ro0,pg0,ur0,ut0,up0,br0,bt0,bp0,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       ed2,mr2,mt2,mp2,
     &                       rr,phi,sain,kosa,sainb,kosab,
     &                       sinph,cosph,ropark0,pgpark0)
                 call axibnd(ii,jj,kk,
     &                       gammax,kstep,omega,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       rr,sain,kosa,sinph,cosph)
                 call lngbnd(ii,jj,kk,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)

                 if (ichrsurf .EQ. 6) then ! in daily-run. no viscosity is added to polytropic cases
                   call prc2nd(ii,jj,kk,gammax,dt,kstep,ncal,
     &                       pg3,ur3,ut3,up3,
     &                       pg9,ur9,ut9,up9,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       rr,theta,phi,sain,     sainb,kosab)
                   call axibnd(ii,jj,kk,
     &                       gammax,kstep,omega,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       rr,sain,kosa,sinph,cosph)
                   call lngbnd(ii,jj,kk,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
                   call surcnd(ii,jj,kk,
     &                       lpert,nstepert,cijkpert,nt,
     &                       omega0,omega,shiftlng,
     &                       iprc2,gammax,gamma0,v0,n0,p0,r0,vcrtrn,
     &                       ro0,pg0,ur0,ut0,up0,br0,bt0,bp0,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       ed2,mr2,mt2,mp2,
     &                       rr,phi,sain,kosa,sainb,kosab,
     &                       sinph,cosph,ropark0,pgpark0)
                   call axibnd(ii,jj,kk,
     &                       gammax,kstep,omega,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       rr,sain,kosa,sinph,cosph)
                   call lngbnd(ii,jj,kk,
     &                       ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
                 endif


               enddo ! end of R-K loop

               if (lshow) then
                 ihead = 1
                 isvld = 0 ! switch ... save var.
                 call svldvar2(nt,rt2,gamma0,isvld,ihead,
     &                         ii,jj,kk,
     &                         dumarray,
     &                         ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
                 lshow = .false.
               endif

               if ((nt .EQ. nn) .AND . (.NOT. lshow)) lstop = .true.
               ntb = nt
             endif

             mhdtime = mhdtime + dt0

           enddo ! ------------------- End of nstep-do-while-loop

           close(19)

         endif

       enddo ! end of ncomp-loop

** added variable..
!       do i = -2, ii + 3
!       do k = -1, kk + 2
!       do j = -2, jj + 2
!         tm2(j,k,i) = pg2(j,k,i) / ro2(j,k,i)
!       enddo
!       enddo
!       enddo
** de-nondimensionalize and prepare export array
!       do i = -2, ii + 3
!       do j = -2, jj + 2
!       do k = -1, kk + 2
!         ro7(k,j,i) = ro2(jj-j,k,i) * n0 / 1.0E+08 ! in count e8 /cc
!         pg7(k,j,i) = pg2(jj-j,k,i) * p0
!         tm7(k,j,i) = tm2(jj-j,k,i)      ! in MK
!         ur7(k,j,i) = ur2(jj-j,k,i) * v0 / 1.0E+05 ! in km/s
!         ut7(k,j,i) = ut2(jj-j,k,i) * v0 / 1.0E+05 ! in km/s
!         up7(k,j,i) = up2(jj-j,k,i) * v0 / 1.0E+05 ! in km/s
!         br7(k,j,i) = br2(jj-j,k,i) * b0 ! in gauss
!         bt7(k,j,i) = bt2(jj-j,k,i) * b0 ! in gauss
!         bp7(k,j,i) = bp2(jj-j,k,i) * b0 ! in gauss
!       enddo
!       enddo
!       enddo

       endsubroutine
*       stop
*       endprogram

*
* -------------------------------------------------------------------------------------------
*
*                                  End of Main Program
*
* -------------------------------------------------------------------------------------------
*

*--------------------------------------------------------------
*
* Save and load data : note subscription order is i,j,k
*
*--------------------------------------------------------------
*
       subroutine svldvar2(nstep,rt,gamma0,isvld,ihead,
     &                     ii,jj,kk,
     &                     dumarray,
     &                     ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
       implicit none
       intent(in)    ::    nstep,rt,gamma0,isvld,ihead
       intent(in)    ::    ii,jj,kk
       intent(out)   ::    dumarray
       intent(inout) ::    ro2,pg2,ur2,ut2,up2,br2,bt2,bp2
* arguments
       integer nstep,rt,isvld,ihead
       real    gamma0
* large arrays
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       real    dumarray(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
* local
       integer i, j, k
       integer idummy1, jdummy1, kdummy1
       integer idummy2, jdummy2, kdummy2
       real    var8(8)
       character*48 flname
       logical ldummy
* save unformatted or ascii
       logical lunform
       parameter(lunform = .true.)

       if (isvld .EQ. 1) then ! load !
         write(flname,'(''d'',i6,''.'',i4)') 300000 + nstep, rt
         inquire(file = flname, exist = ldummy)
         if (ldummy) then
           open(unit=11,file=flname,status='old')
           do k = 0, kk + 1 ! MIND the range ... keep compati. with olds.
           do j = 0, jj
           do i = 0, ii + 1
             read(11,'(3i5,5(1x,e11.5),1x,3i5,3(1x,e11.5))')
     &         idummy1, jdummy1, kdummy1,
     &         var8(1), var8(2), var8(3), var8(4), var8(5),
     &         idummy2, jdummy2, kdummy2,
     &         var8(6), var8(7), var8(8)
             ro2(j,k,i) = var8(1)
             pg2(j,k,i) = var8(2)
             ur2(j,k,i) = var8(3)
             ut2(j,k,i) = var8(4)
             up2(j,k,i) = var8(5)
             br2(j,k,i) = var8(6)
             bt2(j,k,i) = var8(7)
             bp2(j,k,i) = var8(8)
           enddo
           enddo
           enddo
           close(11)
** unconvered part, .... needed ?
!$omp parallel do private(i,j)
           do i = 0, ii + 1
             do j = 0, jj
               ro2(j,-1,i) = ro2(j,kk-1,i)
               pg2(j,-1,i) = pg2(j,kk-1,i)
               ur2(j,-1,i) = ur2(j,kk-1,i)
               ut2(j,-1,i) = ut2(j,kk-1,i)
               up2(j,-1,i) = up2(j,kk-1,i)
               br2(j,-1,i) = br2(j,kk-1,i)
               bt2(j,-1,i) = bt2(j,kk-1,i)
               bp2(j,-1,i) = bp2(j,kk-1,i)
             enddo
             do j = 0, jj
               ro2(j,kk+2,i) = ro2(j,2,i)
               pg2(j,kk+2,i) = pg2(j,2,i)
               ur2(j,kk+2,i) = ur2(j,2,i)
               ut2(j,kk+2,i) = ut2(j,2,i)
               up2(j,kk+2,i) = up2(j,2,i)
               br2(j,kk+2,i) = br2(j,2,i)
               bt2(j,kk+2,i) = bt2(j,2,i)
               bp2(j,kk+2,i) = bp2(j,2,i)
             enddo
           enddo
!$omp end parallel do
         else
           write(flname,'(''x'',i6.6,''.'',i4)') nstep + 300000, rt
           inquire(file = flname, exist = ldummy)
           if (ldummy) then
             open(unit=11,file=flname,status='old',form='unformatted')
             read(11) ro2
             read(11) pg2
             read(11) ur2
             read(11) ut2
             read(11) up2
             read(11) br2
             read(11) bt2
             read(11) bp2
             close(11)
           else
             write(*,*) ' No file to be loaded ',flname
             write(*,*) ' No data file are found for ',rt,nstep
             stop
           endif
         endif

* adjust
!$omp parallel do private(i,k,j)
         do i = 0, ii + 1
           do k = 0, kk + 1
           do j = 0, jj
             pg2(j,k,i) = pg2(j,k,i) / gamma0
           enddo
           enddo
         enddo
!$omp end parallel do

       else ! save

* move global to local
!$omp parallel do private(i,k,j)
         do i = -2, ii + 3
           do k = -1, kk + 2
           do j = -2, jj + 2
             dumarray(j,k,i) = pg2(j,k,i) * gamma0 ! to make unit etc compati. with old prog. by K.H.
           enddo
           enddo
         enddo
!$omp end parallel do
         if (lunform) then
           if (ihead .EQ. 1) then
             write(flname,'(''x'',i6,''.'',i4)') nstep + 300000, rt
           else
             write(flname,'(''y'',i6,''.'',i4)') nstep + 300000, rt
           endif
           open(unit=2,file=flname,status='unknown',form='unformatted')
           write(2) ro2
           write(2) dumarray
           write(2) ur2
           write(2) ut2
           write(2) up2
           write(2) br2
           write(2) bt2
           write(2) bp2
           close(2)
         else
           if (ihead .EQ. 1) then
             write(flname,'(''d'',i6,''.'',i4)') nstep + 300000, rt
           else
             write(flname,'(''b'',i6,''.'',i4)') nstep + 300000, rt
           endif
           open(unit=2,file=flname,status='unknown')
           do k = 0, kk + 1
           do j = 0, jj
           do i = 0, ii + 1
             write(2,'(3i5,5(1x,e11.5),1x,3i5,3(1x,e11.5))')
     &         i+1000, j+2000, k+3000,
     &         ro2(j,k,i),dumarray(j,k,i),
     &         ur2(j,k,i),ut2(j,k,i),up2(j,k,i),
     &         i+1000, j+2000, k+3000,
     &         br2(j,k,i),bt2(j,k,i),bp2(j,k,i)
           enddo
           enddo
           enddo
           close(2)
         endif
       endif

       return
       endsubroutine

*
*--------------------------------------------------------------
*
* load initial plasma, etc... Not used any longer, maybe
*
*--------------------------------------------------------------
*
       subroutine rdintvar(iswitch,rt,gamma0,   ! iswitch = 0 for plasma otherwise for mag.
     &                     ii,jj,kk,
     &                     dumarray,ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
       implicit none
       intent(in)    ::    iswitch,rt,gamma0
       intent(in)    ::    ii,jj,kk
       intent(out)   ::    dumarray
       intent(inout) ::    ro2,pg2,ur2,ut2,up2,br2,bt2,bp2
* arguments
       integer rt, iswitch
       real    gamma0
* large arrays
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64) ! number of grid for each three direction.
       real    dumarray(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
* local
       integer i, j, k
       integer idummy1, jdummy1, kdummy1
       integer idummy2, jdummy2, kdummy2
       character*48 flname
       logical lfile1, lfile2
       real    rdummy1, rdummy2, rdummy3, rdummy4, rdummy5

       if (iswitch .EQ. 0) then ! load plasma
         write(flname,'(''mshplsf.'',i4.4)') rt
         inquire(exist = lfile1, file = flname)
         write(flname,'(''mshplsx.'',i4.4)') rt
         inquire(exist = lfile2, file = flname)
         if (lfile2) then
           write(flname,'(''mshplsx.'',i4.4)') rt
           open(unit=11,file=flname,status='old',form='unformatted')
           read(11) ro2
           read(11) pg2
           read(11) ur2
           read(11) ut2
           read(11) up2
           read(11) dumarray
           read(11) dumarray
           read(11) dumarray
           close(11)
         else if (lfile1) then
           write(flname,'(''mshplsf.'',i4.4)') rt
           open(unit=11,file=flname,status='old')
           do k = 0, kk + 1
           do j = 0, jj
           do i = 0, ii + 1
             read(11,'(3i5,5(1x,e11.5),1x,3i5,3(1x,e11.5))')
     &         idummy1, jdummy1, kdummy1,
     &         ro2(j,k,i),pg2(j,k,i),
     &         ur2(j,k,i),ut2(j,k,i),up2(j,k,i),
     &         idummy2, jdummy2, kdummy2,
     &         rdummy1, rdummy2, rdummy3
           enddo
           enddo
           enddo
           close(11)
         else
           write(*,*) ' No PlasmaInitData found '
           stop
         endif
         write(*,*) ' load .. done ' , flname
* move local to global
!$omp parallel do private(i,k,j)
         do i = 0, ii + 1
           do k = 0, kk + 1
           do j = 0, jj
             pg2(j,k,i) = pg2(j,k,i) / gamma0
           enddo
           enddo
         enddo
!$omp end parallel do

       else ! if not iswitch = 0 ................ load Mag

         write(flname,'(''mshmagf.'',i4.4)') rt
         inquire(exist = lfile1, file = flname)
         write(flname,'(''mshmagx.'',i4.4)') rt
         inquire(exist = lfile2, file = flname)
         if (lfile2) then
           write(flname,'(''mshmagx.'',i4.4)') rt
           open(unit=11,file=flname,status='old',form='unformatted')
           read(11) dumarray
           read(11) dumarray
           read(11) dumarray
           read(11) dumarray
           read(11) dumarray
           read(11) br2
           read(11) bt2
           read(11) bp2
           close(11)
         else if (lfile1) then
           write(flname,'(''mshmagf.'',i4.4)') rt
           open(unit=11,file=flname,status='old')
           do k = 0, kk + 1
           do j = 0, jj
           do i = 0, ii + 1
             read(11,'(3i5,5(1x,e11.5),1x,3i5,3(1x,e11.5))')
     &         idummy1, jdummy1, kdummy1,
     &         rdummy1, rdummy2, rdummy3, rdummy4, rdummy5,
     &         idummy2, jdummy2, kdummy2,
     &         br2(j,k,i),bt2(j,k,i),bp2(j,k,i)
           enddo
           enddo
           enddo
           close(11)
         else
           write(*,*) ' No MagInitData found '
           stop
         endif
         write(*,*) ' load .. done ' , flname
       endif

       return
       endsubroutine

*
*--------------------------------------------------------------
* Read and set the step number for starting and saving
* 1) Step number to start the calculation.
* 2) Step number to stop calculation.
* 3) Step number to save all the valiables.
*--------------------------------------------------------------
*
       subroutine setnn(nn,nref,nsave,nums,flnsave)
       implicit none
* arguments
       intent(out) ::   nn,nref,nsave,nums
       intent(in)  ::                      flnsave
       integer numnsave
       parameter(numnsave = 1000)
       integer nn, nref, nsave(numnsave), nums
       character*48 flnsave
* local
       integer i
       logical lexist

       inquire(file=flnsave,exist=lexist)
       if (lexist) then
         open(unit = 15, file = flnsave, status = 'old')
         read(15,*) nn
         read(15,*) nref
         i = 1
 100     continue
           read(15,*,END=99) nsave(i)
           i = i + 1
           goto 100
 99      continue
         close(15)
         nums = i - 1
       else
*         write(*,*) ' file '// flnsave // ' does not exist !!'
         nref =  0
**
         nn =    4000
         nsave(1) = 2000
         nsave(2) = 4000
         nums = 2

       endif

       return
       endsubroutine

*
*--------------------------------------------------------------
* Display some variable to the standard output device
*         (usually monitor or textfile)
*--------------------------------------------------------------
*
       subroutine showstat(ncal,nn,nt,dt,mhdtime,t0,v0,
     &                     ivmax, jvmax, kvmax, ivtrn, jvtrn, kvtrn,
     &                     ivalf, jvalf, kvalf, ivacs, jvacs, kvacs,
     &                     itmin, jtmin, ktmin,
     &                     vmax, valf, vtrn, vacs)
       implicit none
* arguments
       intent(in) ::       ncal,nn,nt,dt,mhdtime,t0,v0
       intent(in) ::       ivmax, jvmax, kvmax, ivtrn, jvtrn, kvtrn
       intent(in) ::       ivalf, jvalf, kvalf, ivacs, jvacs, kvacs
       intent(in) ::       itmin, jtmin, ktmin
       intent(in) ::       vmax, valf, vtrn, vacs
       integer ncal, nn, nt
       real    dt, mhdtime, t0, v0
       integer ivmax, jvmax, kvmax, ivtrn, jvtrn, kvtrn
       integer ivalf, jvalf, kvalf, ivacs, jvacs, kvacs
       integer itmin, jtmin, ktmin
       real    vmax, valf, vtrn, vacs
*
       write(*,'('' Nt='',i8,''/'',i8,'' Ncal='',i8,'' T='',e11.5,$)')
     &                 nt, nn, ncal,mhdtime
       write(*,'('' / '',e11.5,'' sec'')') mhdtime * t0
       write(*,'('' dT   = '',e11.5,'' / '',e11.5,'' sec  at'',3i4)')
     &                 dt,   dt * t0,             itmin, jtmin, ktmin
       write(*,'('' Vmax = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 vmax, vmax * v0 / 1.0E+05, ivmax, jvmax, kvmax
       write(*,'('' Valf = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 valf, valf * v0 / 1.0E+05, ivalf, jvalf, kvalf
       write(*,'('' Vtrn = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 vtrn, vtrn * v0 / 1.0E+05, ivtrn, jvtrn, kvtrn
       write(*,'('' Vacs = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 vacs, vacs * v0 / 1.0E+05, ivacs, jvacs, kvacs

       write(19,'('' Nt='',i8,''/'',i8,'' Ncal='',i8,'' T='',e11.5,$)')
     &                 nt, nn, ncal,mhdtime
       write(19,'('' / '',e11.5,'' sec'')') mhdtime * t0
       write(19,'('' dT   = '',e11.5,'' / '',e11.5,'' sec  at'',3i4)')
     &                 dt,   dt * t0,             itmin, jtmin, ktmin
       write(19,'('' Vmax = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 vmax, vmax * v0 / 1.0E+05, ivmax, jvmax, kvmax
       write(19,'('' Valf = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 valf, valf * v0 / 1.0E+05, ivalf, jvalf, kvalf
       write(19,'('' Vtrn = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 vtrn, vtrn * v0 / 1.0E+05, ivtrn, jvtrn, kvtrn
       write(19,'('' Vacs = '',e11.5,'' / '',e11.5,'' km/s at'',3i4)')
     &                 vacs, vacs * v0 / 1.0E+05, ivacs, jvacs, kvacs

       return
       end

*
* -----------------------------------------------------------------------
*
* o Set Coordinate System
* o Set Initial Value of Flow and Magnetic field
*
* Magnetic Field
*   1. Pseudo-Dipole Field.
*   2. Potential Field.
*   3. Pseudo-Quadrapole Field.
* Flow
*   1. Hydrostatic Equilibrium
*   2. Parker-Solution of Radial Flow
* -----------------------------------------------------------------------
*
       subroutine setinit(ii,jj,kk,
     &                    pota, potb,
     &                    mgfct,rt,imagdata,iflow,irot,iplasmap,
     &                    denf0,tmpf0,
     &                    omega,omega0,gamma0,gammax,vcrtrn,mcrtrn,
     &                    v0,r0,p0,b0,t0,rhoc,n0,
     &                    ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,dumarray,
     &                    rr,theta,phi,
     &                    sain,kosa,sainb,kosab,
     &                    sinph,cosph,sinphb,cosphb,
     &                    ropark0,pgpark0,roparkb,pgparkb,
     &                    lapsi1,lapsi2,lgndr,dlgndr)
       implicit none
* arguments
       intent(in)  ::     ii,jj,kk
       intent(in)  ::     mgfct,rt,imagdata,iflow,irot,iplasmap
       intent(in)  ::     denf0,tmpf0
       intent(in)  ::                  gamma0,gammax
       intent(out) ::     omega,omega0,              vcrtrn,mcrtrn
       intent(out) ::     v0,r0,p0,b0,t0,rhoc,n0
       intent(out) ::     ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,dumarray
       intent(out) ::     rr,theta,phi
       intent(out) ::     sain,kosa,sainb,kosab,sinph,cosph
       intent(out) ::     ropark0,pgpark0,roparkb,pgparkb
       intent(out) ::     lapsi1,lapsi2,lgndr,dlgndr
**
       real    pota(0:10,0:10), potb(0:10,0:10)
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64) ! number of grid for each three direction.
       real    mgfct
       integer rt, imagdata, iflow, irot, iplasmap
       real    omega, omega0, gamma0, gammax, vcrtrn, mcrtrn
       real    v0, r0, p0, b0, t0, rhoc, n0
       real    denf0(0:jj,0:kk+1) ! Currerntly this is not used.
       real    tmpf0(0:jj,0:kk+1) ! Left to host the density and temp. map from coronal base models.
* large arrays
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    rr(-2:ii+3),theta(0:jj),phi(-1:kk+2)
       real    sain(0:jj), kosa(0:jj)
       real    sainb(0:jj-1), kosab(0:jj-1)
       real    sinphb(0:kk), cosphb(0:kk)
       real    sinph(0:kk+1), cosph(0:kk+1)
       real    lapsi1(0:jj,0:kk+1,-1:ii+1)
       real    lapsi2(0:jj,0:kk+1,-1:ii+1)
       real    pgpark0(-2:ii+3), ropark0(-2:ii+3)
       real    pgparkb(-2:ii+2), roparkb(-2:ii+2)
       real    lgndr(1:10,0:10,0:jj), dlgndr(1:10,0:10,0:jj) ! 10 is given fixed. Maybe extendable.
       real    dumarray(-2:jj+2,-1:kk+2,-2:ii+3) ! dummy working array.
* local
       real    denf(0:jj,0:kk+1) ! this will not be used....
       real    tmpf(0:jj,0:kk+1)
       character*48 flnsave
       real    ro00, vr00, ro00surf
       real    tmpadjst
       integer i, j, k, idummy
       real    beta, rs
       integer nend
       real    tmp0, rbottom, rtop
       real    aa, bb, cc, dd, cosphl, sinphl
       real    r1, vpark, dpark, dprksurf
       real    massflux, bpole, rhosun
       logical lfile1, lfile2
       character*40 nameimag
* choice of whether system is in rest frame or rotating one
       logical lrestsys
       parameter(lrestsys = .false.) ! for coronal model running, this must be false
* density enhancement in accordance with magnetic field strength
!       real    mgstrc
!       parameter(mgstrc = 5.0E+00) ! in gauss
* flag for enforcing v//B condition at i=0,j,k, NOT for i = 0, j+1/2, k+1/2
       logical lvxbzero
       parameter(lvxbzero = .true.) ! maybe true
* universal const.
       real    gm, mh, kb,  pi
       real    rsun
       parameter(pi = 3.1415926535897932385E+00)
       parameter(rsun=6.96E+10)     ! Solar radius
       parameter(gm=6.67259E-08 * 1.989E+33) ! G * M_sun
       parameter(mh=1.6726231E-24)  ! Mass of Proton
       parameter(kb=1.380658E-16)   ! Boltzmann constant
* setting constant (1)
       real     tmpc
       parameter(tmpc=1.0E+06)       ! at r_c = 1.0
       real    tmpfact0, tmpfact
       parameter(tmpfact0 = 1.0E+00)
       real    rbotsr, rtopsr
       parameter(rbotsr = 1.01E+00, rtopsr = 50.0E+00) ! for coronal region
* setting constant (2)
       integer ichrsurf
       real    nsun0
       real    vcrtrn0
       parameter(ichrsurf = 0) ! this may appear at other parts.
       parameter(nsun0 =    2.00E+08) ! in count/cc
       parameter(vcrtrn0 =  5.00E+05) ! cm/s
* local (2)
       real    rri(-2:ii+3),thetai(0:jj),phii(-1:kk+2)
       real    saini(0:jj), kosai(0:jj)
       real    sainbi(0:jj-1), kosabi(0:jj-1)
       real    sinphbi(0:kk), cosphbi(0:kk)
       real    sinphi(0:kk+1), cosphi(0:kk+1)

* open log file
       write(flnsave,'(''init'',i4.4,''.log'')') rt ! flnsave = 'init.dat'
       open(unit = 18, file = flnsave,status = 'unknown')

* corona density chart
       if (iplasmap .EQ. 1) then
         do k = 0, kk + 1
         do j = 0, jj
           denf(j,k) = denf0(j,k)
           tmpf(j,k) = tmpf0(j,k)
         enddo
         enddo
       else
         do k = 0, kk + 1
         do j = 0, jj
           denf(j,k) = 1.0E+00 ! uniform plasma
           tmpf(j,k) = 1.0E+00 ! uniform plasma
         enddo
         enddo
       endif

* calculate normalizing variables (1)
       v0 = sqrt(2.0E+00 * kb * tmpc / mh * gamma0)
       r0 = gm / (2.0E+00 * v0**2)
       t0 = r0 / v0
       rbottom = rbotsr * rsun / r0 ! radius of lower boundary sphere
       rtop    = rtopsr * rsun / r0
       rs = rsun / r0

* set cordinates
       call setcrdnt(ii,jj,kk,
     &               rbottom,rtop,rri,thetai,phii,
     &               saini,kosai,sainbi,kosabi,
     &               sinphi,cosphi,sinphbi,cosphbi)
!$omp parallel do private(i)
       do i = -2, ii + 3
         rr(i) = rri(i)
       enddo
!$omp end parallel do
       do j = 0, jj
         theta(j) = thetai(j)
         sain(j)  = saini(j)
         kosa(j)  = kosai(j)
       enddo
       do j = 0, jj - 1
         sainb(j)  = sainbi(j)
         kosab(j)  = kosabi(j)
       enddo
       do k =-1, kk + 2
         phi(k) = phii(k)
       enddo
       do k = 0, kk
         sinphb(k) = sinphbi(k)
         cosphb(k) = cosphbi(k)
       enddo
       do k = 0, kk + 1
         sinph(k) = sinphi(k)
         cosph(k) = cosphi(k)
       enddo

       write(*,'('' r : '',f9.4,'' =< r =< '',2f9.4)')rri(0),rri(ii+1)
*
       do i = 0, ii + 1
         write(18,*) i+1000, rri(i), rri(i) / rs
       enddo
       do j = 0, jj
         write(18,*) j+2000, thetai(j),9.0E+01-thetai(j)/pi*1.8E+02
       enddo
       do k = 0, kk + 1
         write(18,*) k+3000, phii(k), phii(k) / pi * 1.8E+02
       enddo
       write(18,*) 0, 0.0, 0.0
       write(*,*) 'SET_COORDINATE : completed'
*       write(*,'(A)') '---'

* solar rotation
       if (irot .EQ. 0) then
         omega  = 0.0E+00
         omega0 = 0.0E+00
         write(*,*)  'Rotation NOT Included'
         write(19,*) 'Rotation NOT Included'
       else
         omega0 = t0 * 2.0E+00 * pi
     &         / (25.3E+00 * float(24 * 60 * 60))
         omega0 = omega0 * float(irot)
         if (lrestsys) then
           omega = 0.0E+00
         else
           omega = omega0
         endif
         write(*,*)  'Rotation Included'
         write(19,*) 'Rotation Included'
       endif

* restore Parker solutions at rr(i) and rr(i+1/2)

!$omp parallel do private(i,r1,tmpfact,vpark,dpark,aa)
       do i = -2, ii + 3
         tmpfact = 1.0E+00 ! Tc = 1 MK
         aa = gamma0
         r1 = rr(i)
         call getpark(tmpfact,r1,vpark,dpark,aa)
         ropark0(i) = dpark
         pgpark0(i) = dpark**gamma0 ! MIND gamma
       enddo
!$omp end parallel do
!$omp parallel do private(i,r1,tmpfact,vpark,dpark,aa)
       do i = -2, ii + 2
         tmpfact = 1.0E+00
         aa = gamma0
         r1 = (rr(i) + rr(i+1)) * 0.5E+00
         call getpark(tmpfact,r1,vpark,dpark,aa)
         roparkb(i) = dpark
         pgparkb(i) = dpark**gamma0 ! MIND gamma
       enddo
!$omp end parallel do

* calculate normalizing variables (2)
       r1 = rr(0)
       tmpfact = tmpfact0
       aa = gamma0
       call getpark(tmpfact,r1,vpark,dprksurf,aa)
       rhosun = nsun0 * mh ! surface mass density
       massflux =  rhosun * (vpark * v0)
       rhoc   = rhosun / dprksurf  ! normalized
       p0 = v0**2 * rhoc ! / gamma0
       b0 = sqrt(4.0E+00 * pi * rhoc * v0**2)
       n0 = rhoc / mh
* temperature at surface if T_c = 1 MK, in unit of 1 MK
       tmpadjst = 1.0E+00 / dprksurf**(gamma0 - 1.0E+00)

       if ((ichrsurf .EQ. 0) .OR. (ichrsurf .EQ. -1)) then
         vcrtrn = 1.0E+12 / v0 ! maybe do nothing
       else
         vcrtrn = vcrtrn0 / v0 ! should be same in procst()
       endif
       mcrtrn = rhosun / rhoc * vcrtrn

         if (iflow .EQ. 1) then
           write(*,*)  'Set flow as Hydrostatic'
           write(19,*) 'Set flow as Hydrostatic'

!$omp parallel do private(i,k,j,aa,ro00,vr00)
           do i = -2, ii + 3
             aa = 2.0E+00 * (1.0E+00 / rr(i) - 1.0E+00)
             ro00 = exp(aa)
     &            * 0.3E+00 ! tekitou !!
             vr00 = 0.0E+00
             do k = 0, kk + 1
             do j = 0, jj
               ro2(j,k,i) = ro00
               pg2(j,k,i) = ro00**gamma0 / gamma0 ! adjust !
               ur2(j,k,i) = vr00
               ut2(j,k,i) = 0.0E+00
                 up2(j,k,i) = - omega*sain(j)*(rr(i) - rr(0))
                 if (i .EQ. 0) up2(j,k,i) = 0.0E+00  ! confirm
             enddo
             enddo
           enddo
!$omp end parallel do
         else if (iflow .EQ. -1) then ! test
           write(*,*)  'Set flow linearly'
           write(19,*) 'Set flow linearly'
           aa = 1.0E+06 / v0 ! [cm/s] at inner bound.
           bb = 3.0E+07 / v0 ! [cm/s] at outer bound.
           dd = aa ! i = 0
           ro00surf = massflux / (dd * v0) / rhoc
!$omp parallel do private(i,k,j,cc,dd,ro00,vr00)
           do i = -2, ii + 3
             cc = (rr(i) - rr(0)) / (rr(ii+1) - rr(0))
             cc = cc * 2.0E+00
             cc = cc**2 / (1.0E+00 + cc**2)
             dd = (1.0E+06 / v0) *(1.0E+00 - cc) ! near Sun
     &          + (3.0E+07 / v0) * cc
             vr00 = dd
             ro00 = massflux / (dd * v0) / (rr(i)/rr(0))**2 / rhoc
             do k = 0, kk + 1
             do j = 0, jj
               ro2(j,k,i) = ro00
               pg2(j,k,i) = ro00 * (ro00/ro00surf)**gamma0 /gamma0
               ur2(j,k,i) = vr00
               ut2(j,k,i) = 0.0E+00
                 up2(j,k,i) = - omega*sain(j)*(rr(i) - rr(0))
                 if (i .EQ. 0) up2(j,k,i) = 0.0E+00  ! confirmation
             enddo
             enddo
           enddo
!$omp end parallel do
         else
           write(*,*)  'Set flow as Parker-Solution'
           write(19,*) 'Set flow as Parker-Solution'
           if (iplasmap .NE. 1) then
!$omp parallel do private(i,k,j,r1,aa,bb)
!$omp&            private(tmpfact,ro00,vr00,vpark,dpark)
             do i = -2, ii + 3
!               if (ldistant) then
!               else
                 aa = gamma0
                 tmpfact = tmpfact0
                 r1= rr(i)
                 call getpark(tmpfact,r1,vpark,dpark,aa)
                 bb = 1.0E+00
                 ro00 = dpark
                 vr00 = vpark
                 do k = 0, kk + 1
                 do j = 0, jj
                   ro2(j,k,i) = bb * ro00
                   pg2(j,k,i) = bb * ro00**gamma0 / gamma0 * tmpfact
                   ur2(j,k,i) = vr00
                   ut2(j,k,i) = 0.0E+00
                   up2(j,k,i) = - omega*sain(j)*(rr(i) - rr(0))
                   if (i .EQ. 0) up2(j,k,i) = 0.0E+00  ! confirmation
                 enddo
                 enddo
!               endif
             enddo
!$omp end parallel do
           else
!$omp parallel do private(i,k,j,tmpfact,r1,vpark,dpark,ro00,vr00,aa)
             do i = -2, ii + 3
               do k = 0, kk + 1
               do j = 0, jj
                 aa = gamma0
                 tmpfact = tmpfact0
                 r1= rr(i)
                 call getpark(tmpfact,r1,vpark,dpark,aa)
                 tmpfact = tmpf(j,k) * tmpadjst
                 ro2(j,k,i) = dpark * denf(j,k)
                 pg2(j,k,i) = ro2(j,k,i)**gamma0 / gamma0 * tmpfact
                 ur2(j,k,i) = vpark
                 ut2(j,k,i) = 0.0E+00
                   up2(j,k,i) = - omega*sain(j)*(rr(i) - rr(0))
                   if (i .EQ. 0) up2(j,k,i) = 0.0E+00  ! confirmation
               enddo
               enddo
             enddo
!$omp end parallel do

           endif
         endif  ! end if (iflow .EQ. ......)

       tmp0 = pg2(0,0,0) / ro2(0,0,0) * tmpc * gamma0
       write(*, '('' Temperature at North Pole = '',e11.5)') tmp0
       write(19,'('' Temperature at North Pole = '',e11.5)') tmp0
       write(18,'('' Temperature at North Pole = '',e11.5)') tmp0
       write(*,'(A)') '---'
       rhosun = ro2(0,0,0) * n0 * mh
       write(*, '('' Density at North Pole = '',e11.5)') rhosun
       write(19,'('' Density at North Pole = '',e11.5)') rhosun
       write(18,'('' Density at North Pole = '',e11.5)') rhosun
       write(*,'(A)') '---'
       write(*, '('' Num Dens. at North Pole = '',e11.5)') rhosun/mh
       write(19,'('' Num Dens. at North Pole = '',e11.5)') rhosun/mh
       write(18,'('' Num Dens. at North Pole = '',e11.5)') rhosun/mh
       write(*,'(A)') '---'

* restore Legendre.
       call setlgndr(jj,imagdata,sain,kosa,lgndr,dlgndr)

* calculate potential field in [micro T] or load data
       if (imagdata .LE. 0) then
         write(nameimag,'(''mshmagf.'',i4.4)') rt
         inquire(exist = lfile1, file = nameimag)
         write(nameimag,'(''mshmagx.'',i4.4)') rt
         inquire(exist = lfile2, file = nameimag)
         if (lfile1 .OR. lfile2) then
           idummy  = 1 ! switch to load Mag.
           call rdintvar(idummy,rt,gamma0, ! load plasma/mag
     &                   ii,jj,kk,
     &                   dumarray,ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp parallel do private(i,k,j)
           do i = -2, ii + 3
             do k = 0, kk + 1
             do j = 0, jj
               br2(j,k,i) = br2(j,k,i) * mgfct
               bt2(j,k,i) = bt2(j,k,i) * mgfct
               bp2(j,k,i) = bp2(j,k,i) * mgfct
             enddo
             enddo
           enddo
!$omp end parallel do
         else
           call mshpotmg(ii,jj,kk,rt,
     &                   br2,bt2,bp2,
     &                   rr,theta,phi,sain,lapsi1,lapsi2)
!$omp parallel do private(i,k,j)
           do i = -2, ii + 3
             do k = 0, kk + 1
             do j = 0, jj
               br2(j,k,i) = br2(j,k,i) / 1.0E+02 / b0 * mgfct ! conver mT => normalized unit.
               bt2(j,k,i) = bt2(j,k,i) / 1.0E+02 / b0 * mgfct
               bp2(j,k,i) = bp2(j,k,i) / 1.0E+02 / b0 * mgfct
             enddo
             enddo
           enddo
!$omp end parallel do
         endif
       else
         nend = imagdata
         call expndmag(ii,jj,kk,
     &                 pota, potb,
     &                 nend,rt,rs,
     &                 br2,bt2,bp2,
     &                 rr,phi,sain,kosa,lgndr,dlgndr)
!$omp parallel do private(i,k,j)
         do i = -2, ii + 3
           do k = 0, kk + 1
           do j = 0, jj
             br2(j,k,i) = br2(j,k,i) / 1.0E+02 / b0 * mgfct ! conver mT => normalized unit.
             bt2(j,k,i) = bt2(j,k,i) / 1.0E+02 / b0 * mgfct
             bp2(j,k,i) = bp2(j,k,i) / 1.0E+02 / b0 * mgfct
           enddo
           enddo
         enddo
!$omp end parallel do
       endif


       write(*,*) 'SET_MAG : completed'
       bpole = br2(0,0,0) * b0
       write(*, '('' Br(0,0,0) => '',e11.5,''[gauss]'')') bpole
       write(18,'('' Br(0,0,0) => '',e11.5,''[gauss]'')') bpole
       write(19,'('' Br(0,0,0) => '',e11.5,''[gauss]'')') bpole
       beta = br2(0,0,0)**2 + bt2(0,0,0)**2 + bp2(0,0,0)**2
       beta = beta * 0.5E+00 / pg2(0,0,0)
       write(*, '('' Beta^(-1) at North pole = '',e11.5)') beta
       write(18,'('' Beta^(-1) at North pole = '',e11.5)') beta
       write(19,'('' Beta^(-1) at North pole = '',e11.5)') beta

** make v // b at i = 0, MIND this is very special only for initial value setting.
       if (lvxbzero) then
!$omp parallel do private(i,k,j,aa,bb)
         do i = -2, ii + 3
           if (i .EQ. 0) then
             do k = 0, kk + 1
             do j = 0, jj
               bb = br2(j,k,i)**2 + bt2(j,k,i)**2 + bp2(j,k,i)**2
               if (bb .GT. 1.0E-10) then ! if bb = 0 such as HCS, parker solution may used.....
                 aa = ur2(j,k,i) * br2(j,k,i)
                 ur2(j,k,i) = aa * br2(j,k,i) / bb
                 ut2(j,k,i) = aa * bt2(j,k,i) / bb
                 up2(j,k,i) = aa * bp2(j,k,i) / bb
               endif
             enddo
             enddo
           endif
         enddo ! end i-loop
!$omp end parallel do
       endif

** confirm geometrical requirement(s) at the rotational axis
      if (kk .LE. 3) then
!$omp parallel do private(i,k)
         do i = -2, ii + 3
           do k = 0, kk + 1
             ut2( 0,k,i) = 0.0E+00
             ut2(jj,k,i) = 0.0E+00
             up2( 0,k,i) = 0.0E+00
             up2(jj,k,i) = 0.0E+00
             bt2( 0,k,i) = 0.0E+00
             bt2(jj,k,i) = 0.0E+00
             bp2( 0,k,i) = 0.0E+00
             bp2(jj,k,i) = 0.0E+00
           enddo
         enddo
!$omp end parallel do

       else

!$omp parallel do private(i,k,cosphl,sinphl,aa,bb,cc,dd)
           do i = -2, ii + 3
             aa = 0.0E+00 ! Bt at North
             bb = 0.0E+00 ! Bt at South
             cc = 0.0E+00 ! Bp
             dd = 0.0E+00 ! Bp
             do k = 1, kk
               cosphl = cosph(k)
               sinphl = sinph(k)
               aa = aa + (bt2(   1,k,i)*cosphl
     &                   -bp2(   1,k,i)*sinphl)
               bb = bb - (bt2(jj-1,k,i)*cosphl
     &                   -bp2(jj-1,k,i)*sinphl)
               cc = cc + (bt2(   1,k,i)*sinphl
     &                   +bp2(   1,k,i)*cosphl)
               dd = dd - (bt2(jj-1,k,i)*sinphl
     &                   -up2(jj-1,k,i)*cosphl)
             enddo ! end k-loop
             aa = aa / float(kk)
             bb = bb / float(kk)
             cc = cc / float(kk)
             dd = dd / float(kk)
             do k = 1, kk
               cosphl = cosph(k)
               sinphl = sinph(k)
               ut2( 0,k,i) = 0.0E+00
               ut2(jj,k,i) = 0.0E+00
               up2( 0,k,i) = 0.0E+00
               up2(jj,k,i) = 0.0E+00
               bt2( 0,k,i) =  aa * cosphl + cc * sinphl
               bt2(jj,k,i) = -bb * cosphl - dd * sinphl
               bp2( 0,k,i) = -aa * sinphl + cc * cosphl
               bp2(jj,k,i) = -bb * sinphl + dd * cosphl
             enddo
           enddo ! end of i-loop
!$omp end parallel do

       endif ! end if (kk .LT. 3)

** show variables
       write(18,'(A)') ' PARAMETERs.... '
       write(18,'('' rho_bot   = '',e11.5,''[g/cm^3]  '')')rhosun
       write(18,'('' n_bot     = '',e11.5,''[/cm^3]  '')') rhosun/mh
       write(18,'('' B_pole    = '',e11.5,''[gauss]   '')') bpole
       write(18,'('' omega0    = '',e11.5,''[deg/day] '')')
     &            omega0/ t0 * 3.6E+03 * 2.4E+01 * 1.80E+02 / pi
       write(18,'('' omega     = '',e11.5,''[deg/day] '')')
     &            omega / t0 * 3.6E+03 * 2.4E+01 * 1.80E+02 / pi
       write(18,'('' T_c       = '',e11.5,''[K]       '')') tmpc
       write(18,'('' gammax    = '',e11.5)') gammax
       write(18,'('' gamma     = '',e11.5)') gamma0
       write(18,'('' gamma0    = '',e11.5)') gamma0
       write(18,'('' v0        = '',e11.5,''[cm/sec]  '')') v0
       write(18,'('' r0        = '',e11.5,''[cm]      '')') r0
       write(18,'('' t0        = '',e11.5,''[sec]     '')') t0
       write(18,'('' p0        = '',e11.5,''[erg/cm^3]'')') p0/gamma0 ! keep compati. with old prog.
       write(18,'('' b0        = '',e11.5,''[gauss]   '')') b0
       write(18,'('' rhoc      = '',e11.5,''[g/cm^3]'')') rhoc
       write(18,'('' n0        = '',e11.5,''[/cm^3]'')') n0
       write(18,'('' r_top     = '',e11.5,''[r_sun]'')') rri(ii+1)/rs
       write(18,'('' r_bottom  = '',e11.5,''[r_sun]'')') rri(0   )/rs
**
       write(*,'(A)') '---'
       write(*,'(A)') ' PARAMETERs.... '
       write(*,'('' rho_bot   = '',e11.5,''[g/cm^3]  '')') rhosun
       write(*,'('' n_bot     = '',e11.5,''[/cm^3]  '')') rhosun/mh
       write(*,'('' B_pole    = '',e11.5,''[gauss]   '')') bpole
       write(*,'('' omega     = '',e11.5,''[deg/day] '')')
     &            omega / t0 * 3.6E+03 * 2.4E+01 * 1.80E+02 / pi
       write(*,'('' omega0    = '',e11.5,''[deg/day] '')')
     &            omega0/ t0 * 3.6E+03 * 2.4E+01 * 1.80E+02 / pi
       write(*,'('' T_c       = '',e11.5,''[K]       '')') tmpc
       write(*,'('' gammax    = '',e11.5)') gammax
       write(*,'('' gamma     = '',e11.5)') gamma0
       write(*,'('' gamma0    = '',e11.5)') gamma0
       write(*,'('' v0        = '',e11.5,''[cm/sec]  '')') v0
       write(*,'('' r0        = '',e11.5,''[cm]      '')') r0
       write(*,'('' t0        = '',e11.5,''[sec]     '')') t0
       write(*,'('' p0        = '',e11.5,''[erg/cm^3]'')') p0 / gamma0 ! keep compati. with old prog
       write(*,'('' b0        = '',e11.5,''[gauss]   '')') b0
       write(*,'('' rhoc      = '',e11.5,''[g/cm^3]  '')') rhoc
       write(*,'('' n0        = '',e11.5,''[/cm^3]'')') n0
       write(*,'('' r_top     = '',e11.5,''[r_sun]'')') rri(ii+1)/rs
       write(*,'('' r_bottom  = '',e11.5,''[r_sun]'')') rri(0   )/rs
       write(*,'('' r_surface = '',e11.5)') rs
       write(*,'('' Vcrtrn    = '',e11.5,''[cm/s]'')')    vcrtrn * v0
       write(*,'('' Mcrtrn    = '',e11.5,''[cm/s/cc]'')')
     &                                               mcrtrn * v0 * n0
       write(*,'(A)') '---'

       close(18)

       return
       endsubroutine

*
* -----------------------------------------------------------------------
* set coorinate system
* different from previous version of MHD comp.,
*     dr is uniform within source surface (2.5R) or further
*
* -----------------------------------------------------------------------
*
       subroutine setcrdnt(ii,jj,kk,
     &                     rbottom,rtop,
     &                     rri,thetai,phii,saini,kosai,sainbi,kosabi,
     &                     sinphi,cosphi,sinphbi,cosphbi)
       implicit none
* arguments
       intent(in)  ::      ii,jj,kk
       intent(in)  ::      rbottom,rtop
       intent(out) ::      rri,thetai,phii,saini,kosai,sainbi,kosabi
       intent(out) ::      sinphi,cosphi,sinphbi,cosphbi
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64) ! number of grid for each three direction.
       real    rbottom, rtop
       real    rri(-2:ii+3),thetai(0:jj),phii(-1:kk+2)
       real    saini(0:jj), kosai(0:jj)
       real    sainbi(0:jj-1), kosabi(0:jj-1)
       real    sinphi(0:kk+1), cosphi(0:kk+1)
       real    sinphbi(0:kk), cosphbi(0:kk)
* local
       real    pi
       parameter(pi = 3.1415926535897932385E+00)
*
       logical ldrfine
       parameter(ldrfine = .false.) ! if false, simply extend.
       logical ldrfix
       parameter(ldrfix = .false.) ! must be true when ldistant
       real    rrprov(0:ii+1)
       real    rrprov2(0:ii+1)
       real    thb, dr, dtheta, dphi
       integer i, j, k, j2
       integer iitmp
       real    rwb, aa, bb, cc, dd

* latitudinal
       dtheta = pi / float(jj)
       do j = 0, jj
         thetai(j) = dtheta * float(j) ! .... MIND some other parts will assume uniform grid
       enddo !                                   and seek float(jj) when non-uniform grid will be us

* sinu-funcs.
* Here, jj is assumed to be even, and
* latitudinal grids are assumed to be located in symmetry
*    with respect to the equatorial plane (thus, jj/2 is address for equator)
       j2 = jj / 2
       do j = 0, j2
         saini(j) = sin(thetai(j))
         saini(jj - j) =  saini(j)
         kosai(j) = cos(thetai(j))
         kosai(jj - j) = -kosai(j)
       enddo
       saini(0)   = 0.0E+00
       saini(jj)  = 0.0E+00
       saini(jj/2)= 1.0E+00  ! <== if jj is even
       kosai(0)   = 1.0E+00
       kosai(jj)  =-1.0E+00
       kosai(jj/2)= 0.0E+00

       do j = 0, jj / 2 - 1
         thb = (thetai(j+1) + thetai(j)) * 0.5E+00
         sainbi(j)      = sin(thb)
         sainbi(jj-1-j) = sainbi(j)
         kosabi(j)      = cos(thb)
         kosabi(jj-1-j) =-kosabi(j)
       enddo

* longitudinal
       dphi = 2.0E+00 * pi / float(kk)
       do k =-1, kk + 2
         phii(k) = float(k) * dphi ! MIND some other parts will assume uniform grid,
       enddo !                          and seek float(kk) when non-uniform grid will be used.

       if (kk .GT. 3) then
         do k = 0, kk
           aa = (phii(k+1) + phii(k)) * 0.5E+00
           sinphbi(k) = sin(aa)
           cosphbi(k) = cos(aa)
         enddo
         do k = 0, kk + 1
           aa = phii(k)
           sinphi(k) = sin(aa)
           cosphi(k) = cos(aa)
         enddo
       else
         aa = 0.5E+00 * pi ! dummy grid address.....
         bb = pi / float(jj) ! dummy d(phi)
         do k = 0, kk
           cc = aa + float(k * 2 - 1) * bb * 0.5E+00
           sinphbi(k) = sin(cc)
           cosphbi(k) = cos(cc)
         enddo
         do k = 0, kk + 1
*           cc = aa
           cc = aa + float(k-1) * bb
           sinphi(k) = sin(cc)
           cosphi(k) = cos(cc)
         enddo
       endif

* radial
       if (ldrfix) then ! do again....

         dr = (rtop - rbottom) / float(ii)
         rri(0) = rbottom
         do i = 1, ii
           rri(i) = rri(0) + dr * float(i)
         enddo

       else

         if ((ii .EQ. 72) .OR. (ii .EQ. 144)) then ! special....old rr(i)
           iitmp = 72
           bb = 5.0E-03 ! when r = 1.01 --- 50 Rs ! most frequently used.
           cc = 0.7E+00
           dd = 1.2E+00
           rwb= 1.2E+00
         else if ((ii .EQ. 128) .OR. (ii .EQ. 256)) then
           iitmp = 128
           bb = 2.0E-03 ! when r = 1.01 --- 50 Rs ! most frequently used.
           cc = 0.3E+00
           dd = 1.2E+00
           rwb= 1.3E+00
         else
           write(*,*) 'ii is not expected values, 72, 144, 128 or 256'
           stop
         endif
         rrprov(0) = rbottom
         do i = 1, ii + 1
           if (rrprov(i-1) .LT. rrprov(0) * dd) then
             dr = bb
           else
             aa = (rrprov(i-1) - rrprov(0)*dd)
     &                  / (rwb - rrprov(0)*dd)
             aa = 2.0E+00 / (exp(aa) + 1.0E+00)
             dr = bb * aa
     &          +(1.0E+00 - aa)
     &          * cc * sqrt(rrprov(i-1)*rrprov(0))
           endif
           if (dr .GT. cc) dr = cc
           rrprov(i) = rrprov(i - 1) + dr
         enddo
         aa = (rtop - rbottom) / (rrprov(iitmp) - rrprov(0))
         write(*,*) 'dr-modif factor ', aa

         rri(0) = rbottom

         if (ii .EQ. iitmp) then
           rri(0) = rrprov(0)
           do i = 1, ii + 1
             rri(i) = rri(0) + (rrprov(i) - rrprov(0)) * aa
           enddo
         else ! only OK when II = iitmp x 2

           if (ldrfine) then ! small dr

             do i = 0, ii,  2
               rrprov2(i) = rrprov(i/2)
             enddo
             do i = 1, ii+1, 2
               rrprov2(i) = (rrprov(i/2)+rrprov(i/2+1)) * 0.5E+00
             enddo
             rri(0) = rrprov2(0)
             do i = 1, ii + 1
               rri(i) = rri(0) + (rrprov2(i) - rrprov2(0)) * aa
             enddo

           else ! extent with sqrt(r), with limit

             rri(0) = rrprov(0)
             do i = 1, iitmp + 1
               rri(i) = rri(0) + (rrprov(i) - rrprov(0)) * aa
             enddo
*             iitemp = 49 ! adjust so that ii=144 can reach 500 Rs
             dr = rri(iitmp+1) - rri(iitmp)
             do i = iitmp+1, ii - 1
               aa = rri(i) / rri(iitmp)
               aa = sqrt(aa) * dr
               if (aa/rri(0) .GT. 8.0E+00) aa = rri(0) * 8.0E+00
               rri(i+1) = rri(i) + aa
             enddo

           endif ! end if (ldrfine)

         endif ! end if (ii .EQ....)

       endif ! end if ldrfix

       rri(ii+1) = 2.0E+00 * rri(ii  ) - rri(ii-1) ! outermost
       rri(ii+2) = 2.0E+00 * rri(ii+1) - rri(ii  ) ! outer ghost
       rri(ii+3) = 2.0E+00 * rri(ii+2) - rri(ii+1) ! outer ghost
       rri(  -1) = 2.0E+00 * rri(   0) - rri(   1) ! inner ghost
       rri(  -2) = 2.0E+00 * rri(  -1) - rri(   0) ! inner ghost

       return
       endsubroutine

*
* --------------------------------------------------------------------
*
* Get magnetic field using the Potential Field Approximation
*
* --------------------------------------------------------------------
*
* Fixed mesh method
*
* --------------------------------------------------------------------
*
       subroutine mshpotmg(ii,jj,kk,rt,
     &                     br2,bt2,bp2,
     &                     rr,theta,phi,sain,lapsi1,lapsi2)
       implicit none
* arguments
       intent(in)  ::      ii,jj,kk,rt
       intent(out) ::      br2,bt2,bp2
       intent(in)  ::      rr,theta,phi,sain
       intent(out) ::                        lapsi1,lapsi2
       integer rt
* large array
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64) ! number of grid for each three direction.
       real    br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    rr(-2:ii+3),theta(0:jj),phi(-1:kk+2)
       real    sain(0:jj)
       real    lapsi1(0:jj,0:kk+1,-1:ii+1)
       real    lapsi2(0:jj,0:kk+1,-1:ii+1)
* local
       real    brsurf(0:jj,0:kk+1)
       integer xwso, ywso
       parameter(xwso = 72, ywso = 30)
       real    brwso(0:xwso+1, 1:ywso)
       real    thewso(0:ywso+1), phiwso(0:xwso+1)
* switch
!       logical lcnfpf
!       parameter(lcnfpf = .false.) ! confined-pot.-field. SSPF if false.

!       if (lwsoraw) then
         call wsocord(thewso,phiwso) ! set WSO cooridnates
         call readwso(rt,brwso,kk)   ! read WSO data
         call reformbr(jj,kk,brwso,thewso,phiwso,brsurf, ! remap to match sim. grids
     &                 theta,phi,sain)
!       else
!         call gtmydmmy(jj,kk,brsurf,sain,rt) ! set dummy, for test
!       endif

!       if (lcnfpf) then
!       else
         call slvlapeq(ii,jj,kk,
     &                 brsurf, ! solve Poison's Eq. for SSPF
     &                 br2,bt2,bp2,
     &                 rr,theta,phi,sain,lapsi1,lapsi2)
!       endif

       return
       endsubroutine

*
* solve Laplace Eq. --------------------------------------------------
*
       subroutine slvlapeq(ii,jj,kk,
     &                     brsurf,
     &                     br2,bt2,bp2,
     &                     rr,theta,phi,sain,lapsi1,lapsi2)
       implicit none
* arguments
       intent(in)  ::      ii,jj,kk
       intent(in)  ::      brsurf
       intent(out) ::      br2,bt2,bp2
       intent(in)  ::      rr,theta,phi,sain
       intent(out) ::                        lapsi1,lapsi2
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       real    brsurf(0:jj,0:kk+1)
       real    br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    rr(-2:ii+3),theta(0:jj),phi(-1:kk+2)
       real    sain(0:jj)
       real    lapsi1(0:jj,0:kk+1,-1:ii+1)
       real    lapsi2(0:jj,0:kk+1,-1:ii+1)
* local
       real    rrb(-1:ii+1)
       integer itmax
       real    rw
       parameter(rw = 45.0E+00)
*       parameter(rw = 2.0E+00) ! in some of distant case...
*       parameter(rw = 2.5E+00)
       integer i, j, k, ibsrc, ittime, iescape
       real    mindt
       real    maxdiff
       real    difcrtrn
       real    aa, bb, cc, dd, ee, ff, gg
       real    rhs, lhs
       real    th2
       real    dth, dph

       dth = theta(1) - theta(0)
       dph = phi(1) - phi(0)

!$omp parallel do private(i,k,j)
       do i = -1, ii + 1
         do k = 0, kk + 1
         do j = 0, jj
           lapsi1(j,k,i) = 0.0E+00
           lapsi2(j,k,i) = 0.0E+00
         enddo
         enddo
       enddo
!$omp end parallel do
!$omp parallel do private(i)
       do i = -1, ii + 1
         if (i .EQ. ii + 1) then
*           rrb(i) = rr(i) ! usually not used...
           rrb(i) = rr(i) + (rr(i) - rr(i-1)) * 0.5E+00
         else if (i .EQ. -1) then
           rrb(i) = rr(i+1) - (rr(i+2) - rr(i+1)) * 0.5E+00
         else
           rrb(i) =(rr(i+1) + rr(i)) * 0.5E+00
         endif
       enddo
!$omp end parallel do

       mindt = (rrb(0) * sain(1) * dph)**2
       itmax = jj * 4000 ! 1.0E+01 / mindt
*       itmax = 1.0E+01 / mindt
       write(*,*) ' Min Mesh_dt / Max Iteration count = ', mindt,itmax

* define mesh point(s) of source surface ! assuing rr(0) < rw < rr(ii+1)
       do i = 0, ii
           if ((rw .GE. (rr(i  )/rr(0))) .AND.
     &         (rw .LE. (rr(i+1)/rr(0)))) ibsrc = i
       enddo
       write(*,*) 'Ib SouRCe = ', ibsrc

* initialize psi
!$omp parallel do private(i,k,j)
         do i = -1, ii + 1
           if (i .LE. ibsrc) then
             do k =  0, kk + 1
             do j =  0, jj
               lapsi1(j,k,i) = - brsurf(j,k) * (rrb(i) - rrb(ibsrc))
             enddo
             enddo
           endif
         enddo
!$omp end parallel do

* iterration process
       ittime = 0
       iescape = -1
       do while ((ittime .LE. itmax) .AND. (iescape .EQ. -1)) ! ----
         ittime = ittime + 1
!$omp parallel do private(i,k,j,aa,bb,cc,dd,ee,ff,gg,lhs,rhs,th2)
         do i = 0, ii + 1
           if (i .LE. ibsrc-1) then
             do k =  1, kk
             do j =  1, jj - 1
               th2 = (theta(j+1) + theta(j)) * 0.5E+00
               ee  = sin(th2)
               th2 = (theta(j) + theta(j-1)) * 0.5E+00
               ff  = sin(th2)
               aa  = 2.0E+00 / rrb(i)**2 / (rrb(i+1) - rrb(i-1))
               bb  = (rrb(i+1) + rrb(i))**2 * 0.25E+00
     &             / (rrb(i+1) - rrb(i))
               cc  = (rrb(i) + rrb(i-1))**2 * 0.25E+00
     &             / (rrb(i) - rrb(i-1))
               dd  = 1.0E+00 / (rrb(i)**2 * sain(j) * dth**2)
               gg  = 1.0E+00 / (rrb(i)    * sain(j) * dph   )**2
               rhs = aa * (bb * lapsi1(j,k,i+1) + cc * lapsi1(j,k,i-1))
     &             + dd * (ee * lapsi1(j+1,k,i) + ff * lapsi1(j-1,k,i))
     &             + gg * (     lapsi1(j,k+1,i) +      lapsi1(j,k-1,i))
               lhs = aa * (bb + cc) + dd * (ee + ff) + gg * 2.0E+00
               lapsi2(j,k,i) = rhs / lhs
             enddo
             enddo
           endif
         enddo
!$omp end parallel do

!$omp parallel do private(i,k,j)
           do i = -1, ii + 1
             if (i .EQ. -1) then
               do k =  1, kk
               do j =  1, jj - 1
                 lapsi2(j,k,i) = brsurf(j,k) * (rrb(i+1) - rrb(i))
     &                         + lapsi2(j,k,i+1)
               enddo
               enddo
             endif
             if (i .EQ. ibsrc) then
               do k =  1, kk
               do j =  1, jj - 1
                 lapsi2(j,k,i) = 0.0E+00
               enddo
               enddo
             endif
           enddo
!$omp end parallel do

!$omp parallel do private(i,k,j,aa,bb)
         do i = -1, ii + 1
           if (i .LE. ibsrc) then
* rotational axis
             aa = 0.0E+00
             bb = 0.0E+00
             do k = 1, kk
               aa = aa + lapsi2(   1,k,i)
               bb = bb + lapsi2(jj-1,k,i)
             enddo
             aa = aa / float(kk)
             bb = bb / float(kk)
             do k = 1, kk
               lapsi2( 0,k,i) = aa
               lapsi2(jj,k,i) = bb
             enddo
* overlap along longitude
             do j =  0, jj
               lapsi2(j,   0,i) = lapsi2(j,kk,i)
               lapsi2(j,kk+1,i) = lapsi2(j, 1,i)
             enddo
           endif
         enddo
!$omp end parallel do

* copy Psi and estimate max. difference
         maxdiff = -1.0E+10
!$omp parallel do private(i,k,j,aa)
!$omp&            reduction(max:maxdiff)
         do i = -1, ii + 1
           if (i .LE. ibsrc) then
             do k =  0, kk + 1
             do j =  0, jj
               aa = 1.0E+00
     &            - (abs(lapsi2(j,k,i)) + 1.0E-05)
     &            / (abs(lapsi1(j,k,i)) + 1.0E-05)
               aa = abs(aa)
               if (aa .GT. maxdiff) maxdiff = aa
               lapsi1(j,k,i) = lapsi2(j,k,i)
             enddo
             enddo
           endif
         enddo ! end i-loop
!$omp end parallel do
*
         if (mod(ittime,1000) .EQ. 0)
     &     write(*,'(i6,'' : '',e11.5)') ittime, maxdiff

         if (kk .LT. 2) then
           difcrtrn = 1.0E-05
         else
           difcrtrn = 1.0E-04
*           difcrtrn = 1.0E-03 ! MIND the criterion
         endif
         if (maxdiff .LT. difcrtrn) then
           write(*,*) ' converged.....'
           iescape = 1
           write(*,'(i6,'' : '',e11.5)') ittime, maxdiff
         endif

       enddo ! END of DO-WHILE-LOOP                                ----

       write(*,*) 'iteration .. done '
       write(*,'('' Ittime = '',i6)') ittime

* Obtain potential field
!$omp parallel do private(i,k,j)
       do i = 0, ii + 1
         if (i .LE. ibsrc) then
           do k = 1, kk
           do j = 1, jj - 1
             br2(j,k,i) =- (lapsi2(j  ,k,i  )-lapsi2(j  ,k,i-1))
     &                  / (rrb(i) - rrb(i-1))
             bt2(j,k,i) =-((lapsi2(j+1,k,i  )-lapsi2(j-1,k,i  ))
     &                     /rrb(i  )
     &                    +(lapsi2(j+1,k,i-1)-lapsi2(j-1,k,i-1))
     &                     /rrb(i-1))
     &                  /(4.0E+00 * dth)
             bp2(j,k,i) =-((lapsi2(j,k+1,i  )-lapsi2(j,k-1,i  ))
     &                     /rrb(i  )
     &                    +(lapsi2(j,k+1,i-1)-lapsi2(j,k-1,i-1))
     &                     /rrb(i-1))
     &                  /(4.0E+00 * dph * sain(j))
           enddo
           enddo
         endif
       enddo
!$omp end parallel do

* rotational axis
!$omp parallel do private(i,k,aa,bb)
       do i = 0, ii + 1
         if (i .LE. ibsrc) then
           aa = - (lapsi2( 0,1,i) - lapsi2( 0,1,i-1))
     &          / (rrb(i) - rrb(i-1))
           bb = - (lapsi2(jj,1,i) - lapsi2(jj,1,i-1))
     &          / (rrb(i) - rrb(i-1))
           do k = 1, kk
             br2( 0,k,i) = aa
             br2(jj,k,i) = bb
             bt2( 0,k,i) = 0.0E+00
             bt2(jj,k,i) = 0.0E+00
             bp2( 0,k,i) = 0.0E+00
             bp2(jj,k,i) = 0.0E+00
           enddo
         endif
       enddo
!$omp end parallel do

* beyond the source surface

!$omp parallel do private(i,k,j)
       do i = 0, ii + 1
         if (i .GT. ibsrc) then
           do k = 1, kk
           do j = 0, jj
             br2(j,k,i) = br2(j,k,ibsrc) * (rr(ibsrc)/rr(i))**2
             bt2(j,k,i) = 0.0E+00
             bp2(j,k,i) = 0.0E+00
           enddo
           enddo
         endif
       enddo
!$omp end parallel do

* overlap along longitudinal direction
!$omp parallel do private(i,j)
       do i = 0, ii + 1
         do j = 0, jj
           br2(j,   0,i) = br2(j,kk,i)
           br2(j,kk+1,i) = br2(j, 1,i)
           bt2(j,   0,i) = bt2(j,kk,i)
           bt2(j,kk+1,i) = bt2(j, 1,i)
           bp2(j,   0,i) = bp2(j,kk,i)
           bp2(j,kk+1,i) = bp2(j, 1,i)
         enddo
       enddo
!$omp end parallel do

       return
       endsubroutine

*
** eformat surface field data --------------------------------------
*
       subroutine reformbr(jj,kk,brwso,thewso,phiwso,brsurf,
     &                     theta,phi,sain)
       implicit none
* arguments
       intent(in)  ::      jj,kk
       intent(in)  ::            brwso,thewso,phiwso
       intent(out) ::                                brsurf
       intent(in)  ::      theta,phi,sain
       integer xwso, ywso
       parameter(xwso = 72, ywso = 30)
       real    brwso(0:xwso+1, 1:ywso)
       real    thewso(0:ywso+1), phiwso(0:xwso+1)
       integer jj, kk
*       parameter(jj = 32, kk = 64)
       real    brsurf(0:jj,0:kk+1)
       real    theta(0:jj),phi(-1:kk+2)
       real    sain(0:jj)
* local
       real    brsurfb(0:jj,-1:kk+2), wei(0:jj,-1:kk+2)
       integer j, k, iwso, jwso
       real    sumbr1, sumbr2
       integer j1, j2, k1, k2
       real    aa, bb, cc, dd, ee, ff
       real    x0, y0, z0, x2, y2, z2
       real    th, ph
       real    dth1, dth2, dph1, dph2
       real    domega
* functions defined by myself
       interface
         pure function domeg(jj,kk,j,theta)
           implicit none
           real        domeg
           intent(in) ::     jj,kk,j,theta
           integer j
           integer jj, kk
*           parameter(jj = 32, kk = 64)
           real    theta(0:jj)
         endfunction
       endinterface
* const. and switch
       real    pi
       parameter(pi = 3.1415926535897932385E+00)
!       logical llinear ! true : linear->spread (does not work properly), false : spread directly
!       parameter(llinear = .false.)
       logical loffset
       parameter(loffset = .true.)

* linear extrapolation : Stanford Synoptic Chart => Mesh System
!       if (llinear) then
!       else

         do k = 1, kk
           do j = 1, jj - 1
             x0 = sin(theta(j)) * cos(phi(k))
             y0 = sin(theta(j)) * sin(phi(k))
             z0 = cos(theta(j))
             brsurf(j,k) = 0.0E+00
             wei(j,k)    = 0.0E+00
             do k2 = 1, xwso
             do j2 = 1, ywso
               x2 = sin(thewso(j2)) * cos(phiwso(k2))
               y2 = sin(thewso(j2)) * sin(phiwso(k2))
               z2 = cos(thewso(j2))
               aa = x0 * x2 + y0 * y2 + z0 * z2
               if (aa .GT. 0.8E+00) then
                 aa = aa * 0.999E+00
                 aa = acos(aa) / pi * 180.0E+00 ! angle in degree.
                 bb = cos(thewso(j2))
                 bb = bb**4
                 cc = brwso(k2,j2) / 500.0E+00 ! in 5 gauss
                 cc = cc**2
                 if (cc .GT. 1.0E+00) cc = 1.0E+00
                 dd =  5.0E+00 !      basic half width
     &              + 10.0E+00 * bb ! make large half width near pole
     &              + 10.0E+00 * cc ! make large half width for strong field spot
                 ee = acos(0.8E+00)
                 ee = ee**2
                 ee =(1.0E+00 - exp(-ee)) * dd**2 * 0.5E+00
*                 ff = exp(-aa**2) ! approx. linear
                 aa = aa / dd
                 ff = exp(-aa**2) / ee
               else
                 ff = 0.0E+00
               endif
               brsurf(j,k) = brsurf(j,k)
     &                     + ff * brwso(k2,j2) ! / sin(thewso(j2)) ! radial&area
               wei(j,k)    = wei(j,k) + ff * sin(thewso(j2))
             enddo
             enddo
             brsurf(j,k) = brsurf(j,k) / wei(j,k)
           enddo
           brsurf( 0,k) = brsurf(   1,k)
           brsurf(jj,k) = brsurf(jj-1,k)
         enddo
         do j = 0, jj
           brsurf(j,   0) = brsurf(j,kk)
           brsurf(j,kk+1) = brsurf(j, 1)
         enddo

!       endif

* extract momopole
       sumbr1 = 0.0E+00
       do j = 0, jj
         domega = domeg(jj,kk,j,theta)
         do k = 1, kk
           sumbr1 = sumbr1 + brsurf(j,k) * domega
         enddo
       enddo
*
       if (loffset) then
         sumbr2 = 0.0E+00
         do j = 0, jj
           domega = domeg(jj,kk,j,theta)
           do k = 0, kk + 1
             brsurf(j,k) = brsurf(j,k) - sumbr1 * 0.5E+00
             if ((k .GE. 1) .AND. (k .LE. kk))
     &         sumbr2 = sumbr2 + brsurf(j,k) * domega
           enddo
         enddo
       else
         sumbr2 = sumbr1
       endif
       write(*,'('' SumBr : '',e11.5,'' => '',e11.5)') sumbr1,sumbr2

       return
       endsubroutine

*
** weight at each mesh point for integration of B_r -----------------
*
       pure function domeg(jj,kk,j,theta)
       implicit none
       real          domeg
* arguments
       intent(in) ::       jj,kk,j,theta
       integer j
       integer jj, kk
*       parameter(jj = 32, kk = 64)
       real    theta(0:jj)
* local
       real    th1, th2, aa

       if (j .EQ. 0) then
         th1 = theta( 0)
         th2 =(theta( 1) + theta(   0)) * 0.5E+00
       endif
       if (j .EQ. jj) then
         th1 =(theta(jj) + theta(jj-1)) * 0.5E+00
         th2 = theta(jj)
       endif
       if((j .GT. 0) .AND. (j .LT. jj)) then
         th1 =(theta(j) + theta(j-1)) * 0.5E+00
         th2 =(theta(j) + theta(j+1)) * 0.5E+00
       endif
       aa = (cos(th1) - cos(th2)) / float(kk)
       domeg = aa
       return
       endfunction

*
** find location on Stan. Synoptic Map. -----------------------------
*
       subroutine tp2jkwso(th,ph,iwso,jwso,thewso,phiwso)
       implicit none
* arguments
       intent(in)  ::      th,ph
       intent(in)  ::                      thewso,phiwso
       intent(out) ::            iwso,jwso
       real    th, ph
       integer iwso, jwso
       integer xwso, ywso
       parameter(xwso = 72, ywso = 30)
       real    thewso(0:ywso+1), phiwso(0:xwso+1)
* local
       real    pi
       parameter(pi = 3.1415926535897932385E+00)
       integer m
       real    ph2, th2

       th2 = th
       ph2 = ph

       if (ph2 .GT. 2.0E+00 * pi) ph2 = ph2 - 2.0E+00 * pi
       if (ph2 .LT. 0.0E+00     ) ph2 = ph2 + 2.0E+00 * pi
*       ph2 = ph2 + 2.0E+00 * pi
*       ph2 = mod(ph2, 2.0E+00 * pi)

       do m = 0, xwso
         if ((ph2 .GE .phiwso(m+1)) .AND.
     &       (ph2 .LE .phiwso(m  ))) iwso = m
       enddo

       if  (th2 .LT. thewso(     1)) jwso = 1
       if  (th2 .GT. thewso(ywso-1)) jwso = ywso - 1
       if ((th2 .GE. thewso(     1)) .AND.
     &     (th2 .LE. thewso(ywso-1))) then
         do m = 1, ywso - 1
           if ((th2 .GE .thewso(m  )) .AND.
     &         (th2 .LE .thewso(m+1))) jwso = m
         enddo
       endif
       return
       endsubroutine

*
** read synoptic chart of solar magnetic field ----------------------
*
       subroutine readwso(nrot,brwso,kk)
       implicit none
* arguments
       intent(out) ::          brwso
       intent(in)  ::     nrot,      kk
       integer nrot
       integer xwso, ywso
       parameter(xwso = 72, ywso = 30)
       real    brwso(0:xwso+1, 1:ywso)
* local
       integer kk
*       parameter(kk = 64)
       real    aa
       integer i, j
       character*48 flname

       write(flname,'(''stn'',i4,''.dat'')') nrot
       open(unit = 1, file = flname ,status = 'old')
       write(*,*) 'Now Open File : ', flname
       do i = 1, xwso
         read(1,10) (brwso(i,j), j =  1,  6)
         read(1,20) (brwso(i,j), j =  7, 14)
         read(1,20) (brwso(i,j), j = 15, 22)
         read(1,20) (brwso(i,j), j = 23, 30)
       enddo
   10  format(18x,6f9.3)
   20  format(8f9.3)
       close(1)

* if 2D case, average along longitude
       if (kk .LT. 2) then
         do  j = 1, ywso
           aa = 0.0E+00
           do i = 1, xwso
             aa = aa + brwso(i,j)
           enddo
           aa = aa / 30.0E+00
           do i = 1, xwso
             brwso(i,j) = aa
           enddo
         enddo
       endif

* overlap
       do j = 1, ywso
         brwso(     0,j) = brwso(xwso,j)
         brwso(xwso+1,j) = brwso(   1,j)
       enddo

       return
       endsubroutine

*
** set stanford WSO coordinates  --------------------------------------

       subroutine wsocord(thewso,phiwso)
       implicit none
* arguments
       intent(out) ::     thewso,phiwso
       integer xwso, ywso
       parameter(xwso = 72, ywso = 30)
       real    thewso(0:ywso+1), phiwso(0:xwso+1)
* local
       real    aa
       integer j, k
       real    pi
       parameter(pi = 3.1415926535897932385E+00)

* wso theta : latitude
       do j = 1, ywso
         aa = float(31 - 2 * j) / float(ywso)
         thewso(j) = acos(aa)
       enddo
       thewso(0)       = - thewso(1)
       thewso(ywso+1) = 2.0E+00 * pi - thewso(ywso)
* wso phi   : longitude
       do k = 0, xwso + 1
         phiwso(k) = 2.0E+00 * pi
     &             *(1.0E+00 - float(k - 1) / float(xwso))
       enddo

       write(*,*) 'Stan. Coordinate : completed'

       return
       endsubroutine

*
* --------------------------------------------------------------------
*
*   Legendre function
*
* --------------------------------------------------------------------
*
       subroutine setlgndr(jj,imagdata,sain,kosa,lgndr,dlgndr)
       implicit none
* arguments
       intent(in)  ::      jj,imagdata,sain,kosa
       intent(out) ::                            lgndr,dlgndr
       integer imagdata
       integer jj
*       parameter(jj = 32)
       real    sain(0:jj), kosa(0:jj)
       real    lgndr(1:10,0:10,0:jj), dlgndr(1:10,0:10,0:jj) ! 10 is fixed.....
* local
       integer n, m, j
       real    xx
* functions defined by myself
       interface
         pure function mhdpnmx(l,m,x)
           implicit none
           real        mhdpnmx
           intent(in)  ::   l,m,x
           integer l, m
           real    x
         endfunction
         pure function mhddpnmx(l0,m0,x)
           implicit none
           real        mhddpnmx
           intent(in)  ::    l0,m0,x
           integer l0, m0
           real    x
         endfunction
       endinterface

       xx = sain(1) ! supress warnings

       if (imagdata .LE. 0) then
         do j = 0, jj
         do m = 0, 10
         do n = 1, 10
           lgndr(n,m,j)  =  1.0E+20 ! should not be used
           dlgndr(n,m,j) =  1.0E+20 ! should not be used
         enddo
         enddo
         enddo
       else
         do j = 0, jj
         do m = 0, 10
         do n = 1, 10
           if (m .LE. n) then
             xx = kosa(j)
             lgndr(n,m,j) = mhdpnmx(n,m,xx)
             if ((j .EQ. 0) .OR. (j .EQ. jj)) then
               dlgndr(n,m,j) = 0.0E+00
             else
               dlgndr(n,m,j) = mhddpnmx(n,m,xx)
             endif
           else
             lgndr(n,m,j) = 0.0E+00
             dlgndr(n,m,j) = 0.0E+00
           endif
         enddo
         enddo
         enddo
       endif
       return
       endsubroutine

*
* --------------------------------------------------------------------
*
*   Legendre expansion coefficients -> magnetic field
*
* --------------------------------------------------------------------
*
       subroutine expndmag(ii,jj,kk,
     &                     pota, potb,
     &                     nend,rt,rs,
     &                     br2,bt2,bp2,
     &                     rr,phi,sain,kosa,lgndr,dlgndr)
       implicit none
* arguments
       intent(in)  ::      ii,jj,kk
       intent(in)  ::      nend,rt,rs
       intent(out) ::      br2,bt2,bp2
       intent(in)  ::      rr,phi,sain,kosa,lgndr,dlgndr
       integer nend, rt
       real    rs
* local
       real    rw
       parameter(rw = 45.0E+00)      ! Radius of source surface
*       parameter(rw = 2.0E+00) ! in some of distant case...
*       parameter(rw = 2.5E+00)
*
       real    aa, bb, cc
       real    ra, ph, sin0, cos0
       integer n, m, i, j, k, j2
       real    pota(0:10,0:10), potb(0:10,0:10)
* clone
       real    potal(0:10,0:10), potbl(0:10,0:10)
       integer nendl
* common
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64) ! number of grid for each three direction.
       real    br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    rr(-2:ii+3),phi(-1:kk+2)
       real    sain(0:jj), kosa(0:jj)
       real    lgndr(1:10,0:10,0:jj), dlgndr(1:10,0:10,0:jj) ! 10 is fixed.....

** initial potential magnetic field coefficient

* suppress warnings
       m = rt
       call potcoef(pota, potb, rt) ! load coef.

* convert potential coefficients to field
!$omp parallel do private(i,k,j,n,m,j2)
!$omp&            private(potal,potbl,nendl)
!$omp&            private(aa,bb,cc,ra,ph)
!$omp&            private(sin0,cos0)
       do i = -2, ii + 3 !  -1, ii + 1

         do m = 0, 10
         do n = 0, 10
           potal(n,m) = pota(n,m)
           potbl(n,m) = potb(n,m)
         enddo

         enddo
         do k = 1, kk
         do j = 0, jj
           if (i .EQ. -1) then
             ra = (rr(0) * 2.0E+00 - rr(1)) / rs
           else
             ra = rr(i) / rs
           endif
           nendl = nend
           j2 = j
           sin0 = sain(j)
           cos0 = kosa(j)
           ph = phi(k)
           if (ra .LE. rw) then
             call potco2mg(jj,kk,
     &                     potal,potbl,aa,bb,cc,ra,sin0,ph,nendl,j2,
     &                     lgndr,dlgndr)
           else
             ra = rw
             call potco2mg(jj,kk,
     &                     potal,potbl,aa,bb,cc,ra,sin0,ph,nendl,j2,
     &                     lgndr,dlgndr)
             aa = aa / (rr(i) / rs)**2 * rw**2
             bb = 0.0E+00
             cc = 0.0E+00
           endif
           br2(j,k,i) = aa
           bt2(j,k,i) = bb
           bp2(j,k,i) = cc
         enddo
         enddo
* longi
         do j = 0, jj
           br2(j,   0,i) = br2(j,kk,i)
           br2(j,kk+1,i) = br2(j, 1,i)
           bt2(j,   0,i) = bt2(j,kk,i)
           bt2(j,kk+1,i) = bt2(j, 1,i)
           bp2(j,   0,i) = bp2(j,kk,i)
           bp2(j,kk+1,i) = bp2(j, 1,i)
         enddo

       enddo ! end i-loop
!$omp end parallel do

       return
       endsubroutine

*
* --------------------------------------------------------------------
*
* Potential field coefficients ==> magnetic field [micro T]
*
* --------------------------------------------------------------------
*
       subroutine potco2mg(jj,kk,
     &                     pota,potb,magr,magt,magp,ra,
     &                     sin0,ph,nend,j,lgndr,dlgndr)
       implicit none
* arguments
       intent(in)  ::      jj,kk
       intent(in)  ::      pota,potb,               ra
       intent(in)  ::      sin0,ph,nend,j,lgndr,dlgndr
       intent(out) ::                magr,magt,magp
       real    pota(0:10,0:10), potb(0:10,0:10)
       real    magr, magt, magp    ! magnetic field
       real    ra, ph, sin0 ! , cos0
       integer j, nend   ! address of latitude, and degree of expansion
       integer jj, kk
*       parameter(jj = 32, kk = 64)
       real    lgndr(1:10,0:10,0:jj), dlgndr(1:10,0:10,0:jj) ! 10 is fixed.....
* local
       real    mphi
       integer n, m, me
       real    rrn, drrn, ppm, dppm, ttnm, dttnm
       real    drpot, dtpot, dppot
       real    rw, cn
       parameter(rw = 45.0E+00)
*       parameter(rw = 2.0E+00) ! in some of distant case...
*       parameter(rw = 2.5E+00)
       real    a ! , b

** estimate mag : rad theta phi
       drpot = 0.0E+00
       dtpot = 0.0E+00
       dppot = 0.0E+00
       do n = 1, nend
         cn = -1.0E+00 / (rw**(2 * n + 1) - 1.0E+00)
         rrn  = cn * ra**n
     &        +(1.0E+00 - cn) / ra**(n + 1)
         drrn = cn * ra**(n - 1) * float(n)
     &        -(1.0E+00 - cn) / ra**(n + 2) * float(n + 1)
         if (kk .GT. 3) then
           me = n
         else
           me = 0 ! pick up only primary terms (m = 0)
         endif
         do m = 0, me
           mphi = float(m) * ph
           ppm   = pota(n, m) * cos(mphi) + potb(n, m) * sin(mphi)
           dppm  = float(m)
     &           *(potb(n, m) * cos(mphi) - pota(n, m) * sin(mphi))
           ttnm  = lgndr(n, m, j)
           a = sin0
           a = abs(a)
           if (a .LT. 1.0E-04) then
             dttnm = 0.0E+00
           else
             dttnm = - sin0 * dlgndr(n, m, j)
           endif
           drpot = drpot + drrn *  ttnm *  ppm
           dtpot = dtpot +  rrn * dttnm *  ppm
           dppot = dppot +  rrn *  ttnm * dppm
         enddo
       enddo
       magr = - drpot
       magt = - dtpot / ra
       a = sin0
       a = abs(a)
       if (a .LT. 1.0E-02) magp = 0.0E+00
       if (a .GE. 1.0E-02) magp = - dppot / ra / sin0

       if (kk .LT. 3) then ! adjustment
         magr = magr * 2.0E+00
         magt = magt * 2.0E+00
         magp = magp * 2.0E+00
       endif

       return
       endsubroutine

* --------------------------------------------------------------------
**  The normalized Legendre Polynomials
* --------------------------------------------------------------------
*
       pure function mhdpnmx(l,m,x)
       implicit none
       real          mhdpnmx
* arguments
       intent(in)  ::     l,m,x
       integer l, m
       real    x
* local
       real    pmm, pll, pmmp1, somx2, ans, fact
       real    aan, bbn
       integer i, ll

       pmm = 1.0E+00
       if (m .GT. 0) then
         somx2 = (1.0E+00 - x) * (1.0E+00 + x)
         if (somx2 .GT. 0.0E+00) then
           somx2 = sqrt(somx2)
         else
           somx2 = 0.0E+00
         endif
         fact = 1.0E+00
         do i = 1, m
           pmm = pmm * fact * somx2
           fact = fact + 2.0E+00
         enddo
       endif
       if (l .EQ. m) then
         ans = pmm
       else
         pmmp1 = x * float(2 * m + 1) * pmm
         if (l .EQ. (m + 1)) then
           ans = pmmp1
         else
           do ll = (m + 2), l
             pll = (x * float(ll * 2 - 1) * pmmp1 -
     &               float(ll + m - 1) * pmm) / float(ll - m)
             pmm = pmmp1
             pmmp1 = pll
           enddo
           ans = pll
         endif
       endif
* normalize
       if (m .GT. 0) then
         aan = 1.0E+00
         bbn = 1.0E+00
         do i = (l - m + 1), (l + m)
           aan = aan * float(i)
         enddo
         ans = ans * sqrt(2.0E+00 * bbn / aan)
       endif

       mhdpnmx = ans

       return
       endfunction

*
* --------------------------------------------------------------------
* normalized dPnm(n, m, x) / dx
* --------------------------------------------------------------------
*
       pure function mhddpnmx(l0,m0,x)
       implicit none
       real          mhddpnmx
* arguments
       intent(in)  ::      l0,m0,x
       integer l0, m0
       real    x
* local
       real    pmm, pll, pmmp1, somx2, fact
       real    aan, bbn
       integer i, ll, l, m
       real    ans1, ans2, dans

* get non-normalized mhdpnmx(n, m, x)
       l = l0
       m = m0
       pmm = 1.0E+00
       if (m .GT. 0) then
         somx2 = (1.0E+00 - x) * (1.0E+00 + x)
         if (somx2 .GT. 0.0E+00) then
           somx2 = sqrt(somx2)
         else
           somx2 = 0.0E+00
         endif
         fact = 1.0E+00
         do i = 1, m
           pmm = pmm * fact * somx2
           fact = fact + 2.0E+00
         enddo
       endif
       if (l .EQ. m) then
         ans1 = pmm
       else
         pmmp1 = x * float(2 * m + 1) * pmm
         if (l .EQ. (m + 1)) then
           ans1 = pmmp1
         else
           do ll = (m + 2), l
             pll = (x * float(ll * 2 - 1) * pmmp1 -
     &               float(ll + m - 1) * pmm) / float(ll - m)
             pmm = pmmp1
             pmmp1 = pll
           enddo
           ans1 = pll
         endif
       endif
* get non-normalized mhdpnmx(n + 1, m, x)
       l = l0 + 1
       m = m0
       pmm = 1.0E+00
       if (m .GT. 0) then
         somx2 = (1.0E+00 - x) * (1.0E+00 + x)
         if (somx2 .GT. 0.0E+00) then
           somx2 = sqrt(somx2)
         else
           somx2 = 0.0E+00
         endif
         fact = 1.0E+00
         do i = 1, m
           pmm = pmm * fact * somx2
           fact = fact + 2.0E+00
         enddo
       endif
       if (l .EQ. m) then
         ans2 = pmm
       else
         pmmp1 = x * float(2 * m + 1) * pmm
         if (l .EQ. (m + 1)) then
           ans2 = pmmp1
         else
           do ll = (m + 2), l
             pll = (x * float(ll * 2 - 1) * pmmp1 -
     &               float(ll + m - 1) * pmm) / float(ll - m)
             pmm = pmmp1
             pmmp1 = pll
           enddo
           ans2 = pll
         endif
       endif
* get non-normalized mhddpnmx(n, m, x)/dx
       m = m0
       l = l0
       dans = float(l + 1) * x * ans1 - float(l - m + 1) * ans2
       dans = dans / ((1.0E+00 - x) * (1.0E+00 + x))
* normalize
       if (m .GT. 0) then
         aan = 1.0E+00
         bbn = 1.0E+00
         do i = (l - m + 1), (l + m)
           aan = aan * float(i)
         enddo
         dans = dans * sqrt(2.0E+00 * bbn / aan)
       endif

       mhddpnmx = dans

       return
       endfunction

*
** -----------------------------------------------------------
*
       subroutine getpark(tmpfact,r1,vpark,dpark,gamma0)
       implicit none
* arguments
       intent(in)  ::     tmpfact,r1,            gamma0
       intent(out) ::                vpark,dpark
       real    tmpfact, r1, vpark, dpark, gamma0
* local
       real    r2, vp, vg, dd
* functions defined by myself
       interface
         pure function ff1prk(rad,vr,gamma0)
           implicit none
           real        ff1prk
           intent(in)  ::     rad,vr,gamma0
           real   rad, vr, gamma0
         endfunction
         pure function dff1prk(rad,vr,gamma0)
           implicit none
           real        dff1prk
           intent(in)  ::      rad,vr,gamma0
           real   rad, vr, gamma0
         endfunction
       endinterface

       if ((gamma0 .GE. 1.0E+00) .AND. (gamma0 .LT. 1.6667E+00)) then
         r2 = r1 * tmpfact
         vp = 1.0E-04
         if (r2 .GT. 1.0E+00) vp = 1.0E+02
         dd = 1.0E+02
         do while(dd .GE. 1.0E-04)
           vg = vp - ff1prk(r2,vp,gamma0) / dff1prk(r2,vp,gamma0)
           dd = abs(1.0E+00 - vg / vp)
           vp = vg
         enddo
         dpark = 1.0E+00 / vp / r2**2
         vpark = vp * sqrt(tmpfact)
       else
         dpark = 1.0E+00
         vpark = 1.0E+00
       endif

       return
       endsubroutine

*
** -----------------------------------------------------------
*
       pure function ff1prk(rad,vr,gamma0)
       implicit none
       real          ff1prk
* arguments
       intent(in)  ::       rad,vr,gamma0
       real   rad, vr, gamma0
* local
       real   ans, aa
       real   cc
       parameter(cc = 0.0E+00) ! offset of enthalpy..

       aa = gamma0 - 1.0E+00
       if (abs(aa) .LT. 1.0E-05) then
         ans = 0.5E+00 *(vr**2 - 1.0E+00)
     &       - log(vr) - 2.0E+00 * log(rad)
     &       - 2.0E+00 *(1.0E+00 / rad - 1.0E+00) + cc
       else
         ans = 0.5E+00 *(vr**2 - 1.0E+00)
     &       +((rad**2 * vr)**(-aa) - 1.0E+00) / aa
     &       - 2.0E+00 *(1.0E+00 / rad - 1.0E+00) + cc
       endif
       ff1prk = ans

       return
       endfunction

*
** -----------------------------------------------------------
*
       pure function dff1prk(rad,vr,gamma0)
       implicit none
       real          dff1prk
* arguments
       intent(in)  ::        rad,vr,gamma0
       real   rad, vr, gamma0
* local
       real   ans, aa

       aa = gamma0 - 1.0E+00
       if (abs(aa) .LT. 1.0E-05) then
         ans = vr - 1.0E+00 / vr
       else
         ans = vr - rad**(-2.0E+00*aa) * vr**(-gamma0)
       endif
       dff1prk = ans

       return
       endfunction

*
* ---------------------------------------------------------------------
* READ coefficients for potential field expansion
* ---------------------------------------------------------------------
*
       subroutine potcoef(pota,potb,rt)
       implicit none
* arguments
       intent(in)  ::               rt
       intent(out) ::     pota,potb
       real    pota(0:10,0:10), potb(0:10,0:10)
       integer rt
* local
       integer i, j, nn
       logical ldummy1, ldummy2
       character*64 cdummy
       character*48 flname1, flname2
       do i = 0, 10
         do j = 0, 10
           pota(j,i) = 0.0E+00
           potb(j,i) = 0.0E+00
         enddo
       enddo
* try to open POT file
       write(flname1,'(''pot'',i4,''c.dat'')') rt
       inquire(file = flname1, exist = ldummy1)
       write(flname2,'(''pot'',i4,''d.dat'')') rt
       inquire(file = flname2, exist = ldummy2)
       if (ldummy1) then
         open(unit = 1, file = flname1, status = 'old')
         write(*,*) 'Now opened file = ', flname1
       else if (ldummy2) then
         open(unit = 1, file = flname2, status = 'old')
         write(*,*) 'Now opened file = ', flname2
         flname1 = flname2
       else
         write(*,*) 'file not found'
       endif

       read(1,'(A)') cdummy
       read(1,'(A)') cdummy
       do i = 0, 10
         read(1,33) nn, (pota(j, i), j = 1, 5)
       enddo
       read(1,'(A)') cdummy
       do i = 0, 10
         read(1,33) nn, (pota(j, i), j = 6, 10)
       enddo
       read(1,'(A)') cdummy
       read(1,'(A)') cdummy
       do i = 0, 10
         read(1,33) nn, (potb(j, i), j = 1, 5)
       enddo
       read(1,'(A)') cdummy
       do i = 0, 10
         read(1,33) nn, (potb(j, i), j = 6, 10)
       enddo
 33    format(i3,5e14.7)
       close(1)

       write(*,*) 'Now closed file = ', flname1
       return
       endsubroutine

*
* --------------------------------------------------------------------
*
* --------------------------------------------------------------------
*
* Copy properties
* Estimate velocities of magnetosonic mode
*   and determnine time step 'dt'
*
* --------------------------------------------------------------------
*
       subroutine copyvmax(ii,jj,kk,
     &                     ncal,dt,nt,gammax,
     &                     lstop,kstep,
     &                     ivmax,jvmax,kvmax,ivtrn,jvtrn,kvtrn,
     &                     ivalf,jvalf,kvalf,ivacs,jvacs,kvacs,
     &                     itmin,jtmin,ktmin,
     &                     vmax,valf,vtrn,vacs,
     &                     rr,sain,
     &                     ro1,pg1,ur1,ut1,up1,br1,bt1,bp1,
     &                     ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
       implicit none
* arguments
       intent(in)  ::      ii,jj,kk
       intent(in)  ::      ncal,   nt,gammax
       intent(in)  ::              kstep
       intent(out) ::           dt
       intent(out) ::      lstop
       intent(out) ::      ivmax,jvmax,kvmax,ivtrn,jvtrn,kvtrn
       intent(out) ::      ivalf,jvalf,kvalf,ivacs,jvacs,kvacs
       intent(out) ::      itmin,jtmin,ktmin
       intent(out) ::      vmax,valf,vtrn,vacs
       intent(in)  ::      rr,sain
       intent(out) ::      ro1,pg1,ur1,ut1,up1,br1,bt1,bp1
       intent(in)  ::      ro2,pg2,ur2,ut2,up2,br2,bt2,bp2
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64) ! number of grid for each three direction.
       integer ncal, nt, kstep(0:jj)
       real    dt
       real    gammax
       logical lstop
       integer ivmax, jvmax, kvmax, ivtrn, jvtrn, kvtrn
       integer ivalf, jvalf, kvalf, ivacs, jvacs, kvacs
       integer itmin, jtmin, ktmin
       real    vmax, valf, vtrn, vacs
       real    ro1(-2:jj+2,-1:kk+2,-2:ii+3),pg1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur1(-2:jj+2,-1:kk+2,-2:ii+3),ut1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up1(-2:jj+2,-1:kk+2,-2:ii+3),br1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt1(-2:jj+2,-1:kk+2,-2:ii+3),bp1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    rr(-2:ii+3)
       real    sain(0:jj)
* local
       real    vmax2, valf2, vtrn2, vacs2
       real    valfr, valft, valfp
       real    dt2
       real    aa, bb, cc, vfast
       real    dtr0, dtt0, dtp0
       real    dtra, dtta, dtpa
*       real    dtrv, dttv, dtpv, gradv
       real    dss, ds1, ds2
       integer i, j, k
       integer numnden, numnpre, numneg
* parameter
       real    pi
       parameter(pi = 3.1415926535897932385E+00)
* whether or not check super-alfvenic at upper boundary values
       logical lchkuppr
       parameter(lchkuppr = .true.)
* switches for axi
       logical laxiave
       parameter(laxiave = .true.)

       numnden = 0
       numnpre = 0
       numneg  = 0

!$omp parallel do default(private)
!$omp&            reduction(+:numnden,numnpre,numneg)
!$omp&            shared(ro1,pg1,ur1,ut1,up1,br1,bt1,bp1)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(gammax)
!$omp&            shared(ii,jj,kk)
       do i = -2, ii + 3
         do k =-1, kk + 2
         do j =-2, jj + 2
           ro1(j,k,i) = ro2(j,k,i) ! copy quantities
           pg1(j,k,i) = pg2(j,k,i)
           ur1(j,k,i) = ur2(j,k,i)
           ut1(j,k,i) = ut2(j,k,i)
           up1(j,k,i) = up2(j,k,i)
           br1(j,k,i) = br2(j,k,i)
           bt1(j,k,i) = bt2(j,k,i)
           bp1(j,k,i) = bp2(j,k,i)
           if (i .LE. ii + 1) then
             if (.NOT. (ro2(j,k,i) .GE. -1.0E-30)) then
               numnden = numnden + 1
               write(*,'('' non-positive N  at '',3i4,1x,e11.5)')
     &             i,j,k,ro2(j,k,i)
             endif
             if (.NOT. (pg2(j,k,i) .GE. -1.0E-30)) then
               numnpre = numnpre + 1
               write(*,'('' non-positive P  at '',3i4,1x,e11.5)')
     &             i,j,k,pg2(j,k,i)
             endif
           endif

             if (i .EQ. ii + 1) then
               vacs2 = sqrt(pg2(j,k,i)    / ro2(j,k,i)*gammax)
               valfr = sqrt(br2(j,k,i)**2 / ro2(j,k,i))
               valft = sqrt(bt2(j,k,i)**2 / ro2(j,k,i))
               valfp = sqrt(bp2(j,k,i)**2 / ro2(j,k,i))
               valf2 = sqrt(valfr**2 + valft**2 + valfp**2)
               aa = (valf2**2 + vacs2**2) * 0.5E+00
               bb = aa**2 - vacs2**2 * valfr**2
               if (bb .GT. 0.0E+00) then
                 bb = sqrt(bb)
               else
                 bb = 0.0E+00
               endif
               cc = sqrt(aa + bb)
               cc = cc - ur2(j,k,i)
               if (cc .GE. 0.0E+00) then
                 numneg = numneg + 1
                 write(*,'('' Vr < Vf at '',2i4)') j,k
               endif
             endif

         enddo
         enddo
       enddo ! end of i-loop
!$omp end parallel do

       lstop = .false.
       if (numnden .GT. 0) then
         write(* ,'('' Ro < 0  : Ncal Nt Num ='',3i9)') ncal,nt,numnden
         write(19,'('' Ro < 0  : Ncal Nt Num ='',3i9)') ncal,nt,numnden
         lstop = .true.
       endif
       if (numnpre .GT. 0) then
         write(* ,'('' Pg < 0  : Ncal Nt Num ='',3i9)') ncal,nt,numnpre
         write(19,'('' Pg < 0  : Ncal Nt Num ='',3i9)') ncal,nt,numnpre
         lstop = .true.
       endif
       if (numneg .GT. 0) then
         write(* ,'('' Vr < Vf : Ncal Nt Num ='',3i9)') ncal,nt,numneg
         write(19,'('' Vr < Vf : Ncal Nt Num ='',3i9)') ncal,nt,numneg
         lstop = .true.
       endif
       if (lstop) return

* determine the time step "dt"
       vmax =-1.0E+00
       valf =-1.0E+00
       vtrn =-1.0E+00
       vacs =-1.0E+00
       dt   = 1.0E+30
       itmin = -1
       jtmin = -1
       ktmin = -1
       ivmax = -1
       jvmax = -1
       kvmax = -1
       ivalf = -1
       jvalf = -1
       kvalf = -1
       ivtrn = -1
       jvtrn = -1
       kvtrn = -1
       ivacs = -1
       jvacs = -1
       kvacs = -1

** dare not to include ivalf,,,,etc in neither shared nor private group
!$omp parallel do default(private)
!$omp&            reduction(min:dt)
!$omp&            reduction(max:vmax,valf,vtrn,vacs)
!$omp&            reduction(+:numnden,numnpre,numneg)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(rr,gammax)
!$omp&            shared(kstep,sain)
!$omp&            shared(ii,jj,kk)
       do i = 0, ii + 1
         do k = 1, kk
         do j = 1, jj - 1
* mag.aco.
           vtrn2 = ur2(j,k,i)**2 + ut2(j,k,i)**2 + up2(j,k,i)**2
           vtrn2 = sqrt(vtrn2)
           valfr = sqrt(br2(j,k,i)**2 / ro2(j,k,i))
           valft = sqrt(bt2(j,k,i)**2 / ro2(j,k,i))
           valfp = sqrt(bp2(j,k,i)**2 / ro2(j,k,i))
           valf2 = sqrt(valfr**2 + valft**2 + valfp**2)
           vacs2 = sqrt(pg2(j,k,i) / ro2(j,k,i) * gammax)
           vmax2 = sqrt(vtrn2**2  + valf2**2 + vacs2**2)
           if (valf2 .GT. valf) then
             valf = valf2
             ivalf = i
             jvalf = j
             kvalf = k
           endif
           if (vacs2 .GT. vacs) then
             vacs = vacs2
             ivacs = i
             jvacs = j
             kvacs = k
           endif
           if (vtrn2 .GT. vtrn) then
             vtrn = vtrn2
             ivtrn = i
             jvtrn = j
             kvtrn = k
           endif
           if (vmax2 .GT. vmax) then
             vmax = vmax2
             ivmax = i
             jvmax = j
             kvmax = k
           endif

           aa=(vacs2**2 + valf2**2)**2 - 4.0E+00 * vacs2**2 * valfr**2
           if (aa .LE. 0.0E+00) then
             aa = 0.0E+00
           else
             aa = sqrt(aa)
           endif
           bb = sqrt(0.5E+00 * (vacs2**2 + valf2**2 + aa))
           cc = abs(ur2(j,k,i))
           vfast = bb + cc
*           vfast = vmax2
           if (i .EQ. 0) then
             dss = rr(1) - rr(0)
           else if (i .EQ. ii + 1) then
             dss = rr(ii+1) - rr(ii)
           else
             ds1 = rr(i+1) - rr(i)
             ds2 = rr(i) - rr(i-1)
             dss = min(ds1,ds2)
           endif
           dtra  = vfast / dss
           dtr0  = 1.0E+00 / dtra
**
           aa=(vacs2**2 + valf2**2)**2 - 4.0E+00 * vacs2**2 * valft**2
           if (aa .LE. 0.0E+00) then
             aa = 0.0E+00
           else
             aa = sqrt(aa)
           endif
           bb = sqrt(0.5E+00 * (vacs2**2 + valf2**2 + aa))
           cc = abs(ut2(j,k,i))
           vfast = bb + cc
           dss   = pi / float(jj) * rr(i)
           dtta  = vfast / dss
           dtt0  = 1.0E+00 / dtta
**
           aa=(vacs2**2 + valf2**2)**2 - 4.0E+00 * vacs2**2 * valfp**2
           if (aa .LE. 0.0E+00) then
             aa = 0.0E+00
           else
             aa = sqrt(aa)
           endif
           bb = sqrt(0.5E+00 * (vacs2**2 + valf2**2 + aa))
           cc = abs(up2(j,k,i))
           vfast = bb + cc
           if ((j .EQ. 0) .OR. (j .EQ. jj)) then ! MIND : this is always false in this code
             dss =           pi / float(kk) * rr(i) * sain(1)
           else
             dss = 2.0E+00 * pi / float(kk) * rr(i) * sain(j)
             if ((kk .GT. 3) .AND. (laxiave)) then
               if (j .EQ. 1) then
                 bb = float(kstep(j))
                 cc = float(kstep(j+1))
                 dss = dss * min(bb,cc)
               else if (j .EQ. jj - 1) then
                 aa = float(kstep(j-1))
                 bb = float(kstep(j))
                 dss = dss * min(aa,bb)
               else
                 aa = float(kstep(j-1))
                 bb = float(kstep(j))
                 cc = float(kstep(j+1))
                 dss = dss * min(aa,bb,cc)
               endif
             endif
           endif
           dtpa  = vfast / dss
           dtp0  = 1.0E+00 / dtpa
*
           dt2 = min(dtr0, dtt0, dtp0)
           if (dt2 .LT. dt) then
             dt = dt2
             itmin = i
             jtmin = j
             ktmin = k
           endif
         enddo
         enddo
       enddo ! end i-loop
!$omp end parallel do

       return
       endsubroutine

*
*---------------------------------------------------------------------
*
* I) Boundary conditions
*
* 1) Rotational axis
*
* 2) Outflow boundary
*   a) The deviative d/dr is obtained with the manner
*        "up-wind differencing" of first- or second-order.
*   b) Variables of the next time step is obtained by one-step
*         calclulation which is similar to the first-step
*         in two-step method adapted in the domain of interest.
*
* 3) Inner boundary
*   a) basic parts
*   b) compatibility relation method is NOT treated here.
*
* II) Give perturbation
*    Give perturbation as surface activities which will cause
*       such as CME, spray ejection, and other phenomenon.
*    The location, including the height, where the increase(or decrease) of
*       temperature occur can be arbitrarily set
*       by editing the 'parameter' line in this sub-routine.
*
*---------------------------------------------------------------------
*
       subroutine surcnd(ii,jj,kk,
     &                   lpert,nstepert,cijkpert,nstep,
     &                   omega0,omega,shiftlng,
     &                   iprc,gammax,gamma0,v0,n0,p0,r0,vcrtrn,
     &                   ro0,pg0,ur0,ut0,up0,br0,bt0,bp0,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                   ed2,mr2,mt2,mp2,
     &                   rr,phi,sain,kosa,sainb,kosab,
     &                   sinph,cosph,ropark0,pgpark0)
       implicit none
* arguments
       intent(in)    ::  ii,jj,kk
       intent(inout) ::  lpert
       intent(in)    ::        nstepert,cijkpert,nstep
       intent(in)    ::  omega0,omega,shiftlng
       intent(in)    ::  iprc,gammax,gamma0,v0,n0,p0,r0,vcrtrn
       intent(in)    ::  ro0,pg0,ur0,ut0,up0,br0,bt0,bp0
       intent(inout) ::  ro2,pg2,ur2,ut2,up2,br2,bt2,bp2
       intent(inout) ::  ed2,mr2,mt2,mp2
       intent(in)    ::  rr,phi,sain,kosa,sainb,kosab
       intent(in)    ::  sinph,cosph,ropark0,pgpark0
       logical lpert
       integer nstepert(2), cijkpert(3), nstep, iprc
       real    gammax, gamma0, v0, n0, p0, r0
       real    omega, omega0, shiftlng
       real    vcrtrn
* large arrays
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ed2(-2:jj+2,-1:kk+2,-2:ii+3),mr2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    mt2(-2:jj+2,-1:kk+2,-2:ii+3),mp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    pgpark0(-2:ii+3), ropark0(-2:ii+3)
       integer ipho
       parameter(ipho = 5)
       real    ro0(0:jj,1:kk,-1:ipho),pg0(0:jj,1:kk,-1:ipho)
       real    ur0(0:jj,1:kk,-1:ipho),ut0(0:jj,1:kk,-1:ipho)
       real    up0(0:jj,1:kk,-1:ipho),br0(0:jj,1:kk,-1:ipho)
       real    bt0(0:jj,1:kk,-1:ipho),bp0(0:jj,1:kk,-1:ipho)
       real    rr(-2:ii+3),phi(-1:kk+2)
       real    sain(0:jj), kosa(0:jj), sainb(0:jj-1), kosab(0:jj-1)
       real    sinph(0:kk+1), cosph(0:kk+1)
* choice of compati.
       integer ichrsurf
       parameter(ichrsurf = 0) ! this may appear at other parts.
*
       integer irkstep ! def. Runge-Kutta step
       parameter(irkstep = 1)
* mass flux limits in terms of magnetic field strength etc. : This will be used to host N and T maps from lower coronal model
       logical lmflxvar
       parameter(lmflxvar = .false.)
* smallest possible Mr.
       real     mrzero
       parameter(mrzero = 1.0E-00) ! mind mrzero = 1 / rr**2 or same order.
* criteria of|Br/B|
       real    brabsb
       parameter(brabsb = 1.0E-01)
* impose gradual increase of Br at surface
!       logical ldmagon
!       parameter(ldmagon = .false.)
* local
       integer i, j, k
       real    aa,bb,cc
       real    fext(4), fext2(4) ! lower & upper
       integer ibpower
       parameter(ibpower = 1) ! B r^iblpower for extrapolation
       real    pi
       parameter(pi = 3.1415926535897932385E+00)

* surpress warnings, silly
       lpert = .false.
       i  = iprc
       i  = nstepert(1)
       i  = cijkpert(1)
       i  = nstep
       aa = omega
       aa = omega0
       aa = n0
       aa = r0
       aa = p0
       aa = v0
       aa = vcrtrn
       aa = gamma0
       aa = shiftlng
       aa = phi(1)
       aa = kosa(1)
       aa = sain(1)
       aa = sainb(1)
       aa = kosab(1)
       aa = cosph(1)
       aa = sinph(1)
       aa = sinph(1)
       aa = ropark0(1)
       aa = pgpark0(1)
       aa = ro0(1,1,1)
       aa = pg0(1,1,1)
       aa = ur0(1,1,1)

!$omp parallel do private(i,k,j,aa,bb,cc)
         do i = -2, ii + 3
           if (i .EQ. 0) then

             if (ichrsurf .EQ. -1) then

               do k = 1, kk
               do j = 0, jj
                 br2(j,k,i) = br0(j,k,i) ! unchanged
               enddo
               enddo

             else

               do k = 1, kk
               do j = 0, jj
* Vr
                 if ((.NOT. (lmflxvar)) .AND.
     &               (ichrsurf .NE. 1) .AND.
     &               (ichrsurf .NE. 5) .AND.
     &               (ichrsurf .NE. 6) .AND.
     &               (ur2(j,k,i) .GT. vcrtrn)) ur2(j,k,i) = vcrtrn ! this occur at polar averaging.

                 if  (ro2(j,k,i)*ur2(j,k,i) .LT. mrzero)
     &                                 ur2(j,k,i) = 0.0E+00

                 br2(j,k,i) = br0(j,k,i) ! unchanged
** for v//B

                   cc =(br0(j,k,i)**2+bt0(j,k,i)**2+bp0(j,k,i)**2)
     &                / pg0(j,k,i)
                   if (cc .GT. 1.0E-05) then ! most MHD case
                     aa = br2(j,k,i)**2+bt2(j,k,i)**2+bp2(j,k,i)**2
                     bb = abs(br2(j,k,i)) / sqrt(aa+1.0E-10)
                     if (bb .GT. brabsb) then
                       if (ro2(j,k,i)*ur2(j,k,i) .LT. mrzero) then
                         ur2(j,k,i) = 0.0E+00
                         ut2(j,k,i) = 0.0E+00
                         up2(j,k,i) = 0.0E+00
                       else
                         ut2(j,k,i) = ur2(j,k,i)*bt2(j,k,i)/br2(j,k,i)
                         up2(j,k,i) = ur2(j,k,i)*bp2(j,k,i)/br2(j,k,i)
                       endif
                     else
                       ur2(j,k,i) = 0.0E+00
                       ut2(j,k,i) = 0.0E+00
                       up2(j,k,i) = 0.0E+00
                     endif
                   else ! maybe, only case to test Parker solution or something like.
*                     ur2(j,k,i) = ur0(j,k,i)
                     ut2(j,k,i) = ut0(j,k,i)
                     up2(j,k,i) = up0(j,k,i)
                     bt2(j,k,i) = bt0(j,k,i)
                     bp2(j,k,i) = bp0(j,k,i)
                   endif
               enddo
               enddo
             endif ! end if (ichrsurf=-1) or not
           endif
         enddo ! end-i-loop
!$omp end parallel do

* conservatives...
!$omp parallel do private(i,k,j)
       do i = -2, ii + 3
         do k = -1, kk + 2
         do j = -2, jj + 2
           mr2(j,k,i) = ur2(j,k,i) * ro2(j,k,i)
           mt2(j,k,i) = ut2(j,k,i) * ro2(j,k,i)
           mp2(j,k,i) = up2(j,k,i) * ro2(j,k,i)
           ed2(j,k,i) = pg2(j,k,i) /(gammax - 1.0E+00)
     &                +(ur2(j,k,i)**2 + ut2(j,k,i)**2
     &                 +up2(j,k,i)**2)* ro2(j,k,i) * 0.5E+00
     &                - 2.0E+00 * ro2(j,k,i) / rr(i)
         enddo
         enddo
       enddo
!$omp end parallel do

         fext(1) =  3.0E+00 ! usual 2nd order .... possibly best for most cases
         fext(2) = -3.0E+00
         fext(3) =  1.0E+00
         fext(4) =  0.0E+00

** extrapolation for ghost i > ii + 1
         fext2(1) =  3.0E+00 ! usual 2nd order
         fext2(2) = -3.0E+00
         fext2(3) =  1.0E+00
         fext2(4) =  0.0E+00
**

!$omp parallel do private(i,k,j,aa,bb)
       do i = -2, ii + 3
         if (i .EQ. -1) then
             do k = 1, kk
             do j = 0, jj
* magnetic field...
                 br2(j,k,i) =(br2(j,k,i+1)*rr(i+1)**ibpower * fext(1)
     &                      + br2(j,k,i+2)*rr(i+2)**ibpower * fext(2)
     &                      + br2(j,k,i+3)*rr(i+3)**ibpower * fext(3)
     &                      + br2(j,k,i+4)*rr(i+4)**ibpower * fext(4))
     &                      / rr(i)**ibpower
                 bt2(j,k,i) =(bt2(j,k,i+1)*rr(i+1)**ibpower * fext(1)
     &                      + bt2(j,k,i+2)*rr(i+2)**ibpower * fext(2)
     &                      + bt2(j,k,i+3)*rr(i+3)**ibpower * fext(3)
     &                      + bt2(j,k,i+4)*rr(i+4)**ibpower * fext(4))
     &                      / rr(i)**ibpower
                 bp2(j,k,i) =(bp2(j,k,i+1)*rr(i+1)**ibpower * fext(1)
     &                      + bp2(j,k,i+2)*rr(i+2)**ibpower * fext(2)
     &                      + bp2(j,k,i+3)*rr(i+3)**ibpower * fext(3)
     &                      + bp2(j,k,i+4)*rr(i+4)**ibpower * fext(4))
     &                      / rr(i)**ibpower
* momentum
               mr2(j,k,i) =(mr2(j,k,i+1)*rr(i+1)**2 * fext(1)
     &                    + mr2(j,k,i+2)*rr(i+2)**2 * fext(2)
     &                    + mr2(j,k,i+3)*rr(i+3)**2 * fext(3)
     &                    + mr2(j,k,i+4)*rr(i+4)**2 * fext(4))
     &                    / rr(i)**2
               mt2(j,k,i) =(mt2(j,k,i+1)*rr(i+1)**2 * fext(1)
     &                    + mt2(j,k,i+2)*rr(i+2)**2 * fext(2)
     &                    + mt2(j,k,i+3)*rr(i+3)**2 * fext(3)
     &                    + mt2(j,k,i+4)*rr(i+4)**2 * fext(4))
     &                    / rr(i)**2
               mp2(j,k,i) =(mp2(j,k,i+1)*rr(i+1)**2 * fext(1)
     &                    + mp2(j,k,i+2)*rr(i+2)**2 * fext(2)
     &                    + mp2(j,k,i+3)*rr(i+3)**2 * fext(3)
     &                    + mp2(j,k,i+4)*rr(i+4)**2 * fext(4))
     &                    / rr(i)**2

* gas
                 aa = log(ro2(j,k,i+1)) * fext(1)
     &              + log(ro2(j,k,i+2)) * fext(2)
     &              + log(ro2(j,k,i+3)) * fext(3)
     &              + log(ro2(j,k,i+4)) * fext(4)
                 ro2(j,k,i) = exp(aa)
                 bb = log(pg2(j,k,i+1)) * fext(1)
     &              + log(pg2(j,k,i+2)) * fext(2)
     &              + log(pg2(j,k,i+3)) * fext(3)
     &              + log(pg2(j,k,i+4)) * fext(4)
                 pg2(j,k,i) = exp(bb)
* flow velocity
                     ur2(j,k,i) = mr2(j,k,i) / ro2(j,k,i)
                     ut2(j,k,i) = mt2(j,k,i) / ro2(j,k,i)
                     up2(j,k,i) = mp2(j,k,i) / ro2(j,k,i)
** test
!                   mr2(j,k,i) = 0.0E+00
!                   mt2(j,k,i) = 0.0E+00
!                   mp2(j,k,i) = 0.0E+00
!                   ur2(j,k,i) = 0.0E+00
!                   ut2(j,k,i) = 0.0E+00
!                   up2(j,k,i) = 0.0E+00
**
             enddo ! j,k-loop
             enddo

         endif ! end-if (i .EQ. -1) ! to calculate ghost at i = -1

         if (i .EQ. ii + 2) then
             do k = 1, kk
             do j = 0, jj
               ro2(j,k,i) =(ro2(j,k,i-1)*rr(i-1)**2 * fext2(1)
     &                    + ro2(j,k,i-2)*rr(i-2)**2 * fext2(2)
     &                    + ro2(j,k,i-3)*rr(i-3)**2 * fext2(3)
     &                    + ro2(j,k,i-4)*rr(i-4)**2 * fext2(4))
     &                    / rr(i)**2
               pg2(j,k,i) =(pg2(j,k,i-1)*rr(i-1)**2 * fext2(1)
     &                    + pg2(j,k,i-2)*rr(i-2)**2 * fext2(2)
     &                    + pg2(j,k,i-3)*rr(i-3)**2 * fext2(3)
     &                    + pg2(j,k,i-4)*rr(i-4)**2 * fext2(4))
     &                    / rr(i)**2

               mr2(j,k,i) =(mr2(j,k,i-1)*rr(i-1)**2 * fext2(1)
     &                    + mr2(j,k,i-2)*rr(i-2)**2 * fext2(2)
     &                    + mr2(j,k,i-3)*rr(i-3)**2 * fext2(3)
     &                    + mr2(j,k,i-4)*rr(i-4)**2 * fext2(4))
     &                    / rr(i)**2
               mt2(j,k,i) =(mt2(j,k,i-1)*rr(i-1)**2 * fext2(1)
     &                    + mt2(j,k,i-2)*rr(i-2)**2 * fext2(2)
     &                    + mt2(j,k,i-3)*rr(i-3)**2 * fext2(3)
     &                    + mt2(j,k,i-4)*rr(i-4)**2 * fext2(4))
     &                    / rr(i)**2
               mp2(j,k,i) =(mp2(j,k,i-1)*rr(i-1)**2 * fext2(1)
     &                    + mp2(j,k,i-2)*rr(i-2)**2 * fext2(2)
     &                    + mp2(j,k,i-3)*rr(i-3)**2 * fext2(3)
     &                    + mp2(j,k,i-4)*rr(i-4)**2 * fext2(4))
     &                    / rr(i)**2
               ur2(j,k,i) = mr2(j,k,i) / ro2(j,k,i)
               ut2(j,k,i) = mt2(j,k,i) / ro2(j,k,i)
               up2(j,k,i) = mp2(j,k,i) / ro2(j,k,i)
               br2(j,k,i) =(br2(j,k,i-1)*rr(i-1)    * fext2(1)
     &                    + br2(j,k,i-2)*rr(i-2)    * fext2(2)
     &                    + br2(j,k,i-3)*rr(i-3)    * fext2(3)
     &                    + br2(j,k,i-4)*rr(i-4)    * fext2(4))
     &                    / rr(i)
               bt2(j,k,i) =(bt2(j,k,i-1)*rr(i-1)    * fext2(1)
     &                    + bt2(j,k,i-2)*rr(i-2)    * fext2(2)
     &                    + bt2(j,k,i-3)*rr(i-3)    * fext2(3)
     &                    + bt2(j,k,i-4)*rr(i-4)    * fext2(4))
     &                    / rr(i)
               bp2(j,k,i) =(bp2(j,k,i-1)*rr(i-1)    * fext2(1)
     &                    + bp2(j,k,i-2)*rr(i-2)    * fext2(2)
     &                    + bp2(j,k,i-3)*rr(i-3)    * fext2(3)
     &                    + bp2(j,k,i-4)*rr(i-4)    * fext2(4))
     &                    / rr(i)
             enddo
             enddo

         endif ! at I + 2
       enddo
!$omp end parallel do

!$omp parallel do private(i,k,j,aa,bb)
       do i = -2, ii + 3
         if (i .EQ. -2) then ! xx2(j,k,-2) may not be used
           aa =(rr(i+1)/rr(i))**2
           bb = rr(i+1)/rr(i)
           do k = 1, kk
           do j = 0, jj
             ro2(j,k,i) = ro2(j,k,i+1) * aa
             pg2(j,k,i) = pg2(j,k,i+1) * aa
             ur2(j,k,i) = ur2(j,k,i+1)
             ut2(j,k,i) = ut2(j,k,i+1)
             up2(j,k,i) = up2(j,k,i+1)
             br2(j,k,i) = br2(j,k,i+1) * bb
             bt2(j,k,i) = bt2(j,k,i+1) * bb
             bp2(j,k,i) = bp2(j,k,i+1) * bb
           enddo
           enddo
         endif
         if (i .EQ. ii+3) then ! xx2(j,k,ii+3) may not be used
           aa =(rr(ii+2)/rr(i))**2
           bb = rr(ii+2)/rr(i)
           do k = 1, kk
           do j = 0, jj
             ro2(j,k,i) = ro2(j,k,ii+2) * aa
             pg2(j,k,i) = pg2(j,k,ii+2) * aa
             ur2(j,k,i) = ur2(j,k,ii+2)
             ut2(j,k,i) = ut2(j,k,ii+2)
             up2(j,k,i) = up2(j,k,ii+2)
             br2(j,k,i) = br2(j,k,ii+2) * bb
             bt2(j,k,i) = bt2(j,k,ii+2) * bb
             bp2(j,k,i) = bp2(j,k,ii+2) * bb
           enddo
           enddo
         endif
       enddo
!$omp end parallel do

       return
       endsubroutine

*
* --------------------------------------------------------------------
*
* to paritally-solve Poison eq. for divB-cleaning
*
* --------------------------------------------------------------------
*
       subroutine divbcln(ii,jj,kk,
     &                    gammax,nt,ncal,dt,iprc,kstep,
     &                    pg2,ur2,ut2,up2,br2,bt2,bp2,br3,bt3,bp3,
     &                    rr,theta,phi,sain,kosa,sainb,
     &                    divbpsi1,divbpsi2,divmag)
       implicit none
* arguments
       intent(in)    ::   ii,jj,kk
       intent(in)    ::   gammax,nt,ncal,dt,iprc,kstep
       intent(inout) ::   pg2,            br2,bt2,bp2,br3,bt3,bp3
       intent(in)    ::       ur2,ut2,up2
       intent(in)    ::   rr,theta,phi,sain,kosa,sainb
       intent(out)   ::   divbpsi1,divbpsi2,divmag
       real    gammax
       integer ncal, nt, iprc
       real    dt
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       integer  kstep(0:jj)
       real    rr(-2:ii+3), theta(0:jj), phi(-1:kk+2)
       real    sain(0:jj), kosa(0:jj), sainb(0:jj-1)
       real    pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ut2(-2:jj+2,-1:kk+2,-2:ii+3),up2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    br3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt3(-2:jj+2,-1:kk+2,-2:ii+3),bp3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    divbpsi1(-1:jj,0:kk+1,-1:ii+1)
       real    divbpsi2(-1:jj,0:kk+1,-1:ii+1)
       real    divmag(-1:jj,0:kk+1,-1:ii+1)
* local
       integer i, j, k, k2, j2, i2, kstep2, ittime, idummy
       real    aa, bb, cc, dd
       real    rhs ! , lhs
       real    dth, dph, dtlocal, dtclean
       real    sumsum, sumwei, weiwei
       real    rrb(-1:ii+1)
       integer kstepb(0:jj)
* if distant mode or dmag mode
       logical ldistant
       parameter(ldistant = .false.)
       logical ldmagon
       parameter(ldmagon = .false.)
* R-K
       integer irkstep
       parameter(irkstep = 1)
* choice
       logical ladjstpg
       parameter(ladjstpg = .false.) ! whether adjust Pg after modifying B
       integer idivbflt !        if not 1 nor 2, lnoavebp may be false
       parameter(idivbflt = 0) ! 0=new partial relax, 1=old p.relax, 2=hype-para, 3=modif h-p(test), other=test
       logical lreset0 !         reset values when idivbflt is 0
       parameter(lreset0 = .true.) ! MIND that if false, the extension calculation may be slightly different.
* variables for idivbflt equal to 1
       integer itmax, itmax0
       parameter(itmax0 = 2) ! how many times it will be tried. 1=actually do nothing, 2=only once, 4 is frequent.
* variables for idivbflt equal to 2
       real    ch2, cp2 ! square of velocities (h & p = hyperbolic and parabolic, or propagation and
       parameter(ch2 = 0.1E+00, cp2 = 0.5E+00) ! order of unity may be better, for safty.
* if B_phi will NOT be smoothed, even when laxiave = ON
       logical lnoavebp
       parameter(lnoavebp = .true.) ! MIND : true may not work well when idivflt is not 1 nor 2
* flag whether or not i = 0 will be modified.
       logical ldivbext
       parameter(ldivbext = .false.) ! whether Br at ghost cells are determined from divB=0
* const.
       real    pi
       parameter(pi = 3.1415926535897932385E+00)

       do j = 0, jj
*         kstepb(j) = 1  ! test,
         kstepb(j) = kstep(j)
       enddo

       aa = gammax

       dth = pi / float(jj) ! theta(1) - theta(0)
       dph = pi / float(kk) * 2.0E+00 ! phi(1) - phi(0)

!$omp parallel do private(i)
       do i = -1, ii + 1
         if (i .EQ. ii+1) then
           rrb(i) = rr(i) + (rr(i) - rr(i-1)) * 0.5E+00
         else if (i .EQ. -1) then
           rrb(i) = rr(i+1) - (rr(i+2) - rr(i+1)) * 0.5E+00
         else
           rrb(i) = (rr(i+1) + rr(i)) * 0.5E+00
         endif
       enddo
!$omp end parallel do

* refresh .... needed ???
!$omp parallel do private(i,k,j)
       do i = -1, ii + 1
       do k =  0, kk + 1
       do j = -1, jj
         divmag(j,k,i) = 0.0E+00
       enddo
       enddo
       enddo
!$omp end parallel do

* store divB
!$omp parallel do private(i,k,j,aa)
       do i =  0, ii
         aa  = (rr(i+1)**2 - rr(i)**2)
     &       / (rr(i+1)**3 - rr(i)**3) * 1.5E+00 ! 1/r
         do k = 1, kk
         do j = 0, jj - 1
           divmag(j,k,i)
     &         =(((br2(j  ,k,i+1) + br2(j  ,k+1,i+1)
     &            +br2(j+1,k,i+1) + br2(j+1,k+1,i+1)) * rr(i+1)**2
     &           -(br2(j  ,k,i  ) + br2(j  ,k+1,i  )
     &            +br2(j+1,k,i  ) + br2(j+1,k+1,i  )) * rr(i  )**2)
     &          / (rr(i+1)**3 - rr(i)**3) * 3.0E+00
     &          +((bt2(j+1,k  ,i) + bt2(j+1,k  ,i+1)
     &            +bt2(j+1,k+1,i) + bt2(j+1,k+1,i+1)) * sain(j+1)
     &           -(bt2(j  ,k  ,i) + bt2(j  ,k  ,i+1)
     &            +bt2(j  ,k+1,i) + bt2(j  ,k+1,i+1)) * sain(j  ))
     &          * aa / (kosa(j) - kosa(j+1))
     &          +((bp2(j,k+1,i+1) + bp2(j+1,k+1,i+1)
     &            +bp2(j,k+1,i  ) + bp2(j+1,k+1,i  ))
     &           -(bp2(j,k  ,i+1) + bp2(j+1,k  ,i+1)
     &            +bp2(j,k  ,i  ) + bp2(j+1,k  ,i  )))
     &          * aa / (kosa(j) - kosa(j+1)) * dth / dph)
     &         * 0.25E+00
         enddo
         enddo
       enddo
!$omp end parallel do

!$omp parallel do private(i,k,k2)
       do i = 0, ii
         do k = 1, kk
           k2 = mod(k + kk/2,kk)
           if (k2 .EQ. 0) k2 = kk
           divmag(  -1,k,i) = divmag(   0,k2,i)
           divmag(jj  ,k,i) = divmag(jj-1,k2,i)
         enddo
       enddo
!$omp end parallel do

!$omp parallel do private(i,j)
       do i = 0, ii
         do j = -1, jj
           divmag(j,   0,i) = divmag(j,kk,i)
           divmag(j,kk+1,i) = divmag(j, 1,i)
         enddo
       enddo
!$omp end parallel do

!$omp parallel do private(i,k,j)
       do i = -1, ii + 1
         if (i .EQ. ii + 1) then
           do k =  1, kk
           do j =  0, jj - 1
             divmag(j,k,i) = divmag(j,k,i-1)
           enddo
           enddo
         endif
         if (i .EQ. -1) then
           do k =  1, kk
           do j =  0, jj - 1
             divmag(j,k,i) = divmag(j,k,i+1)
           enddo
           enddo
         endif
       enddo
!$omp end parallel do

!       if (idivbflt .EQ. 0) then ! Gauss-type time-relax ; def. of grad and div identical to other part

* once reset...... (needed or not ?)
         if (lreset0) then
!$omp parallel do private(i,k,j)
           do i = -1, ii + 1
           do k =  0, kk + 1
           do j = -1, jj
             divbpsi1(j,k,i) = 0.0E+00 ! initialize
             divbpsi2(j,k,i) = 0.0E+00 ! initialize
           enddo
           enddo
           enddo
!$omp end parallel do
         else if (ncal .EQ. 0) then
!$omp parallel do private(i,k,j)
           do i = -1, ii + 1
           do k =  0, kk + 1
           do j = -1, jj
             divbpsi1(j,k,i) = 0.0E+00 ! initialize
             divbpsi2(j,k,i) = 0.0E+00 ! initialize
           enddo
           enddo
           enddo
!$omp end parallel do
         endif

* local time-step ..... supposing the smallest grid locate at innermost & pole axis
         aa = rr(1) - rr(0)
         bb = rrb(0) * dth
         cc = rrb(0) * sainb(0) * dph
         dd = min(aa,bb,cc)
         dtlocal = dd**2 * 0.5E+00 * 0.8E+00 ! CFL

** simplified and generally better. Mind do not redefine dtlocal here
         if ((nt .EQ. 0) .AND. (ncal .EQ. 0)) then
           itmax = max(jj,kk)
         else
           itmax = itmax0
         endif

         if (ncal .EQ. 1) then
           write(*,*) 'Grad(divB) --- '
           write(*,*) 'Ntime dt_local(original) dt_global = '
           write(*,*)  itmax, dtlocal, dt
         endif

* iterration process
         do ittime = 1, itmax
           dtclean = 1.0E+00
!$omp parallel do private(i,k,j)
           do i = 0, ii + 1
             if (.NOT. ((ldistant) .AND. (i .EQ. 0))) then
               do k = 1, kk
               do j = 0, jj ! 1, jj - 1
                 if ((j .EQ. 0) .OR. (j .EQ. jj)) then ! rotational axis
                   br3(j,k,i) =
     &                   ((divbpsi1(j,k  ,i  )+divbpsi1(j-1,k  ,i  )
     &                    +divbpsi1(j,k-1,i  )+divbpsi1(j-1,k-1,i  ))
     &                   -(divbpsi1(j,k  ,i-1)+divbpsi1(j-1,k  ,i-1)
     &                    +divbpsi1(j,k-1,i-1)+divbpsi1(j-1,k-1,i-1)))
     &                  / (rrb(i) - rrb(i-1)) * 0.25E+00 * dtclean
                   bt3(j,k,i) = 0.0E+00 ! this is temporal
                   bp3(j,k,i) = 0.0E+00 ! this is temporal
                 else
                   br3(j,k,i) =
     &                    ((divbpsi1(j,k  ,i  )+divbpsi1(j-1,k  ,i  )
     &                     +divbpsi1(j,k-1,i  )+divbpsi1(j-1,k-1,i  ))
     &                    -(divbpsi1(j,k  ,i-1)+divbpsi1(j-1,k  ,i-1)
     &                     +divbpsi1(j,k-1,i-1)+divbpsi1(j-1,k-1,i-1)))
     &                   / (rrb(i) - rrb(i-1)) * 0.25E+00 * dtclean
                   bt3(j,k,i) =
     &                   (((divbpsi1(j  ,k,i  )+divbpsi1(j  ,k-1,i  ))
     &                    -(divbpsi1(j-1,k,i  )+divbpsi1(j-1,k-1,i  )))
     &                    / rrb(i  )
     &                  + ((divbpsi1(j  ,k,i-1)+divbpsi1(j  ,k-1,i-1))
     &                    -(divbpsi1(j-1,k,i-1)+divbpsi1(j-1,k-1,i-1)))
     &                    / rrb(i-1))
     &                   / (4.0E+00 * dth) * dtclean
                   bp3(j,k,i) =
     &                    ((divbpsi1(j  ,k  ,i  )
     &                     -divbpsi1(j  ,k-1,i  ))/rrb(i  )/sainb(j  )
     &                  +  (divbpsi1(j-1,k  ,i  )
     &                     -divbpsi1(j-1,k-1,i  ))/rrb(i  )/sainb(j-1)
     &                  +  (divbpsi1(j  ,k  ,i-1)
     &                     -divbpsi1(j  ,k-1,i-1))/rrb(i-1)/sainb(j  )
     &                  +  (divbpsi1(j-1,k  ,i-1)
     &                     -divbpsi1(j-1,k-1,i-1))/rrb(i-1)/sainb(j-1))
     &                   / (4.0E+00 * dph) * dtclean
                 endif
               enddo
               enddo
* overlap along longitudinal direction
               do j = 0, jj
                 br3(j,   0,i) = br3(j,kk,i)
                 br3(j,kk+1,i) = br3(j, 1,i)
                 bt3(j,   0,i) = bt3(j,kk,i)
                 bt3(j,kk+1,i) = bt3(j, 1,i)
                 bp3(j,   0,i) = bp3(j,kk,i)
                 bp3(j,kk+1,i) = bp3(j, 1,i)
               enddo
               do j = 0, jj
                 br3(j,  -1,i) = br3(j,kk-1,i)
                 br3(j,kk+2,i) = br3(j,   2,i)
                 bt3(j,  -1,i) = bt3(j,kk-1,i)
                 bt3(j,kk+2,i) = bt3(j,   2,i)
                 bp3(j,  -1,i) = bp3(j,kk-1,i)
                 bp3(j,kk+2,i) = bp3(j,   2,i)
               enddo
             endif
* adjust & make sure ....
             if (i .EQ. 0) then
               do k = 0, kk + 1
               do j = 0, jj
                 br3(j,k,i) = 0.0E+00
               enddo
               enddo
             endif
           enddo
!$omp end parallel do

!$omp parallel do private(i,k,j,rhs,aa)
           do i = -1, ii + 1
             if ((i .GE. 0) .AND. (i .LT. ii + 1)) then
               aa  = (rr(i+1)**2 - rr(i)**2)
     &             / (rr(i+1)**3 - rr(i)**3) * 1.5E+00 ! 1/r
               do k =  1, kk
               do j =  0, jj - 1
                 rhs ! should be same as the definition of divB_2
     &            =(((br3(j  ,k,i+1) + br3(j  ,k+1,i+1)
     &               +br3(j+1,k,i+1) + br3(j+1,k+1,i+1)) * rr(i+1)**2
     &              -(br3(j  ,k,i  ) + br3(j  ,k+1,i  )
     &               +br3(j+1,k,i  ) + br3(j+1,k+1,i  )) * rr(i  )**2)
     &             / (rr(i+1)**3 - rr(i)**3) * 3.0E+00
     &             +((bt3(j+1,k  ,i) + bt3(j+1,k  ,i+1)
     &               +bt3(j+1,k+1,i) + bt3(j+1,k+1,i+1)) * sain(j+1)
     &              -(bt3(j  ,k  ,i) + bt3(j  ,k  ,i+1)
     &               +bt3(j  ,k+1,i) + bt3(j  ,k+1,i+1)) * sain(j  ))
     &             * aa / (kosa(j) - kosa(j+1))
     &             +((bp3(j,k+1,i+1) + bp3(j+1,k+1,i+1)
     &               +bp3(j,k+1,i  ) + bp3(j+1,k+1,i  ))
     &              -(bp3(j,k  ,i+1) + bp3(j+1,k  ,i+1)
     &               +bp3(j,k  ,i  ) + bp3(j+1,k  ,i  )))
     &             * aa / (kosa(j) - kosa(j+1)) * dth / dph)
     &            * 0.25E+00
                 divbpsi2(j,k,i) = divbpsi1(j,k,i)
     &                           +(rhs - divmag(j,k,i)) * dtlocal ! Time-relax
               enddo
               enddo
             endif
           enddo
!$omp end parallel do

!$omp parallel do private(i,k,j,bb,cc)
           do i = -1, ii + 1
             if (i .EQ. ii + 1) then
               do k =  1, kk
               do j =  0, jj - 1
                 if ((ldistant) .OR. (ldmagon)) then
                   divbpsi2(j,k,i) = 0.0E+00
                 else
                   divbpsi2(j,k,i) = divbpsi2(j,k,i-1)
                 endif
               enddo
               enddo
             endif
!             if (ldivbext) then
!             else
               if (i .EQ. -1) then
                 do k =  1, kk
                 do j =  0, jj - 1
                   divbpsi2(j,k,i) = divbpsi2(j,k,i+1)
                 enddo
                 enddo
               endif
!             endif
           enddo
!$omp end parallel do

* around axi : ave. along lon.
           if (kk .GT. 3) then
!             if (.NOT. lnoavebp) then
!             else

!$omp parallel do private(i,k,aa,bb) ! old version
               do i = -1, ii + 1
                 aa = 0.0E+00
                 bb = 0.0E+00
                 do k = 1, kk
                   aa = aa + divbpsi2(   0,k,i)
                   bb = bb + divbpsi2(jj-1,k,i)
                 enddo
                 aa = aa / float(kk)
                 bb = bb / float(kk)
                 do k = 1, kk
                   divbpsi2(   0,k,i) = aa
                   divbpsi2(jj-1,k,i) = bb
                 enddo
               enddo
!$omp end parallel do
!             endif ! end if (.NOT. lnoavebp) then

           endif ! end if (kk .GT. 3)...

* rotational axis
!$omp parallel do private(i,k,aa,bb,k2)
           do i = -1, ii + 1
             do k = 1, kk
               k2 = mod(k + kk/2,kk)
               if (k2 .EQ. 0) k2 = kk
               divbpsi2(  -1,k,i) = divbpsi2(   0,k2,i)
               divbpsi2(jj  ,k,i) = divbpsi2(jj-1,k2,i)
             enddo
           enddo
!$omp end parallel do

* overlap along longitude
!$omp parallel do private(i,j)
           do i = -1, ii + 1
             do j = -1, jj
               divbpsi2(j,   0,i) = divbpsi2(j,kk,i)
               divbpsi2(j,kk+1,i) = divbpsi2(j, 1,i)
             enddo
           enddo
!$omp end parallel do

* copy
!$omp parallel do private(i,k,j)
           do i = -1, ii + 1
             do k =  0, kk + 1
             do j = -1, jj
               divbpsi1(j,k,i) = divbpsi2(j,k,i)
             enddo
             enddo
           enddo
!$omp end parallel do

         enddo ! END of DO-iteration loop

         dtclean = 1.0E+00

!       else if (idivbflt .EQ. 1) then ! Gauss - and - Time relax. old but a little bit faster; grad div psi....
!       else if (idivbflt .EQ. 2) then ! telegram. eq.
!       else if (idivbflt .EQ. 3) then ! test
!       endif ! end of if idivbflt .......

* Obtain divB-free(?) field or energy correction.

!$omp parallel do private(i,k,j,aa,bb)
       do i = 0, ii + 1
         if (.NOT. ((ldistant) .AND. (i .EQ. 0))) then

           do k = 1, kk
           do j = 0, jj ! 1, jj - 1
             aa = (br2(j,k,i)**2+bt2(j,k,i)**2+bp2(j,k,i)**2)*0.5E+00 ! mag ene.
             if ((j .EQ. 0) .OR. (j .EQ. jj)) then ! rotational axis
               br2(j,k,i) = br2(j,k,i)
     &                  -((divbpsi2(j,k  ,i  )+divbpsi2(j-1,k  ,i  )
     &                    +divbpsi2(j,k-1,i  )+divbpsi2(j-1,k-1,i  ))
     &                   -(divbpsi2(j,k  ,i-1)+divbpsi2(j-1,k  ,i-1)
     &                    +divbpsi2(j,k-1,i-1)+divbpsi2(j-1,k-1,i-1)))
     &                  / (rrb(i) - rrb(i-1)) * 0.25E+00 * dtclean
             else
               br2(j,k,i) = br2(j,k,i)
     &                  - ((divbpsi2(j,k  ,i  )+divbpsi2(j-1,k  ,i  )
     &                     +divbpsi2(j,k-1,i  )+divbpsi2(j-1,k-1,i  ))
     &                    -(divbpsi2(j,k  ,i-1)+divbpsi2(j-1,k  ,i-1)
     &                     +divbpsi2(j,k-1,i-1)+divbpsi2(j-1,k-1,i-1)))
     &                   / (rrb(i) - rrb(i-1)) * 0.25E+00 * dtclean
               bt2(j,k,i) = bt2(j,k,i)
     &                  -(((divbpsi2(j  ,k,i  )+divbpsi2(j  ,k-1,i  ))
     &                    -(divbpsi2(j-1,k,i  )+divbpsi2(j-1,k-1,i  )))
     &                    / rrb(i  )
     &                  + ((divbpsi2(j  ,k,i-1)+divbpsi2(j  ,k-1,i-1))
     &                    -(divbpsi2(j-1,k,i-1)+divbpsi2(j-1,k-1,i-1)))
     &                    / rrb(i-1))
     &                   / (4.0E+00 * dth) * dtclean
               bp2(j,k,i) = bp2(j,k,i)
     &                  - ((divbpsi2(j  ,k  ,i  )
     &                     -divbpsi2(j  ,k-1,i  ))/rrb(i  )/sainb(j  )
     &                  +  (divbpsi2(j-1,k  ,i  )
     &                     -divbpsi2(j-1,k-1,i  ))/rrb(i  )/sainb(j-1)
     &                  +  (divbpsi2(j  ,k  ,i-1)
     &                     -divbpsi2(j  ,k-1,i-1))/rrb(i-1)/sainb(j  )
     &                  +  (divbpsi2(j-1,k  ,i-1)
     &                     -divbpsi2(j-1,k-1,i-1))/rrb(i-1)/sainb(j-1))
     &                   / (4.0E+00 * dph) * dtclean
             endif
             if (ladjstpg) then
               bb=(br2(j,k,i)**2+bt2(j,k,i)**2+bp2(j,k,i)**2)*0.5E+00 ! Pg adjust, usually off.
               pg2(j,k,i) = pg2(j,k,i) - (bb - aa) * (gammax-1.0E+00)
             endif
           enddo
           enddo

* overlap along longitudinal direction
           do j = 0, jj
             br2(j,   0,i) = br2(j,kk,i)
             br2(j,kk+1,i) = br2(j, 1,i)
             bt2(j,   0,i) = bt2(j,kk,i)
             bt2(j,kk+1,i) = bt2(j, 1,i)
             bp2(j,   0,i) = bp2(j,kk,i)
             bp2(j,kk+1,i) = bp2(j, 1,i)
           enddo
           do j = 0, jj
             br2(j,  -1,i) = br2(j,kk-1,i)
             br2(j,kk+2,i) = br2(j,   2,i)
             bt2(j,  -1,i) = bt2(j,kk-1,i)
             bt2(j,kk+2,i) = bt2(j,   2,i)
             bp2(j,  -1,i) = bp2(j,kk-1,i)
             bp2(j,kk+2,i) = bp2(j,   2,i)
           enddo

         endif
       enddo
!$omp end parallel do

       return
       endsubroutine

*
*--------------------------------------------------------------
*  1) Set periodicity as for the longitudinal direction
*  2) Confirm the theta- and phi-component of vector properties.
*        to be null on the rotational axis.
*  3) Avarage the scalar properties
*        and radial components of vector properties.
*--------------------------------------------------------------
*
       subroutine axibnd(ii,jj,kk,
     &                   gammax,kstep,omega,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                   rr,sain,kosa,sinph,cosph)
       implicit none
* arguments
       intent(in)    ::  ii,jj,kk
       intent(in)    ::  gammax,kstep,omega
       intent(inout) ::  ro2,pg2,ur2,ut2,up2,br2,bt2,bp2
       intent(in)    ::  rr,sain,kosa,sinph,cosph
       real    gammax, omega
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64) ! number of grid for each three direction.
       integer kstep(0:jj)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    rr(-2:ii+3)
       real    sain(0:jj), kosa(0:jj)
       real    sinph(0:kk+1),  cosph(0:kk+1)
* ldistant
       logical ldistant
       parameter(ldistant = .false.)
* local
       integer i, k, l
       real    wsum(16)
       real    cosphi, sinphi
       real    costh1, sinth1, costh2, sinth2
       real    xxxx1, xxyy1, xxzz1, xxrr1, xxtt1, xxpp1
       real    xxxx2, xxyy2, xxzz2, xxrr2, xxtt2, xxpp2
*       logical lchoice
       integer j, k2, jend, kstep2
       integer jendmax
*       parameter(jendmax = jj / 2)
       real    wsum2(16,1:kk) !  wsum1(16,1:kk)
       real    avey ! , avex, avexy, devxx
       real    aa, bb, cc ! , ff, gg ! , cc, dd, ee
* averaging variables at pole
       logical lpoleave
       parameter(lpoleave = .false.)
* averaging variables AROUND pole axis. In 2D, of no meanings.
       logical laxiave
       parameter(laxiave = .true.)
       logical lavecart
       parameter(lavecart = .false.) ! if false, averaging in spherical.
* whether B_phi will NOT be smoothed, even when laxiave = ON
       logical lnoavebp
       parameter(lnoavebp = .true.) ! there is another line to define this parameter.
* criteria of small Mr, to determine to enforce zero
       real     mrzero
       parameter(mrzero = 1.0E-00) ! mind mrzero = 1 / rr**2 or same order.
* criteria of|Br/B|
       real    brabsb
       parameter(brabsb = 1.0E-01)

      jendmax = jj / 2
* supress warnings
       aa = gammax

* averaging around solar rotational axis : the parallel comp. take too much time, here
       if (kk .LE. 3) then ! 2d-test case.

!$omp parallel do private(i,k)
           do i = -2, ii + 3 ! 0, ii + 1
           do k = 1, kk
**
             if (.NOT. lpoleave) then
               ro2( 0,k,i) = ro2(   1,k,i) !  not always needed..
               ro2(jj,k,i) = ro2(jj-1,k,i)
               pg2( 0,k,i) = pg2(   1,k,i)
               pg2(jj,k,i) = pg2(jj-1,k,i)
               ur2( 0,k,i) = ur2(   1,k,i)
               ur2(jj,k,i) = ur2(jj-1,k,i)
             endif
**
             ut2( 0,k,i) = 0.0E+00
             ut2(jj,k,i) = 0.0E+00
             up2( 0,k,i) = 0.0E+00
             up2(jj,k,i) = 0.0E+00
             bt2( 0,k,i) = 0.0E+00
             bt2(jj,k,i) = 0.0E+00
             bp2( 0,k,i) = 0.0E+00
             bp2(jj,k,i) = 0.0E+00
           enddo
           enddo
!$omp end parallel do

         jend = 0

       else

         if (laxiave) then

           jend = 1
           do j = 1, jj / 2
             if (kstep(j) .GT. 1) jend = j
           enddo
           if (jend .GT. jendmax) jend = jendmax

!$omp parallel do default(private)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(sain,kosa,sinph,cosph)
!$omp&            shared(jend,kstep,gammax)
!$omp&            shared(ii,jj,kk)
           do i = -2, ii + 3 ! 0, ii + 1
             do j = 1, jend
               kstep2 = kstep(j)
               do k = 1, kk
                 wsum2( 1,k) = ro2(   j,k,i)
                 wsum2( 2,k) = ro2(jj-j,k,i)
                 cosphi = cosph(k)
                 sinphi = sinph(k)
                 costh1 = kosa(j)
                 sinth1 = sain(j)
                 costh2 =-costh1 ! south
                 sinth2 = sinth1
**
*                 costh1 = 1.0E+00 ! Vr will be isolated...
*                 sinth1 = 0.0E+00
*                 costh2 =-1.0E+00
*                 sinth2 = 0.0E+00
**
!                 if (iavepg .EQ. 1) then ! energy
!                 else if (iavepg .EQ. 2) then
!                 else
                   wsum2( 3,k) = pg2(   j,k,i)
                   wsum2( 4,k) = pg2(jj-j,k,i)
!                 endif

!                 if (lavecart) then ! cart. or spherical
!                 else
                   wsum2( 5,k) = ut2(   j,k,i) * ro2(   j,k,i)
                   wsum2( 6,k) = ut2(jj-j,k,i) * ro2(jj-j,k,i)
                   wsum2( 7,k) = up2(   j,k,i) * ro2(   j,k,i)
                   wsum2( 8,k) = up2(jj-j,k,i) * ro2(jj-j,k,i)
                   wsum2( 9,k) = ur2(   j,k,i) * ro2(   j,k,i)
                   wsum2(10,k) = ur2(jj-j,k,i) * ro2(jj-j,k,i)
!                 endif

                 wsum2(11,k) = bt2(   j,k,i)
                 wsum2(12,k) = bt2(jj-j,k,i)
                 wsum2(13,k) = bp2(   j,k,i)
                 wsum2(14,k) = bp2(jj-j,k,i)
                 wsum2(15,k) = br2(   j,k,i)
                 wsum2(16,k) = br2(jj-j,k,i)
               enddo ! end of k-loop

               if (kstep2 .GT. 1) then
                 do l = 1, 16
* averaging
                   do k = 1, kk, kstep2
                     avey = 0.0E+00
                     do k2 = 0, kstep2 - 1
                       avey = avey + wsum2(l,k+k2)
                     enddo
                     avey = avey / float(kstep2)
                     do k2 = 0, kstep2 - 1
                       wsum2(l,k+k2) = avey
                     enddo ! k-local-loop
                   enddo !   k-loop
                 enddo !     l-loop
               endif

               do k = 1, kk
                 cosphi = cosph(k)
                 sinphi = sinph(k)
                 costh1 = kosa(j)
                 sinth1 = sain(j)
                 costh2 =-costh1 ! south
                 sinth2 = sinth1
**
*                 costh1 = 1.0E+00 ! Vr will be isolated...
*                 sinth1 = 0.0E+00
*                 costh2 =-1.0E+00
*                 sinth2 = 0.0E+00
**
                 ro2(   j,k,i) =  wsum2( 1,k)
                 ro2(jj-j,k,i) =  wsum2( 2,k)

!                 if (lavecart) then ! Cart. or spherical
!                 else
                   ut2(   j,k,i) = wsum2( 5,k) / wsum2(1,k)
                   ut2(jj-j,k,i) = wsum2( 6,k) / wsum2(2,k)
                   up2(   j,k,i) = wsum2( 7,k) / wsum2(1,k)
                   up2(jj-j,k,i) = wsum2( 8,k) / wsum2(2,k)
                   ur2(   j,k,i) = wsum2( 9,k) / wsum2(1,k)
                   ur2(jj-j,k,i) = wsum2(10,k) / wsum2(2,k)
!                 endif

                 bt2(   j,k,i) = wsum2(11,k)
                 bt2(jj-j,k,i) = wsum2(12,k)
!                 if (.NOT. lnoavebp) then
!                   bp2(   j,k,i) = wsum2(13,k)
!                   bp2(jj-j,k,i) = wsum2(14,k)
!                 endif
                 br2(   j,k,i) =  wsum2(15,k)
                 br2(jj-j,k,i) =  wsum2(16,k)
!                 if (iavepg .EQ. 1) then
!                 else if (iavepg .EQ. 2) then
!                 else
                   pg2(   j,k,i) = wsum2( 3,k)
                   pg2(jj-j,k,i) = wsum2( 4,k)
!                  endif
               enddo ! end  k-loop
             enddo ! end j-loop
           enddo ! end i-loop
!$omp end parallel do
         else
           jend = 0
         endif ! end if (laxiave)

!$omp parallel do default(private)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(sain,kosa,sinph,cosph)
!$omp&            shared(ii,jj,kk)
         do i = -2, ii + 3 ! 0, ii + 1
           cc = 1.0E+00 / 3.0E+00
           do l = 1, 16
             wsum(l) = 0.0E+00
           enddo

             do k = 1, kk
               cosphi = cosph(k)
               sinphi = sinph(k)
               costh1 = 1.0E+00 !  at Pole
               sinth1 = 0.0E+00
               costh2 =-1.0E+00
               sinth2 = 0.0E+00

!               if (lpoleave) then
!               else

                 wsum( 1) = wsum( 1) + ro2(   1,k,i)
                 wsum( 2) = wsum( 2) + ro2(jj-1,k,i)
!                 if (iavepg .EQ. 1) then
!                 else if (iavepg .EQ. 2) then
!                 else
                   wsum( 3) = wsum( 3) + pg2(   1,k,i)
                   wsum( 4) = wsum( 4) + pg2(jj-1,k,i)
!                 endif

                 xxxx1 = ur2(   1,k,i) * sinth1 * cosphi
     &                 + ut2(   1,k,i) * costh1 * cosphi
     &                 - up2(   1,k,i)          * sinphi
                 xxyy1 = ur2(   1,k,i) * sinth1 * sinphi
     &                 + ut2(   1,k,i) * costh1 * sinphi
     &                 + up2(   1,k,i)          * cosphi
                 xxzz1 = ur2(   1,k,i) * costh1
     &                 - ut2(   1,k,i) * sinth1
                 wsum( 5) = wsum( 5) + xxxx1 * ro2(   1,k,i)
                 wsum( 7) = wsum( 7) + xxyy1 * ro2(   1,k,i)
                 wsum( 9) = wsum( 9) + xxzz1 * ro2(   1,k,i)
                 xxxx2 = ur2(jj-1,k,i) * sinth2 * cosphi
     &                 + ut2(jj-1,k,i) * costh2 * cosphi
     &                 - up2(jj-1,k,i)          * sinphi
                 xxyy2 = ur2(jj-1,k,i) * sinth2 * sinphi
     &                 + ut2(jj-1,k,i) * costh2 * sinphi
     &                 + up2(jj-1,k,i)          * cosphi
                 xxzz2 = ur2(jj-1,k,i) * costh2
     &                 - ut2(jj-1,k,i) * sinth2
                 wsum( 6) = wsum( 6) + xxxx2 * ro2(jj-1,k,i)
                 wsum( 8) = wsum( 8) + xxyy2 * ro2(jj-1,k,i)
                 wsum(10) = wsum(10) + xxzz2 * ro2(jj-1,k,i)

                 xxxx1 = br2(   0,k,i) * sinth1 * cosphi
     &                 + bt2(   1,k,i) * costh1 * cosphi
     &                 - bp2(   1,k,i)          * sinphi
                 xxyy1 = br2(   0,k,i) * sinth1 * sinphi
     &                 + bt2(   1,k,i) * costh1 * sinphi
     &                 + bp2(   1,k,i)          * cosphi
                 xxzz1 = br2(   0,k,i) * costh1
     &                 - bt2(   1,k,i) * sinth1
                 wsum(11) = wsum(11) + xxxx1
                 wsum(13) = wsum(13) + xxyy1
                 wsum(15) = wsum(15) + xxzz1
                 xxxx2 = br2(jj  ,k,i) * sinth2 * cosphi
     &                 + bt2(jj-1,k,i) * costh2 * cosphi
     &                 - bp2(jj-1,k,i)          * sinphi
                 xxyy2 = br2(jj  ,k,i) * sinth2 * sinphi
     &                 + bt2(jj-1,k,i) * costh2 * sinphi
     &                 + bp2(jj-1,k,i)          * cosphi
                 xxzz2 = br2(jj  ,k,i) * costh2
     &                 - bt2(jj-1,k,i) * sinth2
                 wsum(12) = wsum(12) + xxxx2
                 wsum(14) = wsum(14) + xxyy2
                 wsum(16) = wsum(16) + xxzz2

!               endif ! end-if (lpoleave)

             enddo ! end k-loop


           do l = 1, 16
              wsum(l) = wsum(l) / float(kk)
           enddo
           do k = 1, kk
             cosphi = cosph(k)
             sinphi = sinph(k)
             costh1 = 1.0E+00 !  at Pole
             sinth1 = 0.0E+00
             costh2 =-1.0E+00
             sinth2 = 0.0E+00
             aa = wsum( 1)
             bb = wsum( 2)
             ro2( 0,k,i) = aa
             ro2(jj,k,i) = bb

!             if (iavepg .EQ. 1) then
!             else if (iavepg .EQ. 2) then
!             else
               pg2( 0,k,i) = wsum( 3)
               pg2(jj,k,i) = wsum( 4)
!             endif

             xxxx1 = wsum( 5) / wsum( 1)
             xxyy1 = wsum( 7) / wsum( 1)
             xxzz1 = wsum( 9) / wsum( 1)
             xxrr1 = xxxx1 * sinth1 * cosphi
     &             + xxyy1 * sinth1 * sinphi
     &             + xxzz1 * costh1
             xxtt1 = xxxx1 * costh1 * cosphi
     &             + xxyy1 * costh1 * sinphi
     &             - xxzz1 * sinth1
             xxpp1 =-xxxx1          * sinphi
     &             + xxyy1          * cosphi
             ur2(   0,k,i) = xxrr1
             ut2(   0,k,i) = xxtt1
             up2(   0,k,i) = xxpp1
             xxxx2 = wsum( 6) / wsum( 2)
             xxyy2 = wsum( 8) / wsum( 2)
             xxzz2 = wsum(10) / wsum( 2)
             xxrr2 = xxxx2 * sinth2 * cosphi
     &             + xxyy2 * sinth2 * sinphi
     &             + xxzz2 * costh2
             xxtt2 = xxxx2 * costh2 * cosphi
     &             + xxyy2 * costh2 * sinphi
     &             - xxzz2 * sinth2
             xxpp2 =-xxxx2          * sinphi
     &             + xxyy2          * cosphi
             ur2(jj  ,k,i) = xxrr2
             ut2(jj  ,k,i) = xxtt2
             up2(jj  ,k,i) = xxpp2

             xxxx1 = wsum(11)
             xxyy1 = wsum(13)
             xxzz1 = wsum(15)
             xxrr1 = xxxx1 * sinth1 * cosphi
     &             + xxyy1 * sinth1 * sinphi
     &             + xxzz1 * costh1
             xxtt1 = xxxx1 * costh1 * cosphi
     &             + xxyy1 * costh1 * sinphi
     &             - xxzz1 * sinth1
             xxpp1 =-xxxx1          * sinphi
     &             + xxyy1          * cosphi
             br2(   0,k,i) = xxrr1
             bt2(   0,k,i) = xxtt1
             bp2(   0,k,i) = xxpp1
             xxxx2 = wsum(12)
             xxyy2 = wsum(14)
             xxzz2 = wsum(16)
             xxrr2 = xxxx2 * sinth2 * cosphi
     &             + xxyy2 * sinth2 * sinphi
     &             + xxzz2 * costh2
             xxtt2 = xxxx2 * costh2 * cosphi
     &             + xxyy2 * costh2 * sinphi
     &             - xxzz2 * sinth2
             xxpp2 =-xxxx2          * sinphi
     &             + xxyy2          * cosphi
             br2(jj  ,k,i) = xxrr2
             bt2(jj  ,k,i) = xxtt2
             bp2(jj  ,k,i) = xxpp2
           enddo

         enddo ! end i-loop
!$omp end parallel do
       endif

* for v//B at i = 0
!$omp parallel do default(private)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(jend)
!$omp&            shared(omega,rr,sain)
!$omp&            shared(ii,jj,kk)
       do i = -2, ii + 3
         if (i .EQ. 0) then
           do k = 1, kk
           do j = 0, jj
!             if (ldistant) then
!             else
               if ((j .LE. jend) .OR. (j .GE. jj-jend)) then
                 bb = br2(j,k,i)**2 + bt2(j,k,i)**2 + bp2(j,k,i)**2
                 aa = sqrt(bb)
                 if (aa .GT. 1.0E-05) then
                   bb = abs(br2(j,k,i)) / aa
                   if ((ro2(j,k,i)*ur2(j,k,i) .LT. mrzero) .OR.
     &                 (bb .LT. brabsb)) then ! freeze
                     ur2(j,k,i) = 0.0E+00
                     ut2(j,k,i) = 0.0E+00
                     up2(j,k,i) = 0.0E+00
                   else
                     ut2(j,k,i) = ur2(j,k,i)*bt2(j,k,i)/br2(j,k,i)
                     up2(j,k,i) = ur2(j,k,i)*bp2(j,k,i)/br2(j,k,i)
                   endif
                 endif
               endif
!             endif
           enddo ! end j-loop
           enddo ! end k-loop
         endif ! end if (i .EQ. 0)
       enddo ! end dummy-i-loop
!$omp end parallel do

!$omp parallel do private(i,k,k2)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(ii,jj,kk)
       do i = -2, ii + 3 ! i = 0, ii + 1
         do k = 1, kk
           k2 = mod(k + kk/2,kk)
           if (k2 .EQ. 0) k2 = kk
           ro2(  -2,k,i) = ro2(   2,k2,i)
           pg2(  -2,k,i) = pg2(   2,k2,i)
           ur2(  -2,k,i) = ur2(   2,k2,i)
           ut2(  -2,k,i) =-ut2(   2,k2,i)
           up2(  -2,k,i) =-up2(   2,k2,i)
           br2(  -2,k,i) = br2(   2,k2,i)
           bt2(  -2,k,i) =-bt2(   2,k2,i)
           bp2(  -2,k,i) =-bp2(   2,k2,i)
           ro2(  -1,k,i) = ro2(   1,k2,i)
           pg2(  -1,k,i) = pg2(   1,k2,i)
           ur2(  -1,k,i) = ur2(   1,k2,i)
           ut2(  -1,k,i) =-ut2(   1,k2,i)
           up2(  -1,k,i) =-up2(   1,k2,i)
           br2(  -1,k,i) = br2(   1,k2,i)
           bt2(  -1,k,i) =-bt2(   1,k2,i)
           bp2(  -1,k,i) =-bp2(   1,k2,i)
           ro2(jj+1,k,i) = ro2(jj-1,k2,i)
           pg2(jj+1,k,i) = pg2(jj-1,k2,i)
           ur2(jj+1,k,i) = ur2(jj-1,k2,i)
           ut2(jj+1,k,i) =-ut2(jj-1,k2,i)
           up2(jj+1,k,i) =-up2(jj-1,k2,i)
           br2(jj+1,k,i) = br2(jj-1,k2,i)
           bt2(jj+1,k,i) =-bt2(jj-1,k2,i)
           bp2(jj+1,k,i) =-bp2(jj-1,k2,i)
           ro2(jj+2,k,i) = ro2(jj-2,k2,i)
           pg2(jj+2,k,i) = pg2(jj-2,k2,i)
           ur2(jj+2,k,i) = ur2(jj-2,k2,i)
           ut2(jj+2,k,i) =-ut2(jj-2,k2,i)
           up2(jj+2,k,i) =-up2(jj-2,k2,i)
           br2(jj+2,k,i) = br2(jj-2,k2,i)
           bt2(jj+2,k,i) =-bt2(jj-2,k2,i)
           bp2(jj+2,k,i) =-bp2(jj-2,k2,i)
         enddo
       enddo
!$omp end parallel do

       return
       endsubroutine

*
*---------------------------------------------------------------------
*
       subroutine setkstep(ii,jj,kk,kstep,rr,sain) ! assuming uniform grid along both lat and lon.
       implicit none
* arguments
       intent(in)  ::      ii,jj,kk
       intent(out) ::      kstep
       intent(in)  ::            rr,sain
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       integer kstep(0:jj)
       real    rr(-2:ii+3),sain(0:jj)
* local
       integer j, i, k
       real    aa, bb, cc, dd
* const.
       real    pi
       parameter(pi = 3.1415926535897932385E+00)
       logical laxiave
       parameter(laxiave = .true.)
       logical ldistant
       parameter(ldistant = .false.)

       kstep( 0) = 1
       kstep(jj) = 1
       if ((kk .GT. 1) .AND. (laxiave)) then
         do j = 1, jj - 1
           aa = 1.0E+10
           do i = 0, ii
             bb = rr(i+1) - rr(i)
             cc = rr(i) * pi / float(jj)
             dd = min(bb,cc)
     &          /(rr(i) * sain(j) * (2.0E+00 * pi / float(kk)))

             if (dd .LT. aa) aa = dd
           enddo
           bb = 1.0E+00
           i = -1
           kstep(j) = 0
           do while ((bb .GT. 0.0E+00) .AND. (i .LE. 10))
             i = i + 1
             k = 2**i
             bb = aa - float(k)
           enddo
           if (k .GT. kk / 8) k = kk / 8 ! kk > 8 is assumed.
           kstep(j) = k
         enddo
       else
         do j = 1, jj - 1
           kstep(j) = 1
         enddo
       endif
       write(*,'(A,$)')  'ksteps are.. '
       do j = 0, jj
         write(*,'(i4,''('',i3,'')'',$)') j, kstep(j)
       enddo
       write(*,'(A)')  ' '

       return
       endsubroutine

*
*---------------------------------------------------------------------
*
       subroutine lngbnd(ii,jj,kk,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
       implicit none
* local
       intent(in)    ::  ii,jj,kk
       intent(inout) ::  ro2,pg2,ur2,ut2,up2,br2,bt2,bp2
       integer i, j
* number of grid for each three direction.
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)

* periodicity along longitudinal direction
!$omp parallel do private(i,j)
       do i = -2, ii + 3
         do j =-2, jj + 2
           ro2(j,0,i) = ro2(j,kk,i)
           pg2(j,0,i) = pg2(j,kk,i)
           ur2(j,0,i) = ur2(j,kk,i)
           ut2(j,0,i) = ut2(j,kk,i)
           up2(j,0,i) = up2(j,kk,i)
           br2(j,0,i) = br2(j,kk,i)
           bt2(j,0,i) = bt2(j,kk,i)
           bp2(j,0,i) = bp2(j,kk,i)
         enddo
         do j =-2, jj + 2
           ro2(j,kk+1,i) = ro2(j,1,i)
           pg2(j,kk+1,i) = pg2(j,1,i)
           ur2(j,kk+1,i) = ur2(j,1,i)
           ut2(j,kk+1,i) = ut2(j,1,i)
           up2(j,kk+1,i) = up2(j,1,i)
           br2(j,kk+1,i) = br2(j,1,i)
           bt2(j,kk+1,i) = bt2(j,1,i)
           bp2(j,kk+1,i) = bp2(j,1,i)
         enddo
         do j =-2, jj + 2
           ro2(j,-1,i) = ro2(j,kk-1,i)
           pg2(j,-1,i) = pg2(j,kk-1,i)
           ur2(j,-1,i) = ur2(j,kk-1,i)
           ut2(j,-1,i) = ut2(j,kk-1,i)
           up2(j,-1,i) = up2(j,kk-1,i)
           br2(j,-1,i) = br2(j,kk-1,i)
           bt2(j,-1,i) = bt2(j,kk-1,i)
           bp2(j,-1,i) = bp2(j,kk-1,i)
         enddo
         do j =-2, jj + 2
           ro2(j,kk+2,i) = ro2(j,2,i)
           pg2(j,kk+2,i) = pg2(j,2,i)
           ur2(j,kk+2,i) = ur2(j,2,i)
           ut2(j,kk+2,i) = ut2(j,2,i)
           up2(j,kk+2,i) = up2(j,2,i)
           br2(j,kk+2,i) = br2(j,2,i)
           bt2(j,kk+2,i) = bt2(j,2,i)
           bp2(j,kk+2,i) = bp2(j,2,i)
         enddo
       enddo
!$omp end parallel do

       return
       endsubroutine


*--------------------------------------------------------------------
*
       subroutine prc2nd(ii,jj,kk,gammax,dt,kstep,ncal,
     &                       pg3,ur3,ut3,up3,
     &                       pg9,ur9,ut9,up9,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                   rr,theta,phi,sain,     sainb,kosab)
       implicit none
* arguments
       intent(in)    ::  ii,jj,kk,gammax,dt,kstep,ncal
       intent(out)   ::      pg3,ur3,ut3,up3
       intent(out)   ::      pg9,ur9,ut9,up9
       intent(inout) ::      pg2,ur2,ut2,up2
       intent(in)    ::  ro2,                br2,bt2,bp2
       intent(in)    ::  rr,theta,phi,sain,     sainb,kosab
* interface
       integer ii, jj, kk
       real    gammax,dt
       integer kstep(0:jj)
       integer ncal
* large array
       real    rr(-2:ii+3),theta(0:jj),phi(-1:kk+2)
       real    sain(0:jj)
       real    sainb(0:jj-1), kosab(0:jj-1)
       real    pg3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ut3(-2:jj+2,-1:kk+2,-2:ii+3),up3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    pg9(-2:jj+2,-1:kk+2,-2:ii+9)
       real    ur9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ut9(-2:jj+2,-1:kk+2,-2:ii+3),up9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
* local
       integer kstepb(0:jj)
       integer i, j, k, n, ntime
       real    aa, bb, cc, dd, ee
       integer k2, kstep2
       integer il1, il2, il3
       integer jl1, jl2, jl3
       integer kl1, kl2, kl3
       real    dra0, dra1, dra2
       real    dth0, dth1
       real    dph0, dph2
       real    rrb1, rrb2, rrbc
       real    sinb1, sinb2, cosb1, cosb2, sinbc, cotbc
       real    xxur, xxut, xxup, xxpg
       real    dtlocal
       real    tcndr1, tcndt1, tcndp1
       real    tcndr2, tcndt2, tcndp2
       real    tt1r1, tt1r2, tt1t1, tt1t2, tt1p1, tt1p2
       real    ro1r1, ro1r2, ro1t1, ro1t2, ro1p1, ro1p2
       real    pg1r1, pg1r2, pg1t1, pg1t2, pg1p1, pg1p2
       real    ur1r1, ut1r1, up1r1, ur1r2, ut1r2, up1r2
       real    ur1t1, ut1t1, up1t1, ur1t2, ut1t2, up1t2
       real    ur1p1, ut1p1, up1p1, ur1p2, ut1p2, up1p2
       real    br1r1, bt1r1, bp1r1, br1r2, bt1r2, bp1r2
       real    br1t1, bt1t1, bp1t1, br1t2, bt1t2, bp1t2
       real    br1p1, bt1p1, bp1p1, br1p2, bt1p2, bp1p2
* heat conduction
       logical ltcndct
       parameter(ltcndct = .true.)
       integer itcndtype
       parameter(itcndtype = 0) ! 1=orthodox(a), 2=orthodox(b), otherwise=modified,fast,informal !
       logical  ltcmgpar
       parameter(ltcmgpar = .true.) ! if true, heat may preferentially move along B
       real    cbeta
       parameter(cbeta = 1.0E-02) ! criteria 1/beta old & normal, controlling parallelness
       real    tcndfact
       parameter(tcndfact = 1.00E+00) ! weak for ich6daily
*       parameter(tcndfact = 6.28E+00) ! normalized Spitzer value; 5.0e-07 * tmp0^(7/2) t_0 / l0^2 / p0
*       parameter(tcndfact = 2.00E+01) ! enhanced, for test (factor may ranges 1--10e-07)
*       parameter(tcndfact = 1.00E-02) ! small, for stabilizing, not necessary.
* fluid viscosity,,,, d(v rho d(v)) etc...
       logical lflwvisc
       parameter(lflwvisc = .true.) ! may be false when distant
       real    viscfact
       parameter(viscfact = 5.0E-03)
* constant(s)
       real    pi
       parameter(pi = 3.1415926535897932385E+00)

* show setting
       if (ncal .EQ. 1) then
         write(*,*) ' ltcndct  = ', ltcndct
         write(*,*) ' ltcmgpar = ', ltcmgpar
         write(*,*) ' tcndfact = ', tcndfact
         write(*,*) ' cbeta    = ', cbeta
         write(*,*) ' lflwvisc = ', lflwvisc
         write(*,*) ' viscfact = ', viscfact
         write(19,*) ' ltcndct  = ', ltcndct
         write(19,*) ' ltcmgpar = ', ltcmgpar
         write(19,*) ' tcndfact = ', tcndfact
         write(19,*) ' cbeta    = ', cbeta
         write(19,*) ' lflwvisc = ', lflwvisc
         write(19,*) ' viscfact = ', viscfact
       endif

* copy for sake of safety
       do j = 0, jj
         kstepb(j) = kstep(j)
       enddo

* heat conduction,. ! -----------------------------------------
       if (ltcndct) then

!$omp parallel do private(i,k,j)
!$omp&            shared(pg2,pg3)
!$omp&            shared(ii,jj,kk)
         do i = 0, ii + 1
           do k = -1, kk + 2
           do j = -1, jj + 1
             pg3(j,k,i) = pg2(j,k,i)
           enddo
           enddo
         enddo
!$omp end parallel do

           dtlocal = 1.0E+30
!$omp parallel do private(i,k,j,aa,bb,cc,dd)
!$omp&            shared(gammax)
!$omp&            shared(ro2,pg2)
!$omp&            shared(rr,theta,phi,sain,kstepb)
!$omp&            shared(ii,jj,kk)
!$omp&            reduction(min:dtlocal)
           do i = 1, ii
           do k = 1, kk
           do j = 1, jj - 1
             aa = (rr(i+1) - rr(i-1)) * 0.5E+00
             bb = rr(i) * (theta(j+1) - theta(j-1)) * 0.5E+00
             cc = rr(i) * sain(j) * (phi(k+1) - phi(k-1)) * 0.5E+00
     &          * float(kstepb(j))
             dd = min(aa,bb,cc)
             dd = 0.5E+00 * dd**2
     &          / (tcndfact * (gammax - 1.0E+00)
     &                      * (pg2(j,k,i)/ro2(j,k,i))**2.5E+00)
     &          * 0.80E+00 ! CFL num.
             if (itcndtype .EQ. 1) dd = dd * ro2(j,k,i)
             if (dd .LT. dtlocal) dtlocal = dd
           enddo
           enddo
           enddo
!$omp end parallel do

           ntime = int(dt / dtlocal) + 1
           if (ncal .EQ. 1) then
             write(*,*) 'Heat Condution --- '
             write(*,*) 'Ntime dt_local(original) dt_global = '
             write(*,*)  ntime, dtlocal, dt
           endif
           dtlocal = dt / float(ntime)

           do n = 1, ntime

!$omp parallel do default(private)
!$omp&            shared(dtlocal,gammax)
!$omp&            shared(kstepb)
!$omp&            shared(ro2,br2,bt2,bp2)
!$omp&            shared(pg9,pg3)
!$omp&            shared(rr,theta,phi,sain,     sainb,kosab)
!$omp&            shared(ii,jj,kk)
             do i = 0, ii + 1
             do k = 1, kk
             do j = 1, jj - 1 ! MIND the range...

               il1 = i - 1
               il2 = i
               il3 = i + 1
               if (il1 .LT.    0) il1 =    0
               if (il3 .GT. ii+1) il3 = ii+1
               jl1 = j - 1
               jl2 = j
               jl3 = j + 1
               if (jl1 .LT.   -1) jl1 =   -1
               if (jl3 .GT. jj+1) jl3 = jj+1
               kl1 = k -  kstepb(j) + kk * 2
               kl1 = mod(kl1,kk)
               kl2 = k
               kl3 = k +  kstepb(j) + kk * 2
               kl3 = mod(kl3,kk)

               rrb1 =(rr(il3) + rr(i)) * 0.5E+00
               rrb2 =(rr(il1) + rr(i)) * 0.5E+00
               rrbc =(rrb1**2 + rrb1 * rrb2 + rrb2**2)
     &              /(rrb1 + rrb2) / 1.5E+00
               dra0 = rrb1 - rrb2
               dra2 =(rrb1**3 - rrb2**3) / 3.0E+00
               if (j .EQ. 0) then
                 sinb1 = sainb(j)
                 cosb1 = kosab(j)
                 sinb2 = 0.0E+00
                 cosb2 = 1.0E+00
                 dth0  = pi / float(jj) * 0.5E+00
                 sinbc = (cosb2 - cosb1) / dth0
               else if (j .EQ. jj) then
                 sinb1 = 0.0E+00
                 cosb1 =-1.0E+00
                 sinb2 = sainb(j-1)
                 cosb2 = kosab(j-1)
                 dth0  = pi / float(jj) * 0.5E+00
                 sinbc = (cosb2 - cosb1) / dth0
               else
                 sinb1 = sainb(j)
                 cosb1 = kosab(j)
                 sinb2 = sainb(j-1)
                 cosb2 = kosab(j-1)
                 dth0  = pi / float(jj)
                 sinbc = (cosb2 - cosb1) / dth0
               endif
               dth1 = cosb2 - cosb1

               dph0 = 2.0E+00 * pi / float(kk)
               dph2 = dph0

               dph0 = dph0 * float(kstepb(j))
               dph2 = dph2 * float(kstepb(j))

               if (ltcmgpar) then
                 br1r1 = (br2(j,k,il3) + br2(j,k,i)) * 0.5E+00
                 bt1r1 = (bt2(j,k,il3) + bt2(j,k,i)) * 0.5E+00
                 bp1r1 = (bp2(j,k,il3) + bp2(j,k,i)) * 0.5E+00
                 br1r2 = (br2(j,k,il1) + br2(j,k,i)) * 0.5E+00
                 bt1r2 = (bt2(j,k,il1) + bt2(j,k,i)) * 0.5E+00
                 bp1r2 = (bp2(j,k,il1) + bp2(j,k,i)) * 0.5E+00
                 br1t1 = (br2(jl3,k,i) + br2(j,k,i)) * 0.5E+00
                 bt1t1 = (bt2(jl3,k,i) + bt2(j,k,i)) * 0.5E+00
                 bp1t1 = (bp2(jl3,k,i) + bp2(j,k,i)) * 0.5E+00
                 br1t2 = (br2(jl1,k,i) + br2(j,k,i)) * 0.5E+00
                 bt1t2 = (bt2(jl1,k,i) + bt2(j,k,i)) * 0.5E+00
                 bp1t2 = (bp2(jl1,k,i) + bp2(j,k,i)) * 0.5E+00
                 br1p1 = (br2(j,kl3,i) + br2(j,k,i)) * 0.5E+00
                 bt1p1 = (bt2(j,kl3,i) + bt2(j,k,i)) * 0.5E+00
                 bp1p1 = (bp2(j,kl3,i) + bp2(j,k,i)) * 0.5E+00
                 br1p2 = (br2(j,kl1,i) + br2(j,k,i)) * 0.5E+00
                 bt1p2 = (bt2(j,kl1,i) + bt2(j,k,i)) * 0.5E+00
                 bp1p2 = (bp2(j,kl1,i) + bp2(j,k,i)) * 0.5E+00
               else
                 br1r1 = 0.0E+00
                 br1r2 = 0.0E+00
                 bt1r1 = 0.0E+00
                 bt1r2 = 0.0E+00
                 bp1r1 = 0.0E+00
                 bp1r2 = 0.0E+00
                 br1t1 = 0.0E+00
                 br1t2 = 0.0E+00
                 bt1t1 = 0.0E+00
                 bt1t2 = 0.0E+00
                 bp1t1 = 0.0E+00
                 bp1t2 = 0.0E+00
                 br1p1 = 0.0E+00
                 br1p2 = 0.0E+00
                 bt1p1 = 0.0E+00
                 bt1p2 = 0.0E+00
                 bp1p1 = 0.0E+00
                 bp1p2 = 0.0E+00
               endif
               ro1r1 = (ro2(j,k,il3) + ro2(j,k,i)) * 0.5E+00
               ro1r2 = (ro2(j,k,il1) + ro2(j,k,i)) * 0.5E+00
               ro1t1 = (ro2(jl3,k,i) + ro2(j,k,i)) * 0.5E+00
               ro1t2 = (ro2(jl1,k,i) + ro2(j,k,i)) * 0.5E+00
               ro1p1 = (ro2(j,kl3,i) + ro2(j,k,i)) * 0.5E+00
               ro1p2 = (ro2(j,kl1,i) + ro2(j,k,i)) * 0.5E+00
               pg1r1 = (pg3(j,k,il3) + pg3(j,k,i)) * 0.5E+00
               pg1r2 = (pg3(j,k,il1) + pg3(j,k,i)) * 0.5E+00
               pg1t1 = (pg3(jl3,k,i) + pg3(j,k,i)) * 0.5E+00
               pg1t2 = (pg3(jl1,k,i) + pg3(j,k,i)) * 0.5E+00
               pg1p1 = (pg3(j,kl3,i) + pg3(j,k,i)) * 0.5E+00
               pg1p2 = (pg3(j,kl1,i) + pg3(j,k,i)) * 0.5E+00

               tt1r1 = (pg3(j,k,il3)/ro2(j,k,il3)
     &                 +pg3(j,k,i  )/ro2(j,k,i  ))*0.5E+00
               tt1r2 = (pg3(j,k,il1)/ro2(j,k,il1)
     &                 +pg3(j,k,i  )/ro2(j,k,i  ))*0.5E+00
               tt1t1 = (pg3(jl3,k,i)/ro2(jl3,k,i)
     &                 +pg3(j  ,k,i)/ro2(j  ,k,i))*0.5E+00
               tt1t2 = (pg3(jl1,k,i)/ro2(jl1,k,i)
     &                 +pg3(j  ,k,i)/ro2(j  ,k,i))*0.5E+00
               tt1p1 = (pg3(j,kl3,i)/ro2(j,kl3,i)
     &                 +pg3(j,k  ,i)/ro2(j,k  ,i))*0.5E+00
               tt1p2 = (pg3(j,kl1,i)/ro2(j,kl1,i)
     &                 +pg3(j,k  ,i)/ro2(j,k  ,i))*0.5E+00

               aa = (pg3(jl2,kl2,il3)/ro2(jl2,kl2,il3)
     &              -pg3(jl2,kl2,il2)/ro2(jl2,kl2,il2))
     &             / dra0 * (br1r1**2 + pg1r1*cbeta)
     &            +((pg3(jl3,kl2,il3)/ro2(jl3,kl2,il3)
     &              +pg3(jl3,kl2,il2)/ro2(jl3,kl2,il2)) * 0.5E+00
     &             -(pg3(jl1,kl2,il3)/ro2(jl1,kl2,il3)
     &              +pg3(jl1,kl2,il2)/ro2(jl1,kl2,il2)) * 0.5E+00)
     &             /(dth0 * rrb1 * 2.0E+00)         * bt1r1 * br1r1
     &            +((pg3(jl2,kl3,il3)/ro2(jl2,kl3,il3)
     &              +pg3(jl2,kl3,il2)/ro2(jl2,kl3,il2)) * 0.5E+00
     &             -(pg3(jl2,kl1,il3)/ro2(jl2,kl1,il3)
     &              +pg3(jl2,kl1,il2)/ro2(jl2,kl1,il2)) * 0.5E+00)
     &             /(rrb1 * sinbc * dph0 * 2.0E+00) * bp1r1 * br1r1
               tcndr1 = aa * tt1r1**2.5E+00
     &                /(br1r1**2 + bt1r1**2 + bp1r1**2 + pg1r1*cbeta)
               aa = (pg3(jl2,kl2,il2)/ro2(jl2,kl2,il2)
     &              -pg3(jl2,kl2,il1)/ro2(jl2,kl2,il1))
     &             / dra0 * (br1r2**2 + pg1r2 * cbeta)
     &            +((pg3(jl3,kl2,il1)/ro2(jl3,kl2,il1)
     &              +pg3(jl3,kl2,il2)/ro2(jl3,kl2,il2)) * 0.5E+00
     &             -(pg3(jl1,kl2,il1)/ro2(jl1,kl2,il1)
     &              +pg3(jl1,kl2,il2)/ro2(jl1,kl2,il2)) * 0.5E+00)
     &             /(dth0 * rrb2 * 2.0E+00)         * bt1r2 * br1r2
     &            +((pg3(jl2,kl3,il1)/ro2(jl2,kl3,il1)
     &              +pg3(jl2,kl3,il2)/ro2(jl2,kl3,il2)) * 0.5E+00
     &             -(pg3(jl2,kl1,il1)/ro2(jl2,kl1,il1)
     &              +pg3(jl2,kl1,il2)/ro2(jl2,kl1,il2)) * 0.5E+00)
     &             /(rrb2 * sinbc * dph0 * 2.0E+00) * bp1r2 * br1r2
               tcndr2 = aa * tt1r2**2.5E+00
     &                /(br1r2**2 + bt1r2**2 + bp1r2**2 + pg1r2*cbeta)

               aa = (pg3(jl3,kl2,il2)/ro2(jl3,kl2,il2)
     &              -pg3(jl2,kl2,il2)/ro2(jl2,kl2,il2))
     &             /(rrbc * dth0) * (bt1t1**2 + pg1t1 * cbeta)
     &            +((pg3(jl3,kl3,il2)/ro2(jl3,kl3,il2)
     &              +pg3(jl2,kl3,il2)/ro2(jl2,kl3,il2)) * 0.5E+00
     &             -(pg3(jl3,kl1,il2)/ro2(jl3,kl1,il2)
     &              +pg3(jl2,kl1,il2)/ro2(jl2,kl1,il2)) * 0.5E+00)
     &             /(2.0E+00 * dph0 * rrbc * sinb1) * bp1t1 * bt1t1
     &            +((pg3(jl3,kl2,il3)/ro2(jl3,kl2,il3)
     &              +pg3(jl2,kl2,il3)/ro2(jl2,kl2,il3)) * 0.5E+00
     &             -(pg3(jl3,kl2,il1)/ro2(jl3,kl2,il1)
     &              +pg3(jl2,kl2,il1)/ro2(jl2,kl2,il1)) * 0.5E+00)
     &             /(dra0 * 2.0E+00)                * br1t1 * bt1t1
               tcndt1 = aa * tt1t1**2.5E+00
     &                /(br1t1**2 + bt1t1**2 + bp1t1**2 + pg1t1*cbeta)
               aa = (pg3(jl2,kl2,il2)/ro2(jl2,kl2,il2)
     &              -pg3(jl1,kl2,il2)/ro2(jl1,kl2,il2))
     &             /(rrbc * dth0) * (bt1t2**2 + pg1t2 * cbeta)
     &            +((pg3(jl1,kl3,il2)/ro2(jl1,kl3,il2)
     &              +pg3(jl2,kl3,il2)/ro2(jl2,kl3,il2)) * 0.5E+00
     &             -(pg3(jl1,kl1,il2)/ro2(jl1,kl1,il2)
     &              +pg3(jl2,kl1,il2)/ro2(jl2,kl1,il2)) * 0.5E+00)
     &             /(2.0E+00 * dph0 * rrbc * sinb2) * bp1t2 * bt1t2
     &            +((pg3(jl1,kl2,il3)/ro2(jl1,kl2,il3)
     &              +pg3(jl2,kl2,il3)/ro2(jl2,kl2,il3)) * 0.5E+00
     &             -(pg3(jl1,kl2,il1)/ro2(jl1,kl2,il1)
     &              +pg3(jl2,kl2,il1)/ro2(jl2,kl2,il1)) * 0.5E+00)
     &             /(dra0 * 2.0E+00)                * br1t2 * bt1t2
               tcndt2 = aa * tt1t2**2.5E+00
     &                /(br1t2**2 + bt1t2**2 + bp1t2**2 + pg1t2*cbeta)

               aa = (pg3(jl2,kl3,il2)/ro2(jl2,kl3,il2)
     &              -pg3(jl2,kl2,il2)/ro2(jl2,kl2,il2))
     &             /(dph0 * rrbc * sinbc) * (bp1p1**2 + pg1p1*cbeta)
     &            +((pg3(jl2,kl3,il3)/ro2(jl2,kl3,il3)
     &              +pg3(jl2,kl2,il3)/ro2(jl2,kl2,il3)) * 0.5E+00
     &             -(pg3(jl2,kl3,il1)/ro2(jl2,kl3,il1)
     &              +pg3(jl2,kl2,il1)/ro2(jl2,kl2,il1)) * 0.5E+00)
     &             /(dra0 * 2.0E+00)        * br1p1 * bp1p1
     &            +((pg3(jl3,kl3,il2)/ro2(jl3,kl3,il2)
     &              +pg3(jl3,kl2,il2)/ro2(jl3,kl2,il2)) * 0.5E+00
     &             -(pg3(jl1,kl3,il2)/ro2(jl1,kl3,il2)
     &              +pg3(jl1,kl2,il2)/ro2(jl1,kl2,il2)) * 0.5E+00)
     &             /(2.0E+00 * rrbc * dth0) * bt1p1 * bp1p1
               tcndp1 = aa * tt1p1**2.5E+00
     &                /(br1p1**2 + bt1p1**2 + bp1p1**2 + pg1p1*cbeta)
               aa = (pg3(jl2,kl2,il2)/ro2(jl2,kl2,il2)
     &              -pg3(jl2,kl1,il2)/ro2(jl2,kl1,il2))
     &             / (dph0 * rrbc * sinbc) * (bp1p2**2 + pg1p2*cbeta)
     &            +((pg3(jl2,kl1,il3)/ro2(jl2,kl1,il3)
     &              +pg3(jl2,kl2,il3)/ro2(jl2,kl2,il3)) * 0.5E+00
     &             -(pg3(jl2,kl1,il1)/ro2(jl2,kl1,il1)
     &              +pg3(jl2,kl2,il1)/ro2(jl2,kl2,il1)) * 0.5E+00)
     &             / (dra0 * 2.0E+00)        * br1p2 * bp1p2
     &            +((pg3(jl3,kl1,il2)/ro2(jl3,kl1,il2)
     &              +pg3(jl3,kl2,il2)/ro2(jl3,kl2,il2)) * 0.5E+00
     &             -(pg3(jl1,kl1,il2)/ro2(jl1,kl1,il2)
     &              +pg3(jl1,kl2,il2)/ro2(jl1,kl2,il2)) * 0.5E+00)
     &             / (2.0E+00 * rrbc * dth0) * bt1p2 * bp1p2
               tcndp2 = aa * tt1p2**2.5E+00
     &                / (br1p2**2 + bt1p2**2 + bp1p2**2 + pg1p2*cbeta)

** MIND these two are orthodox !
               if (itcndtype .EQ. 1) then
                 xxpg =((tcndr1*rrb1**2 - tcndr2*rrb2**2)/ dra2
     &                + (tcndt1*sinb1   - tcndt2*sinb2)/(rrbc*dth1)
     &                + (tcndp1         - tcndp2)/ (rrbc*sinbc*dph0))
     &                * tcndfact * (gammax - 1.0E+00)
               else if (itcndtype .EQ. 2) then
                 xxpg =((tcndr1*rrb1**2 - tcndr2*rrb2**2)/ dra2
     &                + (tcndt1*sinb1   - tcndt2*sinb2)/(rrbc*dth1)
     &                + (tcndp1         - tcndp2)/ (rrbc*sinbc*dph0))
     &                * tcndfact * (gammax - 1.0E+00)
     &                * ro2(j,k,i)
               else
** this is not orthodox, but fast and works resonably well.
                 xxpg =((tcndr1 * ro1r1 * rrb1**2
     &                  -tcndr2 * ro1r2 * rrb2**2) / dra2
     &                + (tcndt1 * ro1t1 * sinb1
     &                  -tcndt2 * ro1t2 * sinb2) / (rrbc*dth1)
     &                + (tcndp1 * ro1p1
     &                  -tcndp2 * ro1p2) / (rrbc * sinbc*dph0))
     &                * tcndfact * (gammax - 1.0E+00)
               endif

               if (i .EQ. 0) xxpg = 0.0E+00 ! for the sake of ...
               pg9(j,k,i) = pg3(j,k,i) + dtlocal * xxpg
               
             enddo
             enddo ! end  of j,k-loop,,,, main part
* axis
             do k = 1, kk
               pg9( 0,k,i) = pg9(   1,k,i) ! MIND this is not exactly correct
               pg9(jj,k,i) = pg9(jj-1,k,i)
             enddo
* around axi : ave. along lon.
             if (kk .GT. 3) then
               do j = 0, jj
                 if ((j .EQ. 0) .OR. (j .EQ. jj)) then ! average for all k
                   aa = 0.0E+00
                   do k = 1, kk
                     aa = aa + pg9(j,k,i)
                   enddo
                   aa = aa / float(kk)
                   do k = 1, kk
                     pg9(j,k,i) = aa
                   enddo
                 else
                   kstep2 = kstepb(j)
                   if (kstep2 .GT. 1) then !              average for local k
                     do k = 1, kk, kstep2
                       aa = 0.0E+00
                       do k2 = 0, kstep2 - 1
                         aa = aa + pg9(j,k+k2,i)
                       enddo
                       aa = aa / float(kstep2)
                       do k2 = 0, kstep2 - 1
                         pg9(j,k+k2,i) = aa
                       enddo ! k-local-loop
                     enddo !   k-loop
                   endif
                 endif
               enddo ! end j-loop
             endif ! end if (kk .GT. 3)...
* over rot. axis
             do k = 1, kk
               k2 = mod(k + kk/2,kk)
               if (k2 .EQ. 0) k2 = kk
               pg9(  -1,k,i) = pg9(   1,k2,i)
               pg9(jj+1,k,i) = pg9(jj-1,k2,i)
             enddo
* periodicity along longitude
             do j =-1, jj + 1
               pg9(j,   0,i) = pg9(j,kk  ,i)
               pg9(j,kk+1,i) = pg9(j,   1,i)
             enddo
             do j =-1, jj + 1
               pg9(j,  -1,i) = pg9(j,kk-1,i)
               pg9(j,kk+2,i) = pg9(j,   2,i)
             enddo

             enddo ! end of i - loop
!$omp end parallel do

* copy
!$omp parallel do shared(pg3,pg9)
!$omp&            private(i,k,j)
!$omp&            shared(ii,jj,kk)
             do i = 1, ii + 1
             do k =-1, kk + 2
             do j =-1, jj + 1
               pg3(j,k,i) = pg9(j,k,i)
             enddo
             enddo
             enddo
!$omp end parallel do

           enddo ! end of n-loop, iteration part

       endif ! end if (ltcndct)

* fluid viscosity ! -----------------------------------------
       if (lflwvisc) then
!$omp parallel do private(i,k,j)
!$omp&            shared(ur2,ut2,up2,ur3,ut3,up3)
!$omp&            shared(ii,jj,kk)
         do i = 0, ii + 1
           do k = -1, kk + 2
           do j = -1, jj + 1
             ur3(j,k,i) = ur2(j,k,i)
             ut3(j,k,i) = ut2(j,k,i)
             up3(j,k,i) = up2(j,k,i)
           enddo
           enddo
         enddo
!$omp end parallel do

         dtlocal = 1.0E+30
!$omp parallel do private(i,k,j,aa,bb,cc,dd)
!$omp&            shared(rr,theta,phi,sain,kstepb)
!$omp&            shared(ii,jj,kk)
!$omp&            reduction(min:dtlocal)
         do i = 1, ii
           do k = 1, kk
           do j = 1, jj - 1
             aa = (rr(i+1) - rr(i-1)) * 0.5E+00
             bb = rr(i) * (theta(j+1) - theta(j-1)) * 0.5E+00
             cc = rr(i) * sain(j) * (phi(k+1) - phi(k-1)) * 0.5E+00
     &          * float(kstepb(j))
             dd = min(aa,bb,cc)
             dd = 0.5E+00 * dd**2 / viscfact * 0.80E+00 ! CFL
             if (dd .LT. dtlocal) dtlocal = dd
           enddo
           enddo
         enddo
!$omp end parallel do

         ntime = int(dt / dtlocal) + 1
         if (ncal .EQ. 1) then
           write(*,*) 'Flow Vis --- '
           write(*,*) 'Ntime dt_local(original) dt_global = '
           write(*,*)  ntime, dtlocal, dt
         endif
         dtlocal = dt / float(ntime)

         do n = 1, ntime
!$omp parallel do default(private)
!$omp&            shared(dtlocal)
!$omp&            shared(kstepb)
!$omp&            shared(ro2)
!$omp&            shared(ur3,ut3,up3)
!$omp&            shared(ur9,ut9,up9)
!$omp&            shared(rr,theta,phi,sain,sainb,kosab)
!$omp&            shared(ii,jj,kk)
           do i = 0, ii + 1
             do k = 1, kk
             do j = 1, jj - 1 ! MIND the range...

               il1 = i - 1
               il2 = i
               il3 = i + 1
               if (il1 .LT.    0) il1 =    0
               if (il3 .GT. ii+1) il3 = ii+1
               jl1 = j - 1
               jl2 = j
               jl3 = j + 1
               if (jl1 .LT.   -1) jl1 =   -1
               if (jl3 .GT. jj+1) jl3 = jj+1
               kl1 = k -  kstepb(j) + kk * 2
               kl1 = mod(kl1,kk)
               kl2 = k
               kl3 = k +  kstepb(j) + kk * 2
               kl3 = mod(kl3,kk)
*               kl1 = k - 1
*               kl2 = k
*               kl3 = k + 1

               rrb1 =(rr(il3) + rr(il2)) * 0.5E+00
               rrb2 =(rr(il1) + rr(il2)) * 0.5E+00
               rrbc = rr(il2)
               dra0 =  rrb1 - rrb2
               dra1 = (rrb1**2 - rrb2**2) * 0.5E+00
               dra2 = (rrb1**3 - rrb2**3) / 3.0E+00
               if (j .EQ. 0) then
                 sinb1 = sainb(j)
                 cosb1 = kosab(j)
                 sinb2 = 0.0E+00
                 cosb2 = 1.0E+00
                 sinbc = (cosb2 - cosb1) / (pi / float(jj) * 0.5E+00)
               else if (j .EQ. jj) then
                 sinb1 = 0.0E+00
                 cosb1 =-1.0E+00
                 sinb2 = sainb(j-1)
                 cosb2 = kosab(j-1)
                 sinbc = (cosb2 - cosb1) / (pi / float(jj) * 0.5E+00)
               else
                 sinb1 = sainb(j)
                 cosb1 = kosab(j)
                 sinb2 = sainb(j-1)
                 cosb2 = kosab(j-1)
                 sinbc = (cosb2 - cosb1) / (pi / float(jj))
               endif
               cotbc = (sinb1 - sinb2) / (cosb2 - cosb1)
               dth0 =           pi / float(jj) * rr(i)
               dth1 =  (cosb2 - cosb1)          * rr(i)
               dph0 = 2.0E+00 * pi / float(kk) * rr(i) * sinbc
               dph2 = dph0
               dph0 = dph0 * float(kstepb(j))
               dph2 = dph2 * float(kstepb(j))

               ro1r1 = (ro2(j,k,il3) + ro2(j,k,i)) * 0.5E+00
               ro1r2 = (ro2(j,k,il1) + ro2(j,k,i)) * 0.5E+00
               ro1t1 = (ro2(jl3,k,i) + ro2(j,k,i)) * 0.5E+00
               ro1t2 = (ro2(jl1,k,i) + ro2(j,k,i)) * 0.5E+00
               ro1p1 = (ro2(j,kl3,i) + ro2(j,k,i)) * 0.5E+00
               ro1p2 = (ro2(j,kl1,i) + ro2(j,k,i)) * 0.5E+00

               ur1r1 = (ur3(j,k,il3) + ur3(j,k,i)) * 0.5E+00
               ur1r2 = (ur3(j,k,il1) + ur3(j,k,i)) * 0.5E+00
               ut1r1 = (ut3(j,k,il3) + ut3(j,k,i)) * 0.5E+00
               ut1r2 = (ut3(j,k,il1) + ut3(j,k,i)) * 0.5E+00
               up1r1 = (up3(j,k,il3) + up3(j,k,i)) * 0.5E+00
               up1r2 = (up3(j,k,il1) + up3(j,k,i)) * 0.5E+00
               ur1t1 = (ur3(jl3,k,i) + ur3(j,k,i)) * 0.5E+00
               ur1t2 = (ur3(jl1,k,i) + ur3(j,k,i)) * 0.5E+00
               ut1t1 = (ut3(jl3,k,i) + ut3(j,k,i)) * 0.5E+00
               ut1t2 = (ut3(jl1,k,i) + ut3(j,k,i)) * 0.5E+00
               up1t1 = (up3(jl3,k,i) + up3(j,k,i)) * 0.5E+00
               up1t2 = (up3(jl1,k,i) + up3(j,k,i)) * 0.5E+00
               ur1p1 = (ur3(j,kl3,i) + ur3(j,k,i)) * 0.5E+00
               ur1p2 = (ur3(j,kl1,i) + ur3(j,k,i)) * 0.5E+00
               ut1p1 = (ut3(j,kl3,i) + ut3(j,k,i)) * 0.5E+00
               ut1p2 = (ut3(j,kl1,i) + ut3(j,k,i)) * 0.5E+00
               up1p1 = (up3(j,kl3,i) + up3(j,k,i)) * 0.5E+00
               up1p2 = (up3(j,kl1,i) + up3(j,k,i)) * 0.5E+00

               if ((i .GT. 0) .AND. (i .LE. ii)) then

                 aa=((ur3(jl2,kl2,il3)-ur3(jl2,kl2,il2))*ro1r1*rrb1**2
     &              -(ur3(jl2,kl2,il2)-ur3(jl2,kl2,il1))*ro1r2*rrb2**2)
     &             / dra0 / dra2
                 bb=((ur3(jl3,kl2,il2)-ur3(jl2,kl2,il2))*ro1t1*sinb1
     &              -(ur3(jl2,kl2,il2)-ur3(jl1,kl2,il2))*ro1t2*sinb2)
     &             / dth0 / dth1
                 cc=((ur3(jl2,kl3,il2)-ur3(jl2,kl2,il2))*ro1p1
     &              -(ur3(jl2,kl2,il2)-ur3(jl2,kl1,il2))*ro1p2)
     &             / dph0 / dph2
                 dd= - 2.0E+00 * ur3(jl2,kl2,il2) / rr(i)**2
     &               - 2.0E+00 * ut3(jl2,kl2,il2) / rr(i)**2 * cotbc
                 ee= - 2.0E+00 * ((ut1t1 - ut1t2) / dth0
     &                           +(up1p1 - up1p2) / dph0)
                 if ((j .EQ. 0) .OR. (j .EQ. jj)) then
                   bb = bb * 0.5E+00
                   cc = cc * 0.5E+00
                   ee = ee * 0.5E+00
                 endif
                 xxur = viscfact
     &                *((aa + bb + cc) / ro2(jl2,kl2,il2)
     &                 +(dd + ee))

                 aa=((ut3(jl2,kl2,il3)-ut3(jl2,kl2,il2))*ro1r1*rrb1**2
     &              -(ut3(jl2,kl2,il2)-ut3(jl2,kl2,il1))*ro1r2*rrb2**2)
     &             / dra0 / dra2
                 bb=((ut3(jl3,kl2,il2)-ut3(jl2,kl2,il2))*ro1t1*sinb1
     &              -(ut3(jl2,kl2,il2)-ut3(jl1,kl2,il2))*ro1t2*sinb2)
     &             / dth0 / dth1
                 cc=((ut3(jl2,kl3,il2)-ut3(jl2,kl2,il2))*ro1p1
     &              -(ut3(jl2,kl2,il2)-ut3(jl2,kl1,il2))*ro1p2)
     &             / dph0 / dph2
                 dd= -ut3(jl2,kl2,il2) / rr(i)**2 / sinbc**2
                 ee=   2.0E+00 * ((ur1t1 - ur1t2) / dth0
     &                           -(up1p1 - up1p2) / dph0 * cotbc)
                 if ((j .EQ. 0) .OR. (j .EQ. jj)) then
                   bb = bb * 0.5E+00
                   cc = cc * 0.5E+00
                   ee = ee * 0.5E+00
                 endif
                 xxut = viscfact
     &                *((aa + bb + cc) / ro2(jl2,kl2,il2)
     &                 +(dd + ee))

                 aa=((up3(jl2,kl2,il3)-up3(jl2,kl2,il2))*ro1r1*rrb1**2
     &              -(up3(jl2,kl2,il2)-up3(jl2,kl2,il1))*ro1r2*rrb2**2)
     &             / dra0 / dra2
                 bb=((up3(jl3,kl2,il2)-up3(jl2,kl2,il2))*ro1t1*sinb1
     &              -(up3(jl2,kl2,il2)-up3(jl1,kl2,il2))*ro1t2*sinb2)
     &             / dth0 / dth1
                 cc=((up3(jl2,kl3,il2)-up3(jl2,kl2,il2))*ro1p1
     &              -(up3(jl2,kl2,il2)-up3(jl2,kl1,il2))*ro1p2)
     &             / dph0 / dph2
                 dd= -up3(jl2,kl2,il2) / rr(i)**2 / sinbc**2
                 ee=   2.0E+00 * ((ur1p1 - ur1p2) / dph0
     &                           +(ut1p1 - ut1p2) / dph0 * cotbc)
                 if ((j .EQ. 0) .OR. (j .EQ. jj)) then
                   bb = bb * 0.5E+00
                   cc = cc * 0.5E+00
                   ee = ee * 0.5E+00
                 endif
                 xxup = viscfact
     &                *((aa + bb + cc) / ro2(jl2,kl2,il2)
     &                 +(dd + ee))
               else
                 xxur = 0.0E+00
                 xxut = 0.0E+00
                 xxup = 0.0E+00
               endif
               if (i .EQ. 0) then
                 xxur = 0.0E+00
                 xxut = 0.0E+00
                 xxup = 0.0E+00
               endif
               ur9(j,k,i) = ur3(j,k,i) + dtlocal * xxur
               ut9(j,k,i) = ut3(j,k,i) + dtlocal * xxut
               up9(j,k,i) = up3(j,k,i) + dtlocal * xxup
             enddo
             enddo ! end of j,k-loop

* around axi : ave. along lon.
             if (kk .GT. 3) then
               do j = 1, jj - 1
                 kstep2 = kstepb(j)
                 if (kstep2 .GT. 1) then !              average for local k
                   do k = 1, kk, kstep2
                     aa = 0.0E+00
                     bb = 0.0E+00
                     cc = 0.0E+00
                     do k2 = 0, kstep2 - 1
                       aa = aa + ur9(j,k+k2,i)
                       bb = bb + ut9(j,k+k2,i)
                       cc = cc + up9(j,k+k2,i)
                     enddo
                     aa = aa / float(kstep2)
                     bb = bb / float(kstep2)
                     cc = cc / float(kstep2)
                     do k2 = 0, kstep2 - 1
                       ur9(j,k+k2,i) = aa
                       ut9(j,k+k2,i) = bb
                       up9(j,k+k2,i) = cc
                     enddo ! k-local-loop
                   enddo !   k-loop
                 endif
               enddo ! end j-loop
             endif ! end if (kk .GT. 3)...
* over rot. axis
             do k = 1, kk
               k2 = mod(k + kk/2,kk)
               if (k2 .EQ. 0) k2 = kk
               ur9(  -1,k,i) = ur9(   1,k2,i)
               ur9(jj+1,k,i) = ur9(jj-1,k2,i)
               ut9(  -1,k,i) = ut9(   1,k2,i)
               ut9(jj+1,k,i) = ut9(jj-1,k2,i)
               up9(  -1,k,i) = up9(   1,k2,i)
               up9(jj+1,k,i) = up9(jj-1,k2,i)
             enddo
* axis
             do k = 1, kk
               ur9( 0,k,i) = ur9(   1,k,i) ! MIND this is not exactly correct
               ur9(jj,k,i) = ur9(jj-1,k,i)
               ut9( 0,k,i) = ut9(   1,k,i)
               ut9(jj,k,i) = ut9(jj-1,k,i)
               up9( 0,k,i) = up9(   1,k,i)
               up9(jj,k,i) = up9(jj-1,k,i)
             enddo
* periodicity along longitude
             do j =-1, jj + 1
               ur9(j,   0,i) = ur9(j,kk  ,i)
               ur9(j,kk+1,i) = ur9(j,   1,i)
               ut9(j,   0,i) = ut9(j,kk  ,i)
               ut9(j,kk+1,i) = ut9(j,   1,i)
               up9(j,   0,i) = up9(j,kk  ,i)
               up9(j,kk+1,i) = up9(j,   1,i)
             enddo
             do j =-1, jj + 1
               ur9(j,  -1,i) = ur9(j,kk-1,i)
               ur9(j,kk+2,i) = ur9(j,   2,i)
               ut9(j,  -1,i) = ut9(j,kk-1,i)
               ut9(j,kk+2,i) = ut9(j,   2,i)
               up9(j,  -1,i) = up9(j,kk-1,i)
               up9(j,kk+2,i) = up9(j,   2,i)
             enddo

           enddo ! end of i - loop
!$omp end parallel do

* copy
!$omp parallel do shared(ur3,ut3,up3,ur9,ut9,up9)
!$omp&            private(i,k,j)
!$omp&            shared(ii,jj,kk)
           do i = 1, ii + 1
             do k = -1, kk + 2
             do j = -1, jj + 1
               ur3(j,k,i) = ur9(j,k,i)
               ut3(j,k,i) = ut9(j,k,i)
               up3(j,k,i) = up9(j,k,i)
             enddo
             enddo
           enddo
!$omp end parallel do

         enddo ! end of n-local loop

       endif

* move to common main block.
       if (ltcndct) then
!$omp parallel do private(i,k,j)
!$omp&            shared(pg2,pg3)
!$omp&            shared(ii,jj,kk)
         do i = 1, ii + 1 ! MIND the range
           do k = -1, kk + 2
           do j = -1, jj + 1
             pg2(j,k,i) = pg3(j,k,i)
           enddo
           enddo
         enddo
!$omp end parallel do
       endif

       if (lflwvisc) then
!$omp parallel do private(i,k,j)
!$omp&            shared(ur2,ut2,up2,ur3,ut3,up3)
!$omp&            shared(ii,jj,kk)
         do i = 1, ii + 1 ! MIND the range
           do k = -1, kk + 2
           do j = -1, jj + 1
             ur2(j,k,i) = ur3(j,k,i)
             ut2(j,k,i) = ut3(j,k,i)
             up2(j,k,i) = up3(j,k,i)
           enddo
           enddo
         enddo
!$omp end parallel do
       endif
*
       return
       endsubroutine


*--------------------------------------------------------------------
*
* 1) Scheme in the domain of computation.
*   a) Roe's TVD(?) + MUSCL modificaiton.
*   b) integer "iprc" stands for the first or second step,
*        currently in 1 step method means nothing
*
* 2) Determine variables on the rotational axis
*   a) (pseudo-) Finit Volume Method
*      Use the relation obtained by integrating the differential
*        equation over the sherical segment near the rotational axis
*      The relation obtained by integrating the each differential equation
*        over the small plane parpendicular to the rotational axis.
*   b) extrapolate by avaraging values
*        at the point nearest to the rotational axis.
*   c) Longitudinal and latitudinal component of vector variables
*         are fixed to be null (assumption & geometrial requirements).
*
*--------------------------------------------------------------------
*
       subroutine prc1st(ii,jj,kk,
     &                   dt,iprc,omega,gammax,kstep,vcrtrn,mcrtrn,
     &                   denf0,tmpf0,gamma0,
     &                   v0,r0,p0,b0,t0,rhoc,n0,ncal,lstop,
     &                   ro0,pg0,            br0,bt0,bp0,
     &                   ro1,pg1,ur1,ut1,up1,br1,bt1,bp1,
     &                   ro2,pg2,ur2,ut2,up2,br2,bt2,bp2,
     &                       pg3,            br3,bt3,bp3,
     &                   ro9,pg9,ur9,ut9,up9,br9,bt9,bp9,
     &                   ed1,mr1,mt1,mp1,ed2,mr2,mt2,mp2,
     &                   rr,sain,kosa,sainb,kosab,
     &                   ropark0,pgpark0,roparkb,pgparkb,
     &                   sinph,cosph,sinphb,cosphb,dbr0,dbt0,dbp0)
       implicit none
* arguments
       intent(in)    ::  ii,jj,kk
       intent(in)    ::  dt,iprc,omega,gammax,kstep,vcrtrn,mcrtrn
       intent(in)    ::  denf0,tmpf0,gamma0
       intent(in)    ::  v0,r0,p0,b0,t0,rhoc,n0,ncal
       intent(out)   ::                              lstop
       intent(in)    ::  ro0,pg0,            br0,bt0,bp0
       intent(in)    ::  ro1,pg1,ur1,ut1,up1,br1,bt1,bp1
       intent(inout) ::  ro2,pg2,ur2,ut2,up2,br2,bt2,bp2
       intent(out)   ::      pg3,            br3,bt3,bp3
       intent(out)   ::  ro9,pg9,ur9,ut9,up9,br9,bt9,bp9
       intent(out)   ::  ed1,mr1,mt1,mp1,ed2,mr2,mt2,mp2
       intent(in)    ::  rr,sain,kosa,sainb,kosab
       intent(in)    ::  ropark0,pgpark0,roparkb,pgparkb
       intent(in)    ::  sinph,cosph,sinphb,cosphb,dbr0,dbt0,dbp0
       integer ii, jj, kk
*       parameter(ii = 72, jj = 32, kk = 64)
       real    dt
       integer iprc
       real    omega, gammax, gamma0
       integer kstep(0:jj)
       real    vcrtrn, mcrtrn
       real    v0, r0, p0, b0, t0, rhoc, n0
       integer ncal
       logical lstop
       real    denf0(0:jj,0:kk+1)
       real    tmpf0(0:jj,0:kk+1)
       real    rr(-2:ii+3)
       real    sain(0:jj), kosa(0:jj)
       real    sainb(0:jj-1), kosab(0:jj-1)
       real    sinphb(0:kk), cosphb(0:kk)
       real    sinph(0:kk+1), cosph(0:kk+1)
       real    ro1(-2:jj+2,-1:kk+2,-2:ii+3),pg1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur1(-2:jj+2,-1:kk+2,-2:ii+3),ut1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up1(-2:jj+2,-1:kk+2,-2:ii+3),br1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt1(-2:jj+2,-1:kk+2,-2:ii+3),bp1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ro9(-2:jj+2,-1:kk+2,-2:ii+3),pg9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur9(-2:jj+2,-1:kk+2,-2:ii+3),ut9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up9(-2:jj+2,-1:kk+2,-2:ii+3),br9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt9(-2:jj+2,-1:kk+2,-2:ii+3),bp9(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ro2(-2:jj+2,-1:kk+2,-2:ii+3),pg2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ur2(-2:jj+2,-1:kk+2,-2:ii+3),ut2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    up2(-2:jj+2,-1:kk+2,-2:ii+3),br2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt2(-2:jj+2,-1:kk+2,-2:ii+3),bp2(-2:jj+2,-1:kk+2,-2:ii+3)
       real                                 pg3(-2:jj+2,-1:kk+2,-2:ii+3)
       real                                 br3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    bt3(-2:jj+2,-1:kk+2,-2:ii+3),bp3(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ed1(-2:jj+2,-1:kk+2,-2:ii+3),mr1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    mt1(-2:jj+2,-1:kk+2,-2:ii+3),mp1(-2:jj+2,-1:kk+2,-2:ii+3)
       real    ed2(-2:jj+2,-1:kk+2,-2:ii+3),mr2(-2:jj+2,-1:kk+2,-2:ii+3)
       real    mt2(-2:jj+2,-1:kk+2,-2:ii+3),mp2(-2:jj+2,-1:kk+2,-2:ii+3)
       integer ipho
       parameter(ipho = 5)
       real    ro0(0:jj,1:kk,-1:ipho),pg0(0:jj,1:kk,-1:ipho)
       real    br0(0:jj,1:kk,-1:ipho)
       real    bt0(0:jj,1:kk,-1:ipho),bp0(0:jj,1:kk,-1:ipho)
       real    dbr0(0:jj,0:kk+1,-1:ipho)
       real    dbt0(0:jj,0:kk+1,-1:ipho)
       real    dbp0(0:jj,0:kk+1,-1:ipho)
       real    pgpark0(-2:ii+3), ropark0(-2:ii+3)
       real    pgparkb(-2:ii+2), roparkb(-2:ii+2)
* local
       integer i, j, k, l, m, n, idrct, idrctend, iloc
       integer j2, k2
       real    aa, bb, cc, dd, ee, ff, gg, hh ! , qq, xx
       logical ldummy
       integer idrct2
       integer nabnorm
       integer kl1, kl2, kl3, kl4
* geometrics and volume-controllings
       real    dra0, dra1, dra2
       real    dth0, dth1
       real    dph0, dph2
       real    ddra0, ddra1, ddra2 ! inversed ones
       real    ddth0, ddth1
       real    ddph0, ddph2
       real    sin0, cos0, cotbc
       real    rrb1, rrb2, rrbc, rrbc2nd, drr1, drr2
       real    rr2ndmod, rr1stmod
       real    sinb1, sinb2, sinbc, cosb1, cosb2, cosbc
       real    cosphb1, cosphb2, sinphb1, sinphb2, cosph0, sinph0
       real    sindth2, sindph2
       real    cosdth2, cosdph2
       integer kstepb(0:jj)
* locals for compatibility relation
       real    dtlocal
       real    rhsall(7,0:jj,1:kk), rhsall2(7)
       real    dbrdt(0:jj,1:kk) ! for upper.
       real    matl74(7,4), matl33(3,3)
       real    ansvec1(3), ansvec2(3), ansvec3(3)
       logical lplasmod
       real    dpgdtmod, drodtmod
* mass flux limits in terms of magnetic field strength
       logical lmflxvar
       parameter(lmflxvar = .false.)
       real    vcrtrn2, mcrtrn2
* some flag
       logical lfarmnl  ! Near Magnetic neutral Line ...
       logical lweakmag ! Magnetic field strength is week (compared to Pg)
* functions defined by myself, for MUSCL and other limiters
       interface
         pure function nonmscl3(dy1,dy2,dy3,dx1,dx2,dx3)
           implicit none
           real        nonmscl3
           intent(in) ::        dy1,dy2,dy3,dx1,dx2,dx3
           real    dy1,dy2,dy3
           real    dx1,dx2,dx3
         endfunction
         pure function mscl3rd(dyf,dyc,dx1,dx2)
           implicit none
           real        mscl3rd
           intent(in)  ::      dyf,dyc,dx1,dx2
           real    dyf, dyc
           real    dx1, dx2
         endfunction
         pure function weno3rd(dyf,dyc,dx1,dx2)
           implicit none
           real        weno3rd
           intent(in)  ::      dyf,dyc,dx1,dx2
           real    dyf, dyc
           real    dx1, dx2
         endfunction
       endinterface
* locals for Jacobian Matrix, TVD and MUSCL parts
!       real    dlefrig(7,2)
       real    drestmg(2)
       real    d3r, d2r ! , d1r
       real    d1l, d2l ! , d3l
* TVD, MUSCL for MHD
!       real    modcvar7(7,2)
!       real    reigemat(7,7)
       real    leigemat(7,7)
!       real    reigemtl(7,7), leigemtl(7,7)
!       real    reigemtr(7,7), leigemtr(7,7)
!       real    matabs(7,7), matabsb(7,7), matabsc(7,7)
       real    tvdcr(7,0:1), tvdct(7,0:1), tvdcp(7,0:1)

* numerical flux
       real    tnsrr(-1:0),tnstr(-1:0),tnspr(-1:0) ! diadic tensor
       real    tnsrt(-1:0),tnstt(-1:0),tnspt(-1:0)
       real    tnsrp(-1:0),tnstp(-1:0),tnspp(-1:0)
* loop-locals for globals
       integer ncal2, iprc2
       real    gammax2, dgammax2
* variables
       real    ro1a, ro1b, ro1c, ro1d
       real    mr1a, mr1b, mr1c, mr1d
       real    mt1a, mt1b, mt1c, mt1d
       real    mp1a, mp1b, mp1c, mp1d
       real    br1a, br1b, br1c, br1d
       real    bt1a, bt1b, bt1c, bt1d
       real    bp1a, bp1b, bp1c, bp1d
       real    ed1a, ed1b, ed1c, ed1d
       real    pg1a, pg1b, pg1c, pg1d
       real    ur1a, ur1b, ur1c, ur1d
       real    ut1a, ut1b, ut1c, ut1d
       real    up1a, up1b, up1c, up1d
       real    rod1, rod2, rod3
       real    mrd1, mrd2, mrd3
       real    mtd1, mtd2, mtd3
       real    mpd1, mpd2, mpd3
       real    brd1, brd2, brd3
       real    btd1, btd2, btd3
       real    bpd1, bpd2, bpd3
       real    m1d1, m1d2, m1d3
       real    m2d1, m2d2, m2d3
       real    m3d1, m3d2, m3d3
       real    b1d1, b1d2, b1d3
       real    b2d1, b2d2, b2d3
       real    b3d1, b3d2, b3d3
       real    edd1, edd2, edd3
       real    pgd1, pgd2, pgd3
       real    urd1, urd2, urd3
       real    utd1, utd2, utd3
       real    upd1, upd2, upd3
       real    dx1, dx2, dx3
       real    mra, mta, mpa, ura, uta, upa
       real    bra, bta, bpa, roa, pga, eda
       real    b1, b2, b3, v1, v2, v3
* increments
       real    xx2, xx3, xx4, xx5, xx6
       real    xxxx1, xxyy1, xxzz1, xxrr1, xxtt1, xxpp1
       real    xxxx2, xxyy2, xxzz2, xxrr2, xxtt2, xxpp2
* stored dU/dT : HD and MHD
       real    xxro, xxed, xxpg
       real    xxbr, xxbt, xxbp
       real    xxmr, xxmt, xxmp
       real    xxrol(2), xxedl(2)
       real    xxbrl(2), xxbtl(2), xxbpl(2)
       real    xxmrl(2), xxmtl(2), xxmpl(2)
       real    roltmp, pgltmp, edltmp
       real    urltmp, utltmp, upltmp
       real    brltmp, btltmp, bpltmp
       real    mrltmp, mtltmp, mpltmp
       real    xxpghd
* values at cell interface(s)
       real    rob(0:1,2,3),pgb(0:1,2,3),edb(0:1,2,3)!\pm 1/2*l,L or R,Direction
       real    urb(0:1,2,3),utb(0:1,2,3),upb(0:1,2,3)
       real    mrb(0:1,2,3),mtb(0:1,2,3),mpb(0:1,2,3)
       real    brb(0:1,2,3),btb(0:1,2,3),bpb(0:1,2,3)
* others
       real    ro1r1, ro1r2, ro1t1, ro1t2, ro1p1, ro1p2, ro1cc
       real    pg1r1, pg1r2, pg1t1, pg1t2, pg1p1, pg1p2, pg1cc
       real    mr1r1, mr1r2, mr1t1, mr1t2, mr1p1, mr1p2, mr1cc
       real    mt1r1, mt1r2, mt1t1, mt1t2, mt1p1, mt1p2, mt1cc
       real    mp1r1, mp1r2, mp1t1, mp1t2, mp1p1, mp1p2, mp1cc
       real    ur1r1, ur1r2, ur1t1, ur1t2, ur1p1, ur1p2, ur1cc
       real    ut1r1, ut1r2, ut1t1, ut1t2, ut1p1, ut1p2, ut1cc
       real    up1r1, up1r2, up1t1, up1t2, up1p1, up1p2, up1cc
       real    br1r1, br1r2, br1t1, br1t2, br1p1, br1p2, br1cc
       real    bt1r1, bt1r2, bt1t1, bt1t2, bt1p1, bt1p2, bt1cc
       real    bp1r1, bp1r2, bp1t1, bp1t2, bp1p1, bp1p2, bp1cc
*
       real    mgstr, ss
       real    jr1cc, jt1cc, jp1cc
       real    divmag
       real    gravit
*
* impose gradual increase of Br at surface
!       logical ldmagon
!       parameter(ldmagon = .false.)
*
* switch & option for scheme.
*
* ldistant
       logical ldistant
       parameter(ldistant = .false.)
*
* flags for differenced equations.
       logical lconene !             Conservation form of Energy eq. or not.
       parameter(lconene = .true.) ! If false, axi-inner part will crash....(as of 2008/08/13)
       logical lmodcone
       parameter(lmodcone = .false.) ! modify B-part
       real    modmagen
       logical lvolmom
       parameter(lvolmom = .true.) ! FVM like or old version for surpus term of dm/dt
       logical ldiagmom
       parameter(ldiagmom = .true.) ! Diag parts (Pg + B^2/2) of Mom. will be treated as scalar
       logical ljxb
       parameter(ljxb = .true.) ! B in mom.eq will be treated with non-conservation form, jxB
       logical lvolmag
       parameter(lvolmag = .true.) ! FVM like or not         for surpus term of dB/dt
* flag for Characteristic method at inner surface
!       logical lchkchar
!       parameter(lchkchar = .false.) ! check whetehr used eigen-system be adequate
       logical lchrsurf
       parameter(lchrsurf = .true.) ! characteristic method at surface, should be false for distant
       integer ichrsurf
       parameter(ichrsurf = 0) ! v > v_c (0,1--7)= no special, specified choices
       real     gammab !             polytrope index for V_r < 0 ; usually be same as gammax
       real     gamma6  !            polytrope index for ich is 6
       parameter(gamma6 = 1.0E+00) ! usually be 1.00 for fixed T
!       real     gamma7  !            polytrope index for ich is 7
!       parameter(gamma7 = 1.0E+00) ! usually be 1.00 for fixed T
       real     gamma2  !            polytrope index for ich is 2
       parameter(gamma2 = 1.0E+00) ! usually be 1.00 for fixed T
       logical lnorefbot !     if true, provisional Mr will be evaluated with non-ref
       parameter(lnorefbot = .false.) ! if true, ichrsurf should be -1
       real     mrdummy
       integer  iclass ! 0 ; zeor, 1 : normal, otherwise : V > Vc, or M > Mc
       real     mrzero
       parameter(mrzero = 1.0E-00) ! mind mrzero = 1 / rr**2 or same order.
* dB/dt by CT
!       logical lct
!       parameter(lct = .false.)
* switch to guarantee minimum N and P when stagnant
       integer  imindenp !       0 = positive drho,dpg, 1 = old-type
       parameter(imindenp = 1) ! 2 = new-type, otherwise = None
* switch to change surface plasma map
       logical lmodmap
       parameter(lmodmap = .false.) ! MIND if true, ichrsurf should be 0 or 6.
       real    timemod
       parameter(timemod = 2.0E+00) !  time in hour
* flag for TVD
       logical ltvd
       parameter(ltvd = .true.) ! do tvd or not
       integer itvd
       parameter(itvd = 3) ! 1 : normal TVD, 2 FVS-lile, 3 = LaxFriedrichs, otherwise = Lapidus
* flag for MUSCL interpolation
       logical lmuscl
       parameter(lmuscl = .true.)
!       integer imuscl !        0 = chara, 1=chara with ave,
!       parameter(imuscl = 2) ! 2=cons/prim-Park, otherwize= test. ! if itvd is 3, set this 2
!       logical lmmscl2         ! special choice when imuscl is set 2
!       parameter(lmmscl2 = .true.) ! momentum or bulk speed interpolation : false maybe good for distant
       logical lradpark, lradprk2
       parameter(lradpark = .true.) ! MIND that if itvd is 3, this must be .true.
!       logical lradlog !              additional choice to use log instead of parker-sol.
!       parameter(lradlog = .true.) ! lradpark must be true and then will be overriden.
       integer ib1intp        ! 1= linear, 2= 3rd order non-MUSCL, 3= Moc-like v1 x d(b1) not well
       parameter(ib1intp = 2) ! 4=mix of MUSCL/non-M 3rd, 5= 3rd with limiter, others = weno3rd
* flags for HD plasma
       logical lhdplsm0
       parameter(lhdplsm0 = .false.) ! this must be true when idts is set 1.
!       logical lhdplsm
       integer idts, idte, idt
       parameter(idts = 2) ! MIND this must be set 1 for using HD part
       parameter(idte = 2)
* flag for enforcing EMF and/or v//B condition  at i = 0, j+1/2, k+1/2 etc or not
       logical lvxbzero
       parameter(lvxbzero = .true.) ! may be true unless distant-case (or varying Br is tested).
       logical lvxbzerc
       parameter(lvxbzerc = .true.) ! for average at cell interfaces and some of cell centers
* flag for tweak....
       logical lenfposi ! whether or not enforcing positivity of density and pressure
       parameter(lenfposi = .true.)
* flag for Polytrope..
*       integer ipoly !        MIND when not 0, ichrsurf should be 1 or 0
*       parameter(ipoly = 0) ! 0/1/other : None, P=rho^gamma, rho=P^(1/gamma)
* switch for extrapolation of difference across axi
       integer iextpole
       parameter(iextpole = 3) ! 0=test, 1=0th extrapolation, 2=fixed at zero, otherwize=use ghost
* some parameters for MUSCL : minmod(or similar limiter) or usual 3rd order interplolation
       real    sqrpg
       real    gradlow, gradfact
       parameter(gradlow = 0.01E+00, gradfact = 100.0E+00) ! for 15-zone sine
       real    gradsmll
       parameter(gradsmll = 0.005E+00) ! to prevent dividing by zero
* flag for polar treatment
       logical lpoleave
       parameter(lpoleave = .false.) ! MIND this appear at axibnd()
       real    tempole0
       real    brncc, btncc, bpncc
* criteria of Br/B
       real    brbb
       real    brabsb
       parameter(brabsb = 1.0E-01)
* constant(s)
       real    pi
       parameter(pi = 3.1415926535897932385E+00)

       lstop = .false.

* surpress warning
       aa = v0
       aa = r0
       aa = p0
       aa = b0
       aa = t0
       aa = n0
       aa = rhoc
       aa = denf0(1,1)

* choice !
         gammab = gammax
*       gammab = 1.0E+00 ! if T-map is used, or given (elena's map and hot spot test etc)

       tempole0 = pg0(0,1,0) / ro0(0,1,0)

* show setting
       if (ncal .EQ. 1) then
         write(*,*) ' lconene  = ', lconene
!         write(*,*) ' lmodcone = ', lmodcone
         write(*,*) ' lvolmom  = ', lvolmom
         write(*,*) ' ldiagmom = ', ldiagmom
         write(*,*) ' ljxb     = ', ljxb
         write(*,*) ' lvolmag  = ', lvolmag
         write(*,*) ' lchrsurf = ', lchrsurf
         write(*,*) ' ichrsurf = ', ichrsurf
         write(*,*) ' gammab   = ', gammab
         write(*,*) ' gamma6   = ', gamma6
!         write(*,*) ' gamma7   = ', gamma7
         write(*,*) ' lnonref  = ', lnorefbot
         write(*,*) ' mrzero   = ', mrzero
         write(*,*) ' imindenp = ', imindenp
         write(*,*) ' lmodmap  = ', lmodmap
         write(*,*) ' timemod  = ', timemod
         write(*,*) ' ltvd     = ', ltvd
         write(*,*) ' itvd     = ', itvd
         write(*,*) ' lmuscl   = ', lmuscl
!         write(*,*) ' imuscl   = ', imuscl
         write(*,*) ' lradpark = ', lradpark
         write(*,*) ' lvxbzero = ', lvxbzero
         write(*,*) ' lvxbzerc = ', lvxbzerc
         write(*,*) ' iextpole = ', iextpole
         write(*,*) ' gradlow  = ', gradlow
         write(*,*) ' gradfact = ', gradfact
         write(*,*) ' gradeps  = ', gradsmll
         write(*,*) ' brabsb   = ', brabsb

         write(19,*) ' lconene  = ', lconene
!         write(19,*) ' lmodcone = ', lmodcone
         write(19,*) ' lvolmom  = ', lvolmom
         write(19,*) ' ldiagmom = ', ldiagmom
         write(19,*) ' ljxb     = ', ljxb
         write(19,*) ' lvolmag  = ', lvolmag
         write(19,*) ' lchrsurf = ', lchrsurf
         write(19,*) ' ichrsurf = ', ichrsurf
         write(19,*) ' gammab   = ', gammab
         write(19,*) ' gamma6   = ', gamma6
!         write(19,*) ' gamma7   = ', gamma7
         write(19,*) ' lnonref  = ', lnorefbot
         write(19,*) ' mrzero   = ', mrzero
         write(19,*) ' imindenp = ', imindenp
         write(19,*) ' lmodmap  = ', lmodmap
         write(19,*) ' timemod  = ', timemod
         write(19,*) ' ltvd     = ', ltvd
         write(19,*) ' itvd     = ', itvd
         write(19,*) ' lmuscl   = ', lmuscl
!         write(19,*) ' imuscl   = ', imuscl
         write(19,*) ' lradpark = ', lradpark
         write(19,*) ' lvxbzero = ', lvxbzero
         write(19,*) ' lvxbzerc = ', lvxbzerc
         write(19,*) ' iextpole = ', iextpole
         write(19,*) ' gradlow  = ', gradlow
         write(19,*) ' gradfact = ', gradfact
         write(19,*) ' gradeps  = ', gradsmll
         write(19,*) ' brabsb   = ', brabsb
       endif

       if ((lmodmap) .AND.
     &      (.NOT. ((ichrsurf .EQ. 6) .OR. (ichrsurf .EQ. 0)))) then ! temp.
*     &      (.NOT. ((ichrsurf .EQ. 3) .OR. (ichrsurf .EQ. 0)))) then ! density
         write(*,*) ' choice of ichrsurf and lmodmap is not consistent'
         stop
       endif

!       if (lvolmom) then
         do j = 0, jj
           kstepb(j) = kstep(j)
         enddo
!       else
!       endif

!       if (lvolmom) then
         sindth2 = 0.0E+00
         cosdth2 = 0.0E+00
         sindph2 = 0.0E+00
         cosdph2 = 0.0E+00
!       else !  be sure kstepb() = 1
!       endif

* conservative variables
!$omp parallel do private(i,k,j)
       do i = -2, ii + 3
         do k = -1, kk + 2
         do j = -2, jj + 2
           mr1(j,k,i) = ur1(j,k,i) * ro1(j,k,i)
           mt1(j,k,i) = ut1(j,k,i) * ro1(j,k,i)
           mp1(j,k,i) = up1(j,k,i) * ro1(j,k,i)
           ed1(j,k,i) = pg1(j,k,i) /(gammax - 1.0E+00)
     &                +(ur1(j,k,i)**2 + ut1(j,k,i)**2
     &                 +up1(j,k,i)**2)* ro1(j,k,i) * 0.5E+00
     &                +(br1(j,k,i)**2 + bt1(j,k,i)**2
     &                 +bp1(j,k,i)**2) * 0.5E+00
           mr2(j,k,i) = ur2(j,k,i) * ro2(j,k,i)
           mt2(j,k,i) = ut2(j,k,i) * ro2(j,k,i)
           mp2(j,k,i) = up2(j,k,i) * ro2(j,k,i)
           ed2(j,k,i) = pg2(j,k,i) /(gammax - 1.0E+00)
     &                +(ur2(j,k,i)**2 + ut2(j,k,i)**2
     &                 +up2(j,k,i)**2)* ro2(j,k,i) * 0.5E+00
     &                +(br2(j,k,i)**2 + bt2(j,k,i)**2
     &                 +bp2(j,k,i)**2) * 0.5E+00
           pg3(j,k,i) = pg2(j,k,i) /(gammax - 1.0E+00)
     &                +(ur2(j,k,i)**2 + ut2(j,k,i)**2
     &                 +up2(j,k,i)**2)* ro2(j,k,i) * 0.5E+00 ! MIND this is Hydro !
           br3(j,k,i) = 0.0E+00 ! MIND this is for EMF at i,j+1/2,k+1/2
           bt3(j,k,i) = 0.0E+00 !                         i+1/2,j,k+1/2
           bp3(j,k,i) = 0.0E+00 !                         i+1/2,j+1/2,k
         enddo
         enddo
       enddo
!$omp end parallel do

* counter : abnormal events
       nabnorm = 0
!$omp parallel do default(private)
!$omp&            reduction(+:nabnorm)
!$omp&            shared(ii,jj,kk)
!$omp&            shared(ro0,pg0,            br0,bt0,bp0)
!$omp&            shared(ro1,pg1,ur1,ut1,up1,br1,bt1,bp1)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(    pg3,            br3,bt3,bp3)
!$omp&            shared(ro9,pg9,ur9,ut9,up9,br9,bt9,bp9)
!$omp&            shared(mr1,mt1,mp1,ed1,mr2,mt2,mp2,ed2)
!$omp&            shared(ropark0,pgpark0,roparkb,pgparkb)
!$omp&            shared(rr)
!$omp&            shared(sain,kosa,sainb,kosab)
!$omp&            shared(sinphb,cosphb,sinph,cosph)
!$omp&            shared(kstepb)
!$omp&            shared(denf0,tmpf0)
!$omp&            shared(tempole0)
!$omp&            shared(sindth2,sindph2,cosdth2,cosdph2)
!$omp&            shared(dt,iprc)
!$omp&            shared(omega,gammax,gamma0,gammab)
!$omp&            shared(v0,r0,p0,b0,t0,rhoc,n0,ncal)
!$omp&            shared(vcrtrn,mcrtrn)
       do i = 0, ii + 1
*
         drr1 = 1.0E+00 / rr(i)
         drr2 = 1.0E+00 / rr(i)**2
         do k = 1, kk ! hopefully, this do-loop be parallelized.
         do j = 0, jj

* local-in-loop for globals  ! make clone if shared-variable is called by subroutines
           ncal2 = ncal
           iprc2 = iprc
           gammax2  = gammax
           dgammax2 = 1.0E+00 / (gammax - 1.0E+00)

*** set sinusoidal func. & differencing opperator
           rrb1 =(rr(i+1) + rr(i)) * 0.5E+00
           rrb2 =(rr(i) + rr(i-1)) * 0.5E+00
           dra0 = rrb1    - rrb2
           dra1 =(rrb1**2 - rrb2**2) * 0.5E+00
           dra2 =(rrb1**3 - rrb2**3) / 3.0E+00
           rrbc =(rrb1**2 + rrb1 * rrb2 + rrb2**2)
     &          /(rrb1 + rrb2) / 1.5E+00

*           rrbc2nd = (rrb1**2 + rrb1 * rrb2 + rrb2**2) / 3.0E+00
           rrbc2nd = rr(i)**2
           rr2ndmod = rr(i)**2
     &              /(rrb1**2 + rrb1 * rrb2 + rrb2**2) * 3.0E+00
           rr1stmod = rr(i)/rrbc

           cos0 = kosa(j)
           sin0 = sain(j)
           if (j .EQ. 0) then
             sinb1 =  sainb(0)
             sinb2 =  0.0E+00
             cosb1 =  kosab(0)
             cosb2 =  1.0E+00
             cotbc =   sainb(0) / (1.0E+00 - kosab(0))
             dth0 =          pi / float(jj) * rr(i) * 0.5E+00
             dth1 =     (1.0E+00 - kosab(0)) * rrbc
             sinbc =    (1.0E+00 - kosab(0)) / (pi/float(jj)*0.5E+00)
             cosbc =  sainb(0)               / (pi/float(jj)*0.5E+00)
           else if (j .EQ. jj) then
             sinb1 =  0.0E+00
             sinb2 =  sainb(jj-1)
             cosb1 = -1.0E+00
             cosb2 =  kosab(jj-1)
             cotbc = - sainb(jj-1) / (kosab(jj-1) + 1.0E+00)
             dth0 =          pi / float(jj) * rr(i) * 0.5E+00
             dth1 =  (kosab(jj-1) + 1.0E+00) * rrbc
             sinbc = (kosab(jj-1) + 1.0E+00) / (pi/float(jj)*0.5E+00)
             cosbc = -sainb(jj-1)            / (pi/float(jj)*0.5E+00)
           else
             sinb1 = sainb(j)
             sinb2 = sainb(j-1)
             cosb1 = kosab(j)
             cosb2 = kosab(j-1)
             cotbc = (sainb(j) - sainb(j-1)) / (kosab(j-1) - kosab(j))
             dth0 =           pi / float(jj) * rr(i)
             dth1 =   (kosab(j-1) - kosab(j)) * rrbc
             sinbc =  (kosab(j-1) - kosab(j)) / (pi/ float(jj))
             cosbc =  (sainb(j) - sainb(j-1)) / (pi/ float(jj))
           endif

           if (kk .GT. 3) then
             aa = 2.0E+00 * pi / float(kk)
           else
             aa = pi / float(jj)  ! dummy
           endif
           if ((j .EQ. 0) .OR. (j .EQ. jj)) then
             dph0 = aa * rr(i) * sinbc
           else
             dph0 = aa * rr(i) * sin0
           endif
           dph2 = aa * rrbc * sinbc
           if ((j .NE. 0) .AND. (j .NE. jj)) then
             dph0 = dph0 * float(kstepb(j))
             dph2 = dph2 * float(kstepb(j))
           endif

* here all inversed one are defined
           ddra0 = 1.0E+00 / dra0
           ddra1 = 1.0E+00 / dra1
           ddra2 = 1.0E+00 / dra2
           ddth0 = 1.0E+00 / dth0
           ddth1 = 1.0E+00 / dth1
           ddph0 = 1.0E+00 / dph0
           ddph2 = 1.0E+00 / dph2

           do idt = 1, 2
             xxrol(idt) = 0.0E+00
             xxedl(idt) = 0.0E+00
             xxmrl(idt) = 0.0E+00
             xxmtl(idt) = 0.0E+00
             xxmpl(idt) = 0.0E+00
             xxbrl(idt) = 0.0E+00
             xxbtl(idt) = 0.0E+00
             xxbpl(idt) = 0.0E+00
           enddo

           do idt = idts, idte ! --------------

* pick up variables used to calculate spatial derivatives
             do idrct = 1, 3
               idrct2 = idrct
               do l = -1, 0

                 if (idrct .EQ. 1) then
                   ro1a = ro2(j,k,i+l-1) * rr(i+l-1)**2
                   pg1a = pg2(j,k,i+l-1) * rr(i+l-1)**2
                   ed1a = ed2(j,k,i+l-1) * rr(i+l-1)**2
                   mr1a = mr2(j,k,i+l-1) * rr(i+l-1)**2
                   mt1a = mt2(j,k,i+l-1) * rr(i+l-1)**2
                   mp1a = mp2(j,k,i+l-1) * rr(i+l-1)**2
                   br1a = br2(j,k,i+l-1) * rr(i+l-1)
                   bt1a = bt2(j,k,i+l-1) * rr(i+l-1)
                   bp1a = bp2(j,k,i+l-1) * rr(i+l-1)
                   ur1a = ur2(j,k,i+l-1)
                   ut1a = ut2(j,k,i+l-1)
                   up1a = up2(j,k,i+l-1)
                   ro1b = ro2(j,k,i+l  ) * rr(i+l  )**2
                   pg1b = pg2(j,k,i+l  ) * rr(i+l  )**2
                   ed1b = ed2(j,k,i+l  ) * rr(i+l  )**2
                   mr1b = mr2(j,k,i+l  ) * rr(i+l  )**2
                   mt1b = mt2(j,k,i+l  ) * rr(i+l  )**2
                   mp1b = mp2(j,k,i+l  ) * rr(i+l  )**2
                   br1b = br2(j,k,i+l  ) * rr(i+l  )
                   bt1b = bt2(j,k,i+l  ) * rr(i+l  )
                   bp1b = bp2(j,k,i+l  ) * rr(i+l  )
                   ur1b = ur2(j,k,i+l  )
                   ut1b = ut2(j,k,i+l  )
                   up1b = up2(j,k,i+l  )
                   ro1c = ro2(j,k,i+l+1) * rr(i+l+1)**2
                   pg1c = pg2(j,k,i+l+1) * rr(i+l+1)**2
                   ed1c = ed2(j,k,i+l+1) * rr(i+l+1)**2
                   mr1c = mr2(j,k,i+l+1) * rr(i+l+1)**2
                   mt1c = mt2(j,k,i+l+1) * rr(i+l+1)**2
                   mp1c = mp2(j,k,i+l+1) * rr(i+l+1)**2
                   br1c = br2(j,k,i+l+1) * rr(i+l+1)
                   bt1c = bt2(j,k,i+l+1) * rr(i+l+1)
                   bp1c = bp2(j,k,i+l+1) * rr(i+l+1)
                   ur1c = ur2(j,k,i+l+1)
                   ut1c = ut2(j,k,i+l+1)
                   up1c = up2(j,k,i+l+1)
                   ro1d = ro2(j,k,i+l+2) * rr(i+l+2)**2
                   pg1d = pg2(j,k,i+l+2) * rr(i+l+2)**2
                   ed1d = ed2(j,k,i+l+2) * rr(i+l+2)**2
                   mr1d = mr2(j,k,i+l+2) * rr(i+l+2)**2
                   mt1d = mt2(j,k,i+l+2) * rr(i+l+2)**2
                   mp1d = mp2(j,k,i+l+2) * rr(i+l+2)**2
                   br1d = br2(j,k,i+l+2) * rr(i+l+2)
                   bt1d = bt2(j,k,i+l+2) * rr(i+l+2)
                   bp1d = bp2(j,k,i+l+2) * rr(i+l+2)
                   ur1d = ur2(j,k,i+l+2)
                   ut1d = ut2(j,k,i+l+2)
                   up1d = up2(j,k,i+l+2)
!                   if (idt .EQ. 1) then ! HD part
!                   endif
                 else if (idrct .EQ. 2) then
                   ro1a = ro2(j+l-1,k,i) * rr(i)**2
                   pg1a = pg2(j+l-1,k,i) * rr(i)**2
                   ed1a = ed2(j+l-1,k,i) * rr(i)**2
                   mr1a = mr2(j+l-1,k,i) * rr(i)**2
                   mt1a = mt2(j+l-1,k,i) * rr(i)**2
                   mp1a = mp2(j+l-1,k,i) * rr(i)**2
                   br1a = br2(j+l-1,k,i) * rr(i)
                   bt1a = bt2(j+l-1,k,i) * rr(i)
                   bp1a = bp2(j+l-1,k,i) * rr(i)
                   ur1a = ur2(j+l-1,k,i)
                   ut1a = ut2(j+l-1,k,i)
                   up1a = up2(j+l-1,k,i)
                   ro1b = ro2(j+l  ,k,i) * rr(i)**2
                   pg1b = pg2(j+l  ,k,i) * rr(i)**2
                   ed1b = ed2(j+l  ,k,i) * rr(i)**2
                   mr1b = mr2(j+l  ,k,i) * rr(i)**2
                   mt1b = mt2(j+l  ,k,i) * rr(i)**2
                   mp1b = mp2(j+l  ,k,i) * rr(i)**2
                   br1b = br2(j+l  ,k,i) * rr(i)
                   bt1b = bt2(j+l  ,k,i) * rr(i)
                   bp1b = bp2(j+l  ,k,i) * rr(i)
                   ur1b = ur2(j+l  ,k,i)
                   ut1b = ut2(j+l  ,k,i)
                   up1b = up2(j+l  ,k,i)
                   ro1c = ro2(j+l+1,k,i) * rr(i)**2
                   pg1c = pg2(j+l+1,k,i) * rr(i)**2
                   ed1c = ed2(j+l+1,k,i) * rr(i)**2
                   mr1c = mr2(j+l+1,k,i) * rr(i)**2
                   mt1c = mt2(j+l+1,k,i) * rr(i)**2
                   mp1c = mp2(j+l+1,k,i) * rr(i)**2
                   br1c = br2(j+l+1,k,i) * rr(i)
                   bt1c = bt2(j+l+1,k,i) * rr(i)
                   bp1c = bp2(j+l+1,k,i) * rr(i)
                   ur1c = ur2(j+l+1,k,i)
                   ut1c = ut2(j+l+1,k,i)
                   up1c = up2(j+l+1,k,i)
                   ro1d = ro2(j+l+2,k,i) * rr(i)**2
                   pg1d = pg2(j+l+2,k,i) * rr(i)**2
                   ed1d = ed2(j+l+2,k,i) * rr(i)**2
                   mr1d = mr2(j+l+2,k,i) * rr(i)**2
                   mt1d = mt2(j+l+2,k,i) * rr(i)**2
                   mp1d = mp2(j+l+2,k,i) * rr(i)**2
                   br1d = br2(j+l+2,k,i) * rr(i)
                   bt1d = bt2(j+l+2,k,i) * rr(i)
                   bp1d = bp2(j+l+2,k,i) * rr(i)
                   ur1d = ur2(j+l+2,k,i)
                   ut1d = ut2(j+l+2,k,i)
                   up1d = up2(j+l+2,k,i)
!                   if (idt .EQ. 1) then ! HD part
!                   endif
                 else ! if (idrct .EQ. 3)
                   kl1 = k + (l - 1) * kstepb(j) + kk * 2
                   kl1 = mod(kl1,kk)
                   kl2 = k +  l      * kstepb(j) + kk * 2
                   kl2 = mod(kl2,kk)
                   kl3 = k + (l + 1) * kstepb(j) + kk * 2
                   kl3 = mod(kl3,kk)
                   kl4 = k + (l + 2) * kstepb(j) + kk * 2
                   kl4 = mod(kl4,kk)

                   ro1a = ro2(j,kl1,i) * rr(i)**2
                   pg1a = pg2(j,kl1,i) * rr(i)**2
                   ed1a = ed2(j,kl1,i) * rr(i)**2
                   mr1a = mr2(j,kl1,i) * rr(i)**2
                   mt1a = mt2(j,kl1,i) * rr(i)**2
                   mp1a = mp2(j,kl1,i) * rr(i)**2
                   br1a = br2(j,kl1,i) * rr(i)
                   bt1a = bt2(j,kl1,i) * rr(i)
                   bp1a = bp2(j,kl1,i) * rr(i)
                   ur1a = ur2(j,kl1,i)
                   ut1a = ut2(j,kl1,i)
                   up1a = up2(j,kl1,i)
                   ro1b = ro2(j,kl2,i) * rr(i)**2
                   pg1b = pg2(j,kl2,i) * rr(i)**2
                   ed1b = ed2(j,kl2,i) * rr(i)**2
                   mr1b = mr2(j,kl2,i) * rr(i)**2
                   mt1b = mt2(j,kl2,i) * rr(i)**2
                   mp1b = mp2(j,kl2,i) * rr(i)**2
                   br1b = br2(j,kl2,i) * rr(i)
                   bt1b = bt2(j,kl2,i) * rr(i)
                   bp1b = bp2(j,kl2,i) * rr(i)
                   ur1b = ur2(j,kl2,i)
                   ut1b = ut2(j,kl2,i)
                   up1b = up2(j,kl2,i)
                   ro1c = ro2(j,kl3,i) * rr(i)**2
                   pg1c = pg2(j,kl3,i) * rr(i)**2
                   ed1c = ed2(j,kl3,i) * rr(i)**2
                   mr1c = mr2(j,kl3,i) * rr(i)**2
                   mt1c = mt2(j,kl3,i) * rr(i)**2
                   mp1c = mp2(j,kl3,i) * rr(i)**2
                   br1c = br2(j,kl3,i) * rr(i)
                   bt1c = bt2(j,kl3,i) * rr(i)
                   bp1c = bp2(j,kl3,i) * rr(i)
                   ur1c = ur2(j,kl3,i)
                   ut1c = ut2(j,kl3,i)
                   up1c = up2(j,kl3,i)
                   ro1d = ro2(j,kl4,i) * rr(i)**2
                   pg1d = pg2(j,kl4,i) * rr(i)**2
                   ed1d = ed2(j,kl4,i) * rr(i)**2
                   mr1d = mr2(j,kl4,i) * rr(i)**2
                   mt1d = mt2(j,kl4,i) * rr(i)**2
                   mp1d = mp2(j,kl4,i) * rr(i)**2
                   br1d = br2(j,kl4,i) * rr(i)
                   bt1d = bt2(j,kl4,i) * rr(i)
                   bp1d = bp2(j,kl4,i) * rr(i)
                   ur1d = ur2(j,kl4,i)
                   ut1d = ut2(j,kl4,i)
                   up1d = up2(j,kl4,i)
!                   if (idt .EQ. 1) then ! HD part
!                   endif
                 endif ! end if (idrct .EQ ......)

* restore local array
                 rob(l+1,2,idrct) = ro1b
                 rob(l+1,1,idrct) = ro1c
                 pgb(l+1,2,idrct) = pg1b
                 pgb(l+1,1,idrct) = pg1c
                 edb(l+1,2,idrct) = ed1b
                 edb(l+1,1,idrct) = ed1c
                 mrb(l+1,2,idrct) = mr1b
                 mrb(l+1,1,idrct) = mr1c
                 mtb(l+1,2,idrct) = mt1b
                 mtb(l+1,1,idrct) = mt1c
                 mpb(l+1,2,idrct) = mp1b
                 mpb(l+1,1,idrct) = mp1c
                 brb(l+1,2,idrct) = br1b
                 brb(l+1,1,idrct) = br1c
                 btb(l+1,2,idrct) = bt1b
                 btb(l+1,1,idrct) = bt1c
                 bpb(l+1,2,idrct) = bp1b
                 bpb(l+1,1,idrct) = bp1c
                 urb(l+1,2,idrct) = ur1b
                 urb(l+1,1,idrct) = ur1c
                 utb(l+1,2,idrct) = ut1b
                 utb(l+1,1,idrct) = ut1c
                 upb(l+1,2,idrct) = up1b
                 upb(l+1,1,idrct) = up1c

!                 lradprk2 = lradpark

* MUSCL
                 if (lmuscl) then

*                 if (iprc .EQ. 2) then ! old test... never used for scientific purpose

                  if ((idrct .LE. 2) .OR.
     &                ((idrct .EQ. 3) .AND. (kk .GE. 3))) then ! skip phi-direction if 2D

* local difference for normalized variables
!                   if ((lradprk2) .AND. (idrct .EQ. 1)) then
                   if ((lradpark) .AND. (idrct .EQ. 1)) then
!                     if (lradlog) then
                       rod1 = log(ro1b) - log(ro1a)
                       rod2 = log(ro1c) - log(ro1b)
                       rod3 = log(ro1d) - log(ro1c)
                       pgd1 = log(pg1b) - log(pg1a)
                       pgd2 = log(pg1c) - log(pg1b)
                       pgd3 = log(pg1d) - log(pg1c)
!                     else
!                     endif
                   else
                     rod1 = ro1b - ro1a
                     rod2 = ro1c - ro1b
                     rod3 = ro1d - ro1c
                     pgd1 = pg1b - pg1a
                     pgd2 = pg1c - pg1b
                     pgd3 = pg1d - pg1c
                   endif
                   edd1 = ed1b - ed1a
                   edd2 = ed1c - ed1b
                   edd3 = ed1d - ed1c
                   mrd1 = mr1b - mr1a
                   mrd2 = mr1c - mr1b
                   mrd3 = mr1d - mr1c
                   mtd1 = mt1b - mt1a
                   mtd2 = mt1c - mt1b
                   mtd3 = mt1d - mt1c
                   mpd1 = mp1b - mp1a
                   mpd2 = mp1c - mp1b
                   mpd3 = mp1d - mp1c
                   brd1 = br1b - br1a
                   brd2 = br1c - br1b
                   brd3 = br1d - br1c
                   btd1 = bt1b - bt1a
                   btd2 = bt1c - bt1b
                   btd3 = bt1d - bt1c
                   bpd1 = bp1b - bp1a
                   bpd2 = bp1c - bp1b
                   bpd3 = bp1d - bp1c
                   urd1 = ur1b - ur1a
                   urd2 = ur1c - ur1b
                   urd3 = ur1d - ur1c
                   utd1 = ut1b - ut1a
                   utd2 = ut1c - ut1b
                   utd3 = ut1d - ut1c
                   upd1 = up1b - up1a
                   upd2 = up1c - up1b
                   upd3 = up1d - up1c

* some exception at boundary and adjustment for non-uniform grid..
                   if (idrct .EQ. 1) then
                     if (i+l .EQ. -1) then
                       rod1 = rod2 * 2.0E+00 - rod3 ! 2nd order extrapolation for i=-3/2
                       pgd1 = pgd2 * 2.0E+00 - pgd3
                       edd1 = edd2 * 2.0E+00 - edd3
                       urd1 = urd2 * 2.0E+00 - urd3
                       utd1 = utd2 * 2.0E+00 - utd3
                       upd1 = upd2 * 2.0E+00 - upd3
                       mrd1 = mrd2 * 2.0E+00 - mrd3
                       mtd1 = mtd2 * 2.0E+00 - mtd3
                       mpd1 = mpd2 * 2.0E+00 - mpd3
                       brd1 = brd2 * 2.0E+00 - brd3
                       btd1 = btd2 * 2.0E+00 - btd3
                       bpd1 = bpd2 * 2.0E+00 - bpd3
                     else if (i+l .EQ. ii+1) then ! 0th
                       rod3 = rod1
                       edd3 = edd1
                       mrd3 = mrd1
                       mtd3 = mtd1
                       mpd3 = mpd1
                       brd3 = brd1
                       btd3 = btd1
                       bpd3 = bpd1
                       pgd3 = pgd1
                       urd3 = urd1
                       utd3 = utd1
                       upd3 = upd1
                       rod2 = rod1
                       edd2 = edd1
                       mrd2 = mrd1
                       mtd2 = mtd1
                       mpd2 = mpd1
                       brd2 = brd1
                       btd2 = btd1
                       bpd2 = bpd1
                       pgd2 = pgd1
                       urd2 = urd1
                       utd2 = utd1
                       upd2 = upd1
                     endif
                   else if (idrct .EQ. 2) then
                     if (iextpole .EQ. 0) then
                       if (j+l .EQ. -1) then
                         mtd1 = 0.0E+00
                         mpd1 = 0.0E+00
                         btd1 = 0.0E+00
                         bpd1 = 0.0E+00
                         utd1 = 0.0E+00
                         upd1 = 0.0E+00
                         mtd2 = 0.0E+00
                         mpd2 = 0.0E+00
                         btd2 = 0.0E+00
                         bpd2 = 0.0E+00
                         utd2 = 0.0E+00
                         upd2 = 0.0E+00
                       else if (j+l .EQ. jj) then
                         mtd2 = 0.0E+00
                         mpd2 = 0.0E+00
                         btd2 = 0.0E+00
                         bpd2 = 0.0E+00
                         utd2 = 0.0E+00
                         upd2 = 0.0E+00
                         mtd3 = 0.0E+00
                         mpd3 = 0.0E+00
                         btd3 = 0.0E+00
                         bpd3 = 0.0E+00
                         utd3 = 0.0E+00
                         upd3 = 0.0E+00
                       endif
                     else if (iextpole .EQ. 1) then
                       if (j+l .EQ. -1) then
                         rod1 = rod3
                         edd1 = edd3
                         mrd1 = mrd3
                         mtd1 = mtd3
                         mpd1 = mpd3
                         brd1 = brd3
                         btd1 = btd3
                         bpd1 = bpd3
                         pgd1 = pgd3
                         urd1 = urd3
                         utd1 = utd3
                         upd1 = upd3
                         rod2 = rod3
                         edd2 = edd3
                         mrd2 = mrd3
                         mtd2 = mtd3
                         mpd2 = mpd3
                         brd2 = brd3
                         btd2 = btd3
                         bpd2 = bpd3
                         pgd2 = pgd3
                         urd2 = urd3
                         utd2 = utd3
                         upd2 = upd3
                       else if (j+l .EQ. jj) then
                         rod2 = rod1
                         edd2 = edd1
                         mrd2 = mrd1
                         mtd2 = mtd1
                         mpd2 = mpd1
                         brd2 = brd1
                         btd2 = btd1
                         bpd2 = bpd1
                         pgd2 = pgd1
                         urd2 = urd1
                         utd2 = utd1
                         upd2 = upd1
                         rod3 = rod1
                         edd3 = edd1
                         mrd3 = mrd1
                         mtd3 = mtd1
                         mpd3 = mpd1
                         brd3 = brd1
                         btd3 = btd1
                         bpd3 = bpd1
                         pgd3 = pgd1
                         urd3 = urd1
                         utd3 = utd1
                         upd3 = upd1
                       endif
                     else if (iextpole .EQ. 2) then
                       if (j+l .EQ. -1) then
                         rod1 = 0.0E+00
                         edd1 = 0.0E+00
                         mrd1 = 0.0E+00
                         mtd1 = 0.0E+00
                         mpd1 = 0.0E+00
                         brd1 = 0.0E+00
                         btd1 = 0.0E+00
                         bpd1 = 0.0E+00
                         pgd1 = 0.0E+00
                         urd1 = 0.0E+00
                         utd1 = 0.0E+00
                         upd1 = 0.0E+00
                         rod2 = 0.0E+00
                         edd2 = 0.0E+00
                         mrd2 = 0.0E+00
                         mtd2 = 0.0E+00
                         mpd2 = 0.0E+00
                         brd2 = 0.0E+00
                         btd2 = 0.0E+00
                         bpd2 = 0.0E+00
                         pgd2 = 0.0E+00
                         urd2 = 0.0E+00
                         utd2 = 0.0E+00
                         upd2 = 0.0E+00
                       else if (j+l .EQ. jj) then
                         rod2 = 0.0E+00
                         edd2 = 0.0E+00
                         mrd2 = 0.0E+00
                         mtd2 = 0.0E+00
                         mpd2 = 0.0E+00
                         brd2 = 0.0E+00
                         btd2 = 0.0E+00
                         bpd2 = 0.0E+00
                         pgd2 = 0.0E+00
                         urd2 = 0.0E+00
                         utd2 = 0.0E+00
                         upd2 = 0.0E+00
                         rod3 = 0.0E+00
                         edd3 = 0.0E+00
                         mrd3 = 0.0E+00
                         mtd3 = 0.0E+00
                         mpd3 = 0.0E+00
                         brd3 = 0.0E+00
                         btd3 = 0.0E+00
                         bpd3 = 0.0E+00
                         pgd3 = 0.0E+00
                         urd3 = 0.0E+00
                         utd3 = 0.0E+00
                         upd3 = 0.0E+00
                       endif
                     endif ! end if (iextpole .EQ.....)

                   endif ! end if (idirect .EQ.....)

                   call avehalf(gammax2,! idrct2,
     &                      ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                      ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                      roa,ura,uta,upa,bra,bta,bpa,pga)
                   call vrtp2num(idrct2,
     &                      ura,uta,upa,bra,bta,bpa,
     &                      v1,v2,v3,b1,b2,b3)
                   if (idrct .EQ. 1) then
                     dx1 = rr(i+l  ) - rr(i+l-1)
                     dx2 = rr(i+l+1) - rr(i+l  )
                     dx3 = rr(i+l+2) - rr(i+l+1)
                   else
                     dx1 = 1.0E+00 ! MIND this assume uniform d(th) & d(ph)
                     dx2 = 1.0E+00
                     dx3 = 1.0E+00
                   endif

                     sqrpg = (pgb(l+1,1,idrct) + pgb(l+1,2,idrct))
     &                     * 0.5E+00
                     sqrpg = sqrt(sqrpg)

!                     if (lradprk2) then
                       if (idrct .EQ. 1) then
!                         if (lradlog) then
                           aa = log(rob(l+1,1,idrct))
     &                        - mscl3rd(rod3,rod2,dx3,dx2)
                           bb = log(rob(l+1,2,idrct))
     &                        + mscl3rd(rod1,rod2,dx1,dx2)
                           rob(l+1,1,idrct) = exp(aa)
                           rob(l+1,2,idrct) = exp(bb)
!                         else
!                         endif
                       else
                         aa=rob(l+1,1,idrct)-mscl3rd(rod3,rod2,dx3,dx2)
                         bb=rob(l+1,2,idrct)+mscl3rd(rod1,rod2,dx1,dx2)
                         rob(l+1,1,idrct) = aa
                         rob(l+1,2,idrct) = bb
                       endif
!                     else
!                     endif

!                     if (lradprk2) then
                       if (idrct .EQ. 1) then
!                         if (lradlog) then
                           aa =log(pgb(l+1,1,idrct))
     &                        -mscl3rd(pgd3,pgd2,dx3,dx2)
                           bb =log(pgb(l+1,2,idrct))
     &                        +mscl3rd(pgd1,pgd2,dx1,dx2)
                           pgb(l+1,1,idrct) = exp(aa)
                           pgb(l+1,2,idrct) = exp(bb)
!                         else
!                         endif
                       else
                         aa=pgb(l+1,1,idrct)-mscl3rd(pgd3,pgd2,dx3,dx2)
                         bb=pgb(l+1,2,idrct)+mscl3rd(pgd1,pgd2,dx1,dx2)
                         pgb(l+1,1,idrct) = aa
                         pgb(l+1,2,idrct) = bb
                       endif
!                     else
!                     endif

!                     if (lmmscl2) then
                       mrb(l+1,1,idrct)=mrb(l+1,1,idrct)
     &                                 -mscl3rd(mrd3,mrd2,dx3,dx2)
                       mrb(l+1,2,idrct)=mrb(l+1,2,idrct)
     &                                 +mscl3rd(mrd1,mrd2,dx1,dx2)
                       mtb(l+1,1,idrct)=mtb(l+1,1,idrct)
     &                                 -mscl3rd(mtd3,mtd2,dx3,dx2)
                       mtb(l+1,2,idrct)=mtb(l+1,2,idrct)
     &                                 +mscl3rd(mtd1,mtd2,dx1,dx2)
                       mpb(l+1,1,idrct)=mpb(l+1,1,idrct)
     &                                 -mscl3rd(mpd3,mpd2,dx3,dx2)
                       mpb(l+1,2,idrct)=mpb(l+1,2,idrct)
     &                                 +mscl3rd(mpd1,mpd2,dx1,dx2)
                       urb(l+1,1,idrct)=mrb(l+1,1,idrct)
     &                                 /rob(l+1,1,idrct)
                       urb(l+1,2,idrct)=mrb(l+1,2,idrct)
     &                                 /rob(l+1,2,idrct)
                       utb(l+1,1,idrct)=mtb(l+1,1,idrct)
     &                                 /rob(l+1,1,idrct)
                       utb(l+1,2,idrct)=mtb(l+1,2,idrct)
     &                                 /rob(l+1,2,idrct)
                       upb(l+1,1,idrct)=mpb(l+1,1,idrct)
     &                                 /rob(l+1,1,idrct)
                       upb(l+1,2,idrct)=mpb(l+1,2,idrct)
     &                                 /rob(l+1,2,idrct)
!                     else
!                     endif

                     if (idrct .EQ. 1)  then ! this part is different from origial q30a1f

!                       if (ib1intp .EQ. 1) then
!                       else if (ib1intp .EQ. 2) then
                         drestmg(1)=nonmscl3(brd3,brd2,brd1,dx3,dx2,dx1)
                         drestmg(2)=nonmscl3(brd1,brd2,brd3,dx1,dx2,dx3)
!                       else if (ib1intp .EQ. 3) then
!                       else if (ib1intp .EQ. 4) then
!                       else if (ib1intp .EQ. 5) then
!                       else
!                       endif
                       brb(l+1,1,idrct) = brb(l+1,1,idrct) - drestmg(1)
                       brb(l+1,2,idrct) = brb(l+1,2,idrct) + drestmg(2)

                     else
                       dd = brb(l+1,1,idrct)-mscl3rd(brd3,brd2,dx3,dx2)
                       ee = brb(l+1,2,idrct)+mscl3rd(brd1,brd2,dx1,dx2)
                       ff = brb(l+1,1,idrct)
     &                    - nonmscl3(brd3,brd2,brd1,dx3,dx2,dx1)
                       gg = brb(l+1,2,idrct)
     &                    + nonmscl3(brd1,brd2,brd3,dx1,dx2,dx3)
                       hh= abs(brd1 - 2.0E+00 * brd2 + brd3)
     &                   /(abs(brd2) + gradsmll * sqrpg) * 0.5E+00
                       hh=min(1.0E+00,
     &                        max(0.0E+00,gradfact*(hh-gradlow)))
                       brb(l+1,1,idrct) = hh * dd + (1.0E+00 - hh) * ff
                       brb(l+1,2,idrct) = hh * ee + (1.0E+00 - hh) * gg
                     endif

                     if (idrct .EQ. 2) then

!                       if (ib1intp .EQ. 1) then
!                       else if (ib1intp .EQ. 2) then
                         drestmg(1)=nonmscl3(btd3,btd2,btd1,dx3,dx2,dx1)
                         drestmg(2)=nonmscl3(btd1,btd2,btd3,dx1,dx2,dx3)
!                       else if (ib1intp .EQ. 3) then
!                       else if (ib1intp .EQ. 4) then
!                       else if (ib1intp .EQ. 5) then
!                       else
!                       endif
                       btb(l+1,1,idrct) = btb(l+1,1,idrct) - drestmg(1)
                       btb(l+1,2,idrct) = btb(l+1,2,idrct) + drestmg(2)

                     else
                       dd = btb(l+1,1,idrct)-mscl3rd(btd3,btd2,dx3,dx2)
                       ee = btb(l+1,2,idrct)+mscl3rd(btd1,btd2,dx1,dx2)
                       ff = btb(l+1,1,idrct)
     &                    - nonmscl3(btd3,btd2,btd1,dx3,dx2,dx1)
                       gg = btb(l+1,2,idrct)
     &                    + nonmscl3(btd1,btd2,btd3,dx1,dx2,dx3)
                       hh= abs(btd1 - 2.0E+00 * btd2 + btd3)
     &                   /(abs(btd2) + gradsmll * sqrpg) * 0.5E+00
                       hh=min(1.0E+00,
     &                        max(0.0E+00,gradfact*(hh-gradlow)))
                       btb(l+1,1,idrct) = hh * dd + (1.0E+00 - hh) * ff
                       btb(l+1,2,idrct) = hh * ee + (1.0E+00 - hh) * gg
                     endif

                     if (idrct .EQ. 3) then

!                       if (ib1intp .EQ. 1) then
!                       else if (ib1intp .EQ. 2) then
                         drestmg(1)=nonmscl3(bpd3,bpd2,bpd1,dx3,dx2,dx1)
                         drestmg(2)=nonmscl3(bpd1,bpd2,bpd3,dx1,dx2,dx3)
!                       else if (ib1intp .EQ. 3) then
!                       else if (ib1intp .EQ. 4) then
!                       else if (ib1intp .EQ. 5) then
!                       else
!                       endif
                       bpb(l+1,1,idrct) = bpb(l+1,1,idrct) - drestmg(1)
                       bpb(l+1,2,idrct) = bpb(l+1,2,idrct) + drestmg(2)

                     else
                       dd = bpb(l+1,1,idrct)-mscl3rd(bpd3,bpd2,dx3,dx2)
                       ee = bpb(l+1,2,idrct)+mscl3rd(bpd1,bpd2,dx1,dx2)
                       ff = bpb(l+1,1,idrct)
     &                    - nonmscl3(bpd3,bpd2,bpd1,dx3,dx2,dx1)
                       gg = bpb(l+1,2,idrct)
     &                    + nonmscl3(bpd1,bpd2,bpd3,dx1,dx2,dx3)
                       hh=abs(bpd1 - 2.0E+00 * bpd2 + bpd3)
     &                   /(abs(bpd2) + gradsmll * sqrpg) * 0.5E+00
                       hh=min(1.0E+00,
     &                        max(0.0E+00,gradfact*(hh-gradlow)))
                       bpb(l+1,1,idrct) = hh * dd + (1.0E+00 - hh) * ff
                       bpb(l+1,2,idrct) = hh * ee + (1.0E+00 - hh) * gg
                     endif

                     edb(l+1,1,idrct)
     &                  = pgb(l+1,1,idrct) * dgammax2
     &                  +(urb(l+1,1,idrct)**2 + utb(l+1,1,idrct)**2
     &                   +upb(l+1,1,idrct)**2)
     &                   *rob(l+1,1,idrct) * 0.5E+00
     &                  +(brb(l+1,1,idrct)**2 + btb(l+1,1,idrct)**2
     &                   +bpb(l+1,1,idrct)**2) * 0.5E+00
                     edb(l+1,2,idrct)
     &                  = pgb(l+1,2,idrct) * dgammax2
     &                  +(urb(l+1,2,idrct)**2 + utb(l+1,2,idrct)**2
     &                   +upb(l+1,2,idrct)**2)
     &                   *rob(l+1,2,idrct) * 0.5E+00
     &                  +(brb(l+1,2,idrct)**2 + btb(l+1,2,idrct)**2
     &                   +bpb(l+1,2,idrct)**2) * 0.5E+00

                  endif ! if (idrct ......)

                 endif ! if (lmuscl)

               enddo ! end of l = -1, 0

* adjustment at i = 0, if chosen
               if (lvxbzero) then
                 if ((i .EQ. 0) .AND. (idrct .NE. 1)) then
                   do l = -1, 0
                     bb = brb(l+1,1,idrct)**2 + btb(l+1,1,idrct)**2
     &                  + bpb(l+1,1,idrct)**2
                     bb = sqrt(bb + 1.0E-20)
                     aa  = abs(brb(l+1,1,idrct)) / (bb + 1.0E-05)
                     if (aa .GT. brabsb) then
                       aa = mrb(l+1,1,idrct) / brb(l+1,1,idrct)
                       mtb(l+1,1,idrct) = aa * btb(l+1,1,idrct)
                       mpb(l+1,1,idrct) = aa * bpb(l+1,1,idrct)
                     else
                       mrb(l+1,1,idrct) = 0.0E+00 ! stagnant
                       mtb(l+1,1,idrct) = 0.0E+00
                       mpb(l+1,1,idrct) = 0.0E+00
                     endif

                     bb = brb(l+1,2,idrct)**2 + btb(l+1,2,idrct)**2
     &                  + bpb(l+1,2,idrct)**2
                     bb = sqrt(bb + 1.0E-20)
                     aa  = abs(brb(l+1,2,idrct)) / (bb + 1.0E-05)
                     if (aa .GT. brabsb) then
                       aa = mrb(l+1,2,idrct) / brb(l+1,2,idrct)
                       mtb(l+1,2,idrct) = aa * btb(l+1,2,idrct)
                       mpb(l+1,2,idrct) = aa * bpb(l+1,2,idrct)
                     else ! stagnant
                       mrb(l+1,2,idrct) = 0.0E+00
                       mtb(l+1,2,idrct) = 0.0E+00
                       mpb(l+1,2,idrct) = 0.0E+00
                     endif

                     urb(l+1,1,idrct)=mrb(l+1,1,idrct)/rob(l+1,1,idrct)
                     utb(l+1,1,idrct)=mtb(l+1,1,idrct)/rob(l+1,1,idrct)
                     upb(l+1,1,idrct)=mpb(l+1,1,idrct)/rob(l+1,1,idrct)
                     urb(l+1,2,idrct)=mrb(l+1,2,idrct)/rob(l+1,2,idrct)
                     utb(l+1,2,idrct)=mtb(l+1,2,idrct)/rob(l+1,2,idrct)
                     upb(l+1,2,idrct)=mpb(l+1,2,idrct)/rob(l+1,2,idrct)
                   enddo
                 endif
               endif

             enddo ! end of idrct = 1, 3

* calculate aux. term
             idrct2 = 1
             l = 0
             ro1c = rob(l+1,1,idrct2)
             ro1b = rob(l+1,2,idrct2)
             mr1c = mrb(l+1,1,idrct2)
             mr1b = mrb(l+1,2,idrct2)
             mt1c = mtb(l+1,1,idrct2)
             mt1b = mtb(l+1,2,idrct2)
             mp1c = mpb(l+1,1,idrct2)
             mp1b = mpb(l+1,2,idrct2)
             br1c = brb(l+1,1,idrct2)
             br1b = brb(l+1,2,idrct2)
             bt1c = btb(l+1,1,idrct2)
             bt1b = btb(l+1,2,idrct2)
             bp1c = bpb(l+1,1,idrct2)
             bp1b = bpb(l+1,2,idrct2)
             ed1c = edb(l+1,1,idrct2)
             ed1b = edb(l+1,2,idrct2)
             call avehalf(gammax2, ! idrct2,
     &                    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                    roa,ura,uta,upa,bra,bta,bpa,pga)
             ro1r1 = roa
             pg1r1 = pga
             ur1r1 = ura
             ut1r1 = uta
             up1r1 = upa
             br1r1 = bra
             bt1r1 = bta
             bp1r1 = bpa
             mr1r1 = ura * roa
             mt1r1 = uta * roa
             mp1r1 = upa * roa

             idrct2 = 1
             l = -1
             ro1c = rob(l+1,1,idrct2)
             ro1b = rob(l+1,2,idrct2)
             mr1c = mrb(l+1,1,idrct2)
             mr1b = mrb(l+1,2,idrct2)
             mt1c = mtb(l+1,1,idrct2)
             mt1b = mtb(l+1,2,idrct2)
             mp1c = mpb(l+1,1,idrct2)
             mp1b = mpb(l+1,2,idrct2)
             br1c = brb(l+1,1,idrct2)
             br1b = brb(l+1,2,idrct2)
             bt1c = btb(l+1,1,idrct2)
             bt1b = btb(l+1,2,idrct2)
             bp1c = bpb(l+1,1,idrct2)
             bp1b = bpb(l+1,2,idrct2)
             ed1c = edb(l+1,1,idrct2)
             ed1b = edb(l+1,2,idrct2)
             call avehalf(gammax2, ! idrct2,
     &                    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                    roa,ura,uta,upa,bra,bta,bpa,pga)
             ro1r2 = roa
             pg1r2 = pga
             ur1r2 = ura
             ut1r2 = uta
             up1r2 = upa
             br1r2 = bra
             bt1r2 = bta
             bp1r2 = bpa
             mr1r2 = ura * roa
             mt1r2 = uta * roa
             mp1r2 = upa * roa

             idrct2 = 2
             l = 0
             ro1c = rob(l+1,1,idrct2)
             ro1b = rob(l+1,2,idrct2)
             mr1c = mrb(l+1,1,idrct2)
             mr1b = mrb(l+1,2,idrct2)
             mt1c = mtb(l+1,1,idrct2)
             mt1b = mtb(l+1,2,idrct2)
             mp1c = mpb(l+1,1,idrct2)
             mp1b = mpb(l+1,2,idrct2)
             br1c = brb(l+1,1,idrct2)
             br1b = brb(l+1,2,idrct2)
             bt1c = btb(l+1,1,idrct2)
             bt1b = btb(l+1,2,idrct2)
             bp1c = bpb(l+1,1,idrct2)
             bp1b = bpb(l+1,2,idrct2)
             ed1c = edb(l+1,1,idrct2)
             ed1b = edb(l+1,2,idrct2)
             call avehalf(gammax2, ! idrct2,
     &                    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                    roa,ura,uta,upa,bra,bta,bpa,pga)
             ro1t1 = roa
             pg1t1 = pga
             ur1t1 = ura
             ut1t1 = uta
             up1t1 = upa
             br1t1 = bra
             bt1t1 = bta
             bp1t1 = bpa
             mr1t1 = ura * roa
             mt1t1 = uta * roa
             mp1t1 = upa * roa

             idrct2 = 2
             l = -1
             ro1c = rob(l+1,1,idrct2)
             ro1b = rob(l+1,2,idrct2)
             mr1c = mrb(l+1,1,idrct2)
             mr1b = mrb(l+1,2,idrct2)
             mt1c = mtb(l+1,1,idrct2)
             mt1b = mtb(l+1,2,idrct2)
             mp1c = mpb(l+1,1,idrct2)
             mp1b = mpb(l+1,2,idrct2)
             br1c = brb(l+1,1,idrct2)
             br1b = brb(l+1,2,idrct2)
             bt1c = btb(l+1,1,idrct2)
             bt1b = btb(l+1,2,idrct2)
             bp1c = bpb(l+1,1,idrct2)
             bp1b = bpb(l+1,2,idrct2)
             ed1c = edb(l+1,1,idrct2)
             ed1b = edb(l+1,2,idrct2)
             call avehalf(gammax2, ! idrct2,
     &                    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                    roa,ura,uta,upa,bra,bta,bpa,pga)
             ro1t2 = roa
             pg1t2 = pga
             ur1t2 = ura
             ut1t2 = uta
             up1t2 = upa
             br1t2 = bra
             bt1t2 = bta
             bp1t2 = bpa
             mr1t2 = ura * roa
             mt1t2 = uta * roa
             mp1t2 = upa * roa

             idrct2 = 3
             l = 0
             ro1c = rob(l+1,1,idrct2)
             ro1b = rob(l+1,2,idrct2)
             mr1c = mrb(l+1,1,idrct2)
             mr1b = mrb(l+1,2,idrct2)
             mt1c = mtb(l+1,1,idrct2)
             mt1b = mtb(l+1,2,idrct2)
             mp1c = mpb(l+1,1,idrct2)
             mp1b = mpb(l+1,2,idrct2)
             br1c = brb(l+1,1,idrct2)
             br1b = brb(l+1,2,idrct2)
             bt1c = btb(l+1,1,idrct2)
             bt1b = btb(l+1,2,idrct2)
             bp1c = bpb(l+1,1,idrct2)
             bp1b = bpb(l+1,2,idrct2)
             ed1c = edb(l+1,1,idrct2)
             ed1b = edb(l+1,2,idrct2)
             call avehalf(gammax2, ! idrct2,
     &                    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                    roa,ura,uta,upa,bra,bta,bpa,pga)
             ro1p1 = roa
             pg1p1 = pga
             ur1p1 = ura
             ut1p1 = uta
             up1p1 = upa
             br1p1 = bra
             bt1p1 = bta
             bp1p1 = bpa
             mr1p1 = ura * roa
             mt1p1 = uta * roa
             mp1p1 = upa * roa

             idrct2 = 3
             l = -1
             ro1c = rob(l+1,1,idrct2)
             ro1b = rob(l+1,2,idrct2)
             mr1c = mrb(l+1,1,idrct2)
             mr1b = mrb(l+1,2,idrct2)
             mt1c = mtb(l+1,1,idrct2)
             mt1b = mtb(l+1,2,idrct2)
             mp1c = mpb(l+1,1,idrct2)
             mp1b = mpb(l+1,2,idrct2)
             br1c = brb(l+1,1,idrct2)
             br1b = brb(l+1,2,idrct2)
             bt1c = btb(l+1,1,idrct2)
             bt1b = btb(l+1,2,idrct2)
             bp1c = bpb(l+1,1,idrct2)
             bp1b = bpb(l+1,2,idrct2)
             ed1c = edb(l+1,1,idrct2)
             ed1b = edb(l+1,2,idrct2)
             call avehalf(gammax2, ! idrct2,
     &                    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                    roa,ura,uta,upa,bra,bta,bpa,pga)
             ro1p2 = roa
             pg1p2 = pga
             ur1p2 = ura
             ut1p2 = uta
             up1p2 = upa
             br1p2 = bra
             bt1p2 = bta
             bp1p2 = bpa
             mr1p2 = ura * roa
             mt1p2 = uta * roa
             mp1p2 = upa * roa

             ro1cc = (ro1r1+ro1r2+ro1t1+ro1t2+ro1p1+ro1p2)/6.0E+00
             pg1cc = (pg1r1+pg1r2+pg1t1+pg1t2+pg1p1+pg1p2)/6.0E+00
             mr1cc = (mr1r1+mr1r2+mr1t1+mr1t2+mr1p1+mr1p2)/6.0E+00
             mt1cc = (mt1r1+mt1r2+mt1t1+mt1t2+mt1p1+mt1p2)/6.0E+00
             mp1cc = (mp1r1+mp1r2+mp1t1+mp1t2+mp1p1+mp1p2)/6.0E+00
             ur1cc = (ur1r1+ur1r2+ur1t1+ur1t2+ur1p1+ur1p2)/6.0E+00
             ut1cc = (ut1r1+ut1r2+ut1t1+ut1t2+ut1p1+ut1p2)/6.0E+00
             up1cc = (up1r1+up1r2+up1t1+up1t2+up1p1+up1p2)/6.0E+00
             br1cc = (br1r1+br1r2+br1t1+br1t2+br1p1+br1p2)/6.0E+00
             bt1cc = (bt1r1+bt1r2+bt1t1+bt1t2+bt1p1+bt1p2)/6.0E+00
             bp1cc = (bp1r1+bp1r2+bp1t1+bp1t2+bp1p1+bp1p2)/6.0E+00

* make sure vxb=0 at i=0
             if (lvxbzerc) then
               if (i .EQ. 0) then
                 bb = br1cc**2 + bt1cc**2 + bp1cc**2
                 bb = sqrt(bb + 1.0E-20)
                 if (bb .GT. 1.0E-05) then
                   aa  = abs(br1cc) / (bb + 1.0E-05)
                   if (aa .GT. brabsb) then
                     aa = mr1cc / br1cc
                     mt1cc = aa * bt1cc
                     mp1cc = aa * bp1cc
                   else
                     mr1cc = 0.0E+00 ! stagnant
                     mt1cc = 0.0E+00
                     mp1cc = 0.0E+00
                   endif
                 else ! if B ~ 0
                   mt1cc = 0.0E+00
                   mp1cc = 0.0E+00
                 endif
               endif
             endif

*TVD modification term
             do l = 0, 1
             do m = 1, 7
               tvdcr(m,l) = 0.0E+00
               tvdct(m,l) = 0.0E+00
               tvdcp(m,l) = 0.0E+00
             enddo
             enddo
             if (ltvd) then
               if (kk .GE. 3) then
                 idrctend = 3
               else
                 idrctend = 2
               endif
               do idrct = 1, idrctend
                 idrct2 = idrct
                 do l = -1, 0
* pick up L & R : "b" and "c" for each side
                   ro1c = rob(l+1,1,idrct)
                   ro1b = rob(l+1,2,idrct)
                   ur1c = urb(l+1,1,idrct)
                   ur1b = urb(l+1,2,idrct)
                   ut1c = utb(l+1,1,idrct)
                   ut1b = utb(l+1,2,idrct)
                   up1c = upb(l+1,1,idrct)
                   up1b = upb(l+1,2,idrct)
                   mr1c = mrb(l+1,1,idrct)
                   mr1b = mrb(l+1,2,idrct)
                   mt1c = mtb(l+1,1,idrct)
                   mt1b = mtb(l+1,2,idrct)
                   mp1c = mpb(l+1,1,idrct)
                   mp1b = mpb(l+1,2,idrct)
                   br1c = brb(l+1,1,idrct)
                   br1b = brb(l+1,2,idrct)
                   bt1c = btb(l+1,1,idrct)
                   bt1b = btb(l+1,2,idrct)
                   bp1c = bpb(l+1,1,idrct)
                   bp1b = bpb(l+1,2,idrct)
                   ed1c = edb(l+1,1,idrct)
                   ed1b = edb(l+1,2,idrct)
                   rod2 = ro1c - ro1b
                   urd2 = ur1c - ur1b
                   utd2 = ut1c - ut1b
                   upd2 = up1c - up1b
                   mrd2 = mr1c - mr1b
                   mtd2 = mt1c - mt1b
                   mpd2 = mp1c - mp1b
                   brd2 = br1c - br1b
                   btd2 = bt1c - bt1b
                   bpd2 = bp1c - bp1b
                   edd2 = ed1c - ed1b
                   call drtp2num(idrct2,
     &                    mrd1,mrd2,mrd3,mtd1,mtd2,mtd3,mpd1,mpd2,mpd3,
     &                    brd1,brd2,brd3,btd1,btd2,btd3,bpd1,bpd2,bpd3,
     &                    m1d1,m1d2,m1d3,m2d1,m2d2,m2d3,m3d1,m3d2,m3d3,
     &                    b1d1,b1d2,b1d3,b2d1,b2d2,b2d3,b3d1,b3d2,b3d3)

* calculate Matrix(-rices)
                     call avehalf(gammax2, ! idrct2,
     &                        ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b,
     &                        ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c,
     &                        roa,ura,uta,upa,bra,bta,bpa,pga)
                     call vrtp2num(idrct2,
     &                        ura,uta,upa,bra,bta,bpa,
     &                        v1,v2,v3,b1,b2,b3)
                     aa = pga / roa * gammax2
                     bb =(b1**2 + b2**2 + b3**2) / roa
                     cc = b1**2 / roa
                     dd = (aa + bb)**2 - 4.0E+00 * aa * cc
                     if (dd .LE. 0.0E+00) then
                       dd = 0.0E+00
                     else
                       dd = sqrt(dd)
                     endif
                     dd = 0.5E+00 * (aa + bb + dd)
                     dd = abs(v1) + sqrt(dd) ! |V_1| + |V_f|

                     if (idrct .EQ. 1) then
                       tvdcr(1,l+1) = dd * rod2 * 0.5E+00
                       tvdcr(2,l+1) = dd * m1d2 * 0.5E+00
                       tvdcr(3,l+1) = dd * m2d2 * 0.5E+00
                       tvdcr(4,l+1) = dd * m3d2 * 0.5E+00
                       tvdcr(5,l+1) = dd * b2d2 * 0.5E+00
                       tvdcr(6,l+1) = dd * b3d2 * 0.5E+00
                       tvdcr(7,l+1) = dd * edd2 * 0.5E+00
                     else if (idrct .EQ. 2) then
                       tvdct(1,l+1) = dd * rod2 * 0.5E+00
                       tvdct(2,l+1) = dd * m1d2 * 0.5E+00
                       tvdct(3,l+1) = dd * m2d2 * 0.5E+00
                       tvdct(4,l+1) = dd * m3d2 * 0.5E+00
                       tvdct(5,l+1) = dd * b2d2 * 0.5E+00
                       tvdct(6,l+1) = dd * b3d2 * 0.5E+00
                       tvdct(7,l+1) = dd * edd2 * 0.5E+00
                     else
                       tvdcp(1,l+1) = dd * rod2 * 0.5E+00
                       tvdcp(2,l+1) = dd * m1d2 * 0.5E+00
                       tvdcp(3,l+1) = dd * m2d2 * 0.5E+00
                       tvdcp(4,l+1) = dd * m3d2 * 0.5E+00
                       tvdcp(5,l+1) = dd * b2d2 * 0.5E+00
                       tvdcp(6,l+1) = dd * b3d2 * 0.5E+00
                       tvdcp(7,l+1) = dd * edd2 * 0.5E+00
                     endif

                 enddo ! end of l = -1, 0
               enddo ! end of idrct = 1, 3
             endif ! end of (ltvd)

             if (kk .LT. 3) then ! not be defined in 2D
               tvdcp(1,0) = 0.0E+00
               tvdcp(2,0) = 0.0E+00
               tvdcp(3,0) = 0.0E+00
               tvdcp(4,0) = 0.0E+00
               tvdcp(5,0) = 0.0E+00
               tvdcp(6,0) = 0.0E+00
               tvdcp(7,0) = 0.0E+00
               tvdcp(1,1) = 0.0E+00
               tvdcp(2,1) = 0.0E+00
               tvdcp(3,1) = 0.0E+00
               tvdcp(4,1) = 0.0E+00
               tvdcp(5,1) = 0.0E+00
               tvdcp(6,1) = 0.0E+00
               tvdcp(7,1) = 0.0E+00
             endif

* auxi.
             divmag = (br1r1 * rrb1  - br1r2 * rrb2)  * ddra1
     &              + (bt1t1 * sinb1 - bt1t2 * sinb2) * ddth1
     &              + (bp1p1 - bp1p2) * ddph2

             jr1cc = (bp1t1 * sinb1 - bp1t2 * sinb2) * ddth1
     &             - (bt1p1 - bt1p2) * ddph2
             jt1cc = (br1p1 - br1p2) * ddph0
     &             - (bp1r1 - bp1r2) * ddra0
             jp1cc = (bt1r1 - bt1r2) * ddra0
     &             - (br1t1 - br1t2) * ddth0

             gravit  = 2.0E+00 * (1.0E+00/rrb1 - 1.0E+00/rrb2) * ddra0

** differencing equations, with modification terms

* RO
             aa =-urb(1,1,1) * rob(1,1,1)
             bb =-urb(1,2,1) * rob(1,2,1)
             dd = (aa + bb) * 0.5E+00 + tvdcr(1,1)
             aa =-urb(0,1,1) * rob(0,1,1)
             bb =-urb(0,2,1) * rob(0,2,1)
             ee = (aa + bb) * 0.5E+00 + tvdcr(1,0)
             xx2= (dd - ee) * ddra0

             aa =-utb(1,1,2) * rob(1,1,2)
             bb =-utb(1,2,2) * rob(1,2,2)
             dd = (aa + bb) * 0.5E+00 + tvdct(1,1)
             aa =-utb(0,1,2) * rob(0,1,2)
             bb =-utb(0,2,2) * rob(0,2,2)
             ee = (aa + bb) * 0.5E+00 + tvdct(1,0)
             xx3= (dd * sinb1 - ee * sinb2) * ddth1

             aa =-upb(1,1,3) * rob(1,1,3)
             bb =-upb(1,2,3) * rob(1,2,3)
             dd = (aa + bb) * 0.5E+00 + tvdcp(1,1)
             aa =-upb(0,1,3) * rob(0,1,3)
             bb =-upb(0,2,3) * rob(0,2,3)
             ee = (aa + bb) * 0.5E+00 + tvdcp(1,0)
             xx4= (dd - ee) * ddph2

             xxrol(idt) = xx2 + xx3 + xx4

* URTP, MRTP
!             if (ljxb) then
               aa=- urb(1,1,1) * mrb(1,1,1) - pgb(1,1,1)
               bb=- urb(1,2,1) * mrb(1,2,1) - pgb(1,2,1)
               dd = (aa + bb) * 0.5E+00 + tvdcr(2,1)
               tnsrr( 0) = dd    ! Trr at i + 1/2
               aa=- urb(0,1,1) * mrb(0,1,1) - pgb(0,1,1)
               bb=- urb(0,2,1) * mrb(0,2,1) - pgb(0,2,1)
               ee = (aa + bb) * 0.5E+00 + tvdcr(2,0)
               tnsrr(-1) = ee

               aa =-utb(1,1,2) * mrb(1,1,2)
               bb =-utb(1,2,2) * mrb(1,2,2)
               dd = (aa + bb) * 0.5E+00 + tvdct(4,1)
               tnstr( 0) = dd
               aa =-utb(0,1,2) * mrb(0,1,2)
               bb =-utb(0,2,2) * mrb(0,2,2)
               ee = (aa + bb) * 0.5E+00 + tvdct(4,0)
               tnstr(-1) = ee

               aa =-upb(1,1,3) * mrb(1,1,3)
               bb =-upb(1,2,3) * mrb(1,2,3)
               dd = (aa + bb) * 0.5E+00 + tvdcp(3,1)
               tnspr( 0) = dd
               aa =-upb(0,1,3) * mrb(0,1,3)
               bb =-upb(0,2,3) * mrb(0,2,3)
               ee = (aa + bb) * 0.5E+00 + tvdcp(3,0)
               tnspr(-1) = ee
***
               aa=- urb(1,1,1) * mtb(1,1,1)
               bb=- urb(1,2,1) * mtb(1,2,1)
               dd = (aa + bb) * 0.5E+00 + tvdcr(3,1)
               tnsrt( 0) = dd
               aa=- urb(0,1,1) * mtb(0,1,1)
               bb=- urb(0,2,1) * mtb(0,2,1)
               ee = (aa + bb) * 0.5E+00 + tvdcr(3,0)
               tnsrt(-1) = ee

               aa =- utb(1,1,2) * mtb(1,1,2) - pgb(1,1,2)
               bb =- utb(1,2,2) * mtb(1,2,2) - pgb(1,2,2)
               dd = (aa + bb) * 0.5E+00 + tvdct(2,1)
               tnstt( 0) = dd
               aa =- utb(0,1,2) * mtb(0,1,2) - pgb(0,1,2)
               bb =- utb(0,2,2) * mtb(0,2,2) - pgb(0,2,2)
               ee = (aa + bb) * 0.5E+00 + tvdct(2,0)
               tnstt(-1) = ee

               aa =-upb(1,1,3) * mtb(1,1,3)
               bb =-upb(1,2,3) * mtb(1,2,3)
               dd = (aa + bb) * 0.5E+00 + tvdcp(4,1)
               tnspt( 0) = dd
               aa =-upb(0,1,3) * mtb(0,1,3)
               bb =-upb(0,2,3) * mtb(0,2,3)
               ee = (aa + bb) * 0.5E+00 + tvdcp(4,0)
               tnspt(-1) = ee
**
               aa=- urb(1,1,1) * mpb(1,1,1)
               bb=- urb(1,2,1) * mpb(1,2,1)
               dd = (aa + bb) * 0.5E+00 + tvdcr(4,1)
               tnsrp( 0) = dd
               aa=- urb(0,1,1) * mpb(0,1,1)
               bb=- urb(0,2,1) * mpb(0,2,1)
               ee = (aa + bb) * 0.5E+00 + tvdcr(4,0)
               tnsrp(-1) = ee

               aa =-utb(1,1,2) * mpb(1,1,2)
               bb =-utb(1,2,2) * mpb(1,2,2)
               dd = (aa + bb) * 0.5E+00 + tvdct(3,1)
               tnstp( 0) = dd
               aa =-utb(0,1,2) * mpb(0,1,2)
               bb =-utb(0,2,2) * mpb(0,2,2)
               ee = (aa + bb) * 0.5E+00 + tvdct(3,0)
               tnstp(-1) = ee

               aa =- upb(1,1,3) * mpb(1,1,3) - pgb(1,1,3)
               bb =- upb(1,2,3) * mpb(1,2,3) - pgb(1,2,3)
               dd = (aa + bb) * 0.5E+00 + tvdcp(2,1)
               tnspp( 0) = dd
               aa =- upb(0,1,3) * mpb(0,1,3) - pgb(0,1,3)
               bb =- upb(0,2,3) * mpb(0,2,3) - pgb(0,2,3)
               ee = (aa + bb) * 0.5E+00 + tvdcp(2,0)
               tnspp(-1) = ee

!             else
!             endif

!             if (lvolmom) then

               if (ldiagmom) then
                 if (ljxb) then
                   aa=- urb(1,1,1) * mrb(1,1,1)
                   bb=- urb(1,2,1) * mrb(1,2,1)
                   dd = (aa + bb) * 0.5E+00 + tvdcr(2,1)
                   tnsrr( 0) = dd
                   aa=- urb(0,1,1) * mrb(0,1,1)
                   bb=- urb(0,2,1) * mrb(0,2,1)
                   ee = (aa + bb) * 0.5E+00 + tvdcr(2,0)
                   tnsrr(-1) = ee

                   aa=- utb(1,1,2) * mtb(1,1,2)
                   bb=- utb(1,2,2) * mtb(1,2,2)
                   dd = (aa + bb) * 0.5E+00 + tvdct(2,1)
                   tnstt( 0) = dd
                   aa=- utb(0,1,2) * mtb(0,1,2)
                   bb=- utb(0,2,2) * mtb(0,2,2)
                   ee = (aa + bb) * 0.5E+00 + tvdct(2,0)
                   tnstt(-1) = ee

                   aa=- upb(1,1,3) * mpb(1,1,3)
                   bb=- upb(1,2,3) * mpb(1,2,3)
                   dd = (aa + bb) * 0.5E+00 + tvdcp(2,1)
                   tnspp( 0) = dd
                   aa=- upb(0,1,3) * mpb(0,1,3)
                   bb=- upb(0,2,3) * mpb(0,2,3)
                   ee = (aa + bb) * 0.5E+00 + tvdcp(2,0)
                   tnspp(-1) = ee
                  else
                   aa=- urb(1,1,1) * mrb(1,1,1) + brb(1,1,1)**2
                   bb=- urb(1,2,1) * mrb(1,2,1) + brb(1,2,1)**2
                   dd = (aa + bb) * 0.5E+00 + tvdcr(2,1)
                   tnsrr( 0) = dd
                   aa=- urb(0,1,1) * mrb(0,1,1) + brb(0,1,1)**2
                   bb=- urb(0,2,1) * mrb(0,2,1) + brb(0,2,1)**2
                   ee = (aa + bb) * 0.5E+00 + tvdcr(2,0)
                   tnsrr(-1) = ee

                   aa=- utb(1,1,2) * mtb(1,1,2) + btb(1,1,2)**2
                   bb=- utb(1,2,2) * mtb(1,2,2) + btb(1,2,2)**2
                   dd = (aa + bb) * 0.5E+00 + tvdct(2,1)
                   tnstt( 0) = dd
                   aa=- utb(0,1,2) * mtb(0,1,2) + btb(0,1,2)**2
                   bb=- utb(0,2,2) * mtb(0,2,2) + btb(0,2,2)**2
                   ee = (aa + bb) * 0.5E+00 + tvdct(2,0)
                   tnstt(-1) = ee

                   aa=- upb(1,1,3) * mpb(1,1,3) + bpb(1,1,3)**2
                   bb=- upb(1,2,3) * mpb(1,2,3) + bpb(1,2,3)**2
                   dd = (aa + bb) * 0.5E+00 + tvdcp(2,1)
                   tnspp( 0) = dd
                   aa=- upb(0,1,3) * mpb(0,1,3) + bpb(0,1,3)**2
                   bb=- upb(0,2,3) * mpb(0,2,3) + bpb(0,2,3)**2
                   ee = (aa + bb) * 0.5E+00 + tvdcp(2,0)
                   tnspp(-1) = ee
                 endif ! end if (ljxb)
               endif

               sinph0 = sinph(k)
               cosph0 = cosph(k)
               if (kk .LT. 3) then
                 cosphb1 = cosphb(k)
                 sinphb1 = sinphb(k)
                 cosphb2 = cosphb(k-1)
                 sinphb2 = sinphb(k-1)
               else
                 if (kstepb(j) .GT. 1) then
                   kl2 = k - kstepb(j) / 2 + kk * 2
                   kl2 = mod(kl2,kk)
                   kl3 = k + kstepb(j) / 2 + kk * 2
                   kl3 = mod(kl3,kk)
                   cosphb1 = cosph(kl3)
                   sinphb1 = sinph(kl3)
                   cosphb2 = cosph(kl2)
                   sinphb2 = sinph(kl2)
                 else
                   cosphb1 = cosphb(k)
                   sinphb1 = sinphb(k)
                   cosphb2 = cosphb(k-1)
                   sinphb2 = sinphb(k-1)
                 endif
               endif

               xxmrl(idt) = (tnsrr(0) - tnsrr(-1)) * ddra0
               xxmtl(idt) = (tnsrt(0) - tnsrt(-1)) * ddra0
               xxmpl(idt) = (tnsrp(0) - tnsrp(-1)) * ddra0

               xxxx1= tnstr( 0) * sinb1 * cosph0
     &              + tnstt( 0) * cosb1 * cosph0
     &              - tnstp( 0)         * sinph0
               xxyy1= tnstr( 0) * sinb1 * sinph0
     &              + tnstt( 0) * cosb1 * sinph0
     &              + tnstp( 0)         * cosph0
               xxzz1= tnstr( 0) * cosb1
     &              - tnstt( 0) * sinb1
               xxrr1= xxxx1 * sin0 * cosph0
     &              + xxyy1 * sin0 * sinph0
     &              + xxzz1 * cos0
               xxtt1= xxxx1 * cos0 * cosph0
     &              + xxyy1 * cos0 * sinph0
     &              - xxzz1 * sin0
               xxpp1=-xxxx1        * sinph0
     &              + xxyy1        * cosph0
               xxxx2= tnstr(-1) * sinb2 * cosph0
     &              + tnstt(-1) * cosb2 * cosph0
     &              - tnstp(-1)         * sinph0
               xxyy2= tnstr(-1) * sinb2 * sinph0
     &              + tnstt(-1) * cosb2 * sinph0
     &              + tnstp(-1)         * cosph0
               xxzz2= tnstr(-1) * cosb2
     &              - tnstt(-1) * sinb2
               xxrr2= xxxx2 * sin0 * cosph0
     &              + xxyy2 * sin0 * sinph0
     &              + xxzz2 * cos0
               xxtt2= xxxx2 * cos0 * cosph0
     &              + xxyy2 * cos0 * sinph0
     &              - xxzz2 * sin0
               xxpp2=-xxxx2        * sinph0
     &              + xxyy2        * cosph0
               xxmrl(idt) =xxmrl(idt)+(xxrr1*sinb1-xxrr2*sinb2) * ddth1
               xxmtl(idt) =xxmtl(idt)+(xxtt1*sinb1-xxtt2*sinb2) * ddth1
               xxmpl(idt) =xxmpl(idt)+(xxpp1*sinb1-xxpp2*sinb2) * ddth1

               xxxx1= tnspr( 0) * sin0 * cosphb1
     &              + tnspt( 0) * cos0 * cosphb1
     &              - tnspp( 0)        * sinphb1
               xxyy1= tnspr( 0) * sin0 * sinphb1
     &              + tnspt( 0) * cos0 * sinphb1
     &              + tnspp( 0)        * cosphb1
               xxzz1= tnspr( 0) * cos0
     &              - tnspt( 0) * sin0
               xxrr1= xxxx1 * sin0 * cosph0
     &              + xxyy1 * sin0 * sinph0
     &              + xxzz1 * cos0
               xxtt1= xxxx1 * cos0 * cosph0
     &              + xxyy1 * cos0 * sinph0
     &              - xxzz1 * sin0
               xxpp1=-xxxx1        * sinph0
     &              + xxyy1        * cosph0
               xxxx2= tnspr(-1) * sin0 * cosphb2
     &              + tnspt(-1) * cos0 * cosphb2
     &              - tnspp(-1)        * sinphb2
               xxyy2= tnspr(-1) * sin0 * sinphb2
     &              + tnspt(-1) * cos0 * sinphb2
     &              + tnspp(-1)        * cosphb2
               xxzz2= tnspr(-1) * cos0
     &              - tnspt(-1) * sin0
               xxrr2= xxxx2 * sin0 * cosph0
     &              + xxyy2 * sin0 * sinph0
     &              + xxzz2 * cos0
               xxtt2= xxxx2 * cos0 * cosph0
     &              + xxyy2 * cos0 * sinph0
     &              - xxzz2 * sin0
               xxpp2=-xxxx2        * sinph0
     &              + xxyy2        * cosph0
               xxmrl(idt) = xxmrl(idt) + (xxrr1 - xxrr2) * ddph2
               xxmtl(idt) = xxmtl(idt) + (xxtt1 - xxtt2) * ddph2
               xxmpl(idt) = xxmpl(idt) + (xxpp1 - xxpp2) * ddph2

* diag. part
               if (ldiagmom) then
!                 if (ljxb) then
                   aa= - pgb(1,1,1)
                   bb= - pgb(1,2,1)
                   dd= (aa + bb) * 0.5E+00
                   aa= - pgb(0,1,1)
                   bb= - pgb(0,2,1)
                   ee= (aa + bb) * 0.5E+00
                   xx2 = (dd/rrb1**2 - ee/rrb2**2) * ddra0 * rrbc2nd
                   xxmrl(idt) = xxmrl(idt) + xx2

                   aa= - pgb(1,1,2)
                   bb= - pgb(1,2,2)
                   dd= (aa + bb) * 0.5E+00
                   aa= - pgb(0,1,2)
                   bb= - pgb(0,2,2)
                   ee= (aa + bb) * 0.5E+00
                   xx2 = (dd - ee) * ddth0
                   xxmtl(idt) = xxmtl(idt) + xx2

                   aa= - pgb(1,1,3)
                   bb= - pgb(1,2,3)
                   dd= (aa + bb) * 0.5E+00
                   aa= - pgb(0,1,3)
                   bb= - pgb(0,2,3)
                   ee= (aa + bb) * 0.5E+00
                   xx2 = (dd - ee) * ddph0
                   xxmpl(idt) = xxmpl(idt) + xx2
!                 else
!                 endif ! end if (ljxb)
               endif

!             else ! if .NOT. lvolmom Mind : check kstep() = 1
!             endif ! end of if lvolmom or not

* source term
             bb  = 2.0E+00 * omega * sin0 * up1cc ! Coliori_r
             cc  = rr(i) * omega**2 * sin0**2 !     Centri_r
             dd = gravit !                          Gravitational
             xx5 = (bb + cc + dd) * ro1cc
             xxmrl(idt) = xxmrl(idt) + xx5
             bb  = 2.0E+00 * omega * cos0 * up1cc ! Coliori_t
             cc  = rr(i) * omega**2 * sin0 * cos0 ! Centri_t
             xx5 = (bb + cc) * ro1cc
             xxmtl(idt) = xxmtl(idt) + xx5
             xx5 = -2.0E+00 * omega
     &           * (sin0 * ur1cc + cos0 * ut1cc) * ro1cc
             xxmpl(idt) = xxmpl(idt) + xx5

* add jxB
!             if (ljxb) then ! add jxb term
               xxmrl(idt) = xxmrl(idt) + jt1cc * bp1cc - jp1cc * bt1cc
               xxmtl(idt) = xxmtl(idt) + jp1cc * br1cc - jr1cc * bp1cc
               xxmpl(idt) = xxmpl(idt) + jr1cc * bt1cc - jt1cc * br1cc
!             else
* projection divB sweeper
!               xxmrl(idt) = xxmrl(idt) - divmag * br1cc
!               xxmtl(idt) = xxmtl(idt) - divmag * bt1cc
!               xxmpl(idt) = xxmpl(idt) - divmag * bp1cc
!             endif

* calculating EMF
             aa =-utb(1,1,2) * brb(1,1,2) + btb(1,1,2) * urb(1,1,2)
             bb =-utb(1,2,2) * brb(1,2,2) + btb(1,2,2) * urb(1,2,2)
             dd = (aa + bb) * 0.5E+00 + tvdct(6,1)
             tnstr( 0) = dd
             aa =-utb(0,1,2) * brb(0,1,2) + btb(0,1,2) * urb(0,1,2)
             bb =-utb(0,2,2) * brb(0,2,2) + btb(0,2,2) * urb(0,2,2)
             ee = (aa + bb) * 0.5E+00 + tvdct(6,0)
             tnstr(-1) = ee

             aa =-upb(1,1,3) * brb(1,1,3) + bpb(1,1,3) * urb(1,1,3)
             bb =-upb(1,2,3) * brb(1,2,3) + bpb(1,2,3) * urb(1,2,3)
             dd = (aa + bb) * 0.5E+00 + tvdcp(5,1)
             tnspr( 0) = dd
             aa =-upb(0,1,3) * brb(0,1,3) + bpb(0,1,3) * urb(0,1,3)
             bb =-upb(0,2,3) * brb(0,2,3) + bpb(0,2,3) * urb(0,2,3)
             ee = (aa + bb) * 0.5E+00 + tvdcp(5,0)
             tnspr(-1) = ee
***
             aa =-urb(1,1,1) * btb(1,1,1) + brb(1,1,1) * utb(1,1,1)
             bb =-urb(1,2,1) * btb(1,2,1) + brb(1,2,1) * utb(1,2,1)
             dd = (aa + bb) * 0.5E+00 + tvdcr(5,1)
             tnsrt( 0) = dd
             aa =-urb(0,1,1) * btb(0,1,1) + brb(0,1,1) * utb(0,1,1)
             bb =-urb(0,2,1) * btb(0,2,1) + brb(0,2,1) * utb(0,2,1)
             ee = (aa + bb) * 0.5E+00 + tvdcr(5,0)
             tnsrt(-1) = ee

             aa =-upb(1,1,3) * btb(1,1,3) + bpb(1,1,3) * utb(1,1,3)
             bb =-upb(1,2,3) * btb(1,2,3) + bpb(1,2,3) * utb(1,2,3)
             dd = (aa + bb) * 0.5E+00 + tvdcp(6,1)
             tnspt( 0) = dd
             aa =-upb(0,1,3) * btb(0,1,3) + bpb(0,1,3) * utb(0,1,3)
             bb =-upb(0,2,3) * btb(0,2,3) + bpb(0,2,3) * utb(0,2,3)
             ee = (aa + bb) * 0.5E+00 + tvdcp(6,0)
             tnspt(-1) = ee
**
             aa =-urb(1,1,1) * bpb(1,1,1) + brb(1,1,1) * upb(1,1,1)
             bb =-urb(1,2,1) * bpb(1,2,1) + brb(1,2,1) * upb(1,2,1)
             dd = (aa + bb) * 0.5E+00 + tvdcr(6,1)
             tnsrp( 0) = dd
             aa =-urb(0,1,1) * bpb(0,1,1) + brb(0,1,1) * upb(0,1,1)
             bb =-urb(0,2,1) * bpb(0,2,1) + brb(0,2,1) * upb(0,2,1)
             ee = (aa + bb) * 0.5E+00 + tvdcr(6,0)
             tnsrp(-1) = ee

             aa =-utb(1,1,2) * bpb(1,1,2) + btb(1,1,2) * upb(1,1,2)
             bb =-utb(1,2,2) * bpb(1,2,2) + btb(1,2,2) * upb(1,2,2)
             dd = (aa + bb) * 0.5E+00 + tvdct(5,1)
             tnstp( 0) = dd
             aa =-utb(0,1,2) * bpb(0,1,2) + btb(0,1,2) * upb(0,1,2)
             bb =-utb(0,2,2) * bpb(0,2,2) + btb(0,2,2) * upb(0,2,2)
             ee = (aa + bb) * 0.5E+00 + tvdct(5,0)
             tnstp(-1) = ee

** EMF at i = 0
!             if (lemfzero) then
!             endif

* BR, BT, BP
!             if (lvolmag) then

               tnsrr( 0) = 0.0E+00 ! None !
               tnsrr(-1) = 0.0E+00
               tnstt( 0) = 0.0E+00
               tnstt(-1) = 0.0E+00
               tnspp( 0) = 0.0E+00
               tnspp(-1) = 0.0E+00

               xxbrl(idt) = (tnsrr(0)*rrb1 - tnsrr(-1)*rrb2) * ddra1 ! <=== added
               xxbtl(idt) = (tnsrt(0)*rrb1 - tnsrt(-1)*rrb2) * ddra1
               xxbpl(idt) = (tnsrp(0)*rrb1 - tnsrp(-1)*rrb2) * ddra1

               sinph0 = sinph(k)
               cosph0 = cosph(k)
               if (kk .LT. 3) then
                 cosphb1 = cosphb(k)
                 sinphb1 = sinphb(k)
                 cosphb2 = cosphb(k-1)
                 sinphb2 = sinphb(k-1)
               else
                 if (kstepb(j) .GT. 1) then
                   kl2 = k - kstepb(j) / 2 + kk * 2
                   kl2 = mod(kl2,kk)
                   kl3 = k + kstepb(j) / 2 + kk * 2
                   kl3 = mod(kl3,kk)
                   cosphb1 = cosph(kl3)
                   sinphb1 = sinph(kl3)
                   cosphb2 = cosph(kl2)
                   sinphb2 = sinph(kl2)
                 else
                   cosphb1 = cosphb(k)
                   sinphb1 = sinphb(k)
                   cosphb2 = cosphb(k-1)
                   sinphb2 = sinphb(k-1)
                 endif
               endif

               xxxx1= tnstr( 0) * sinb1 * cosph0
     &              + tnstt( 0) * cosb1 * cosph0
     &              - tnstp( 0)         * sinph0
               xxyy1= tnstr( 0) * sinb1 * sinph0
     &              + tnstt( 0) * cosb1 * sinph0
     &              + tnstp( 0)         * cosph0
               xxzz1= tnstr( 0) * cosb1
     &              - tnstt( 0) * sinb1
               xxrr1= xxxx1 * sin0 * cosph0
     &              + xxyy1 * sin0 * sinph0
     &              + xxzz1 * cos0
               xxtt1= xxxx1 * cos0 * cosph0
     &              + xxyy1 * cos0 * sinph0
     &              - xxzz1 * sin0
               xxpp1=-xxxx1        * sinph0
     &              + xxyy1        * cosph0
               xxxx2= tnstr(-1) * sinb2 * cosph0
     &              + tnstt(-1) * cosb2 * cosph0
     &              - tnstp(-1)         * sinph0
               xxyy2= tnstr(-1) * sinb2 * sinph0
     &              + tnstt(-1) * cosb2 * sinph0
     &              + tnstp(-1)         * cosph0
               xxzz2= tnstr(-1) * cosb2
     &              - tnstt(-1) * sinb2
               xxrr2= xxxx2 * sin0 * cosph0
     &              + xxyy2 * sin0 * sinph0
     &              + xxzz2 * cos0
               xxtt2= xxxx2 * cos0 * cosph0
     &              + xxyy2 * cos0 * sinph0
     &              - xxzz2 * sin0
               xxpp2=-xxxx2        * sinph0
     &              + xxyy2        * cosph0
               xxbrl(idt) = xxbrl(idt)+(xxrr1*sinb1-xxrr2*sinb2)*ddth1
               xxbtl(idt) = xxbtl(idt)+(xxtt1*sinb1-xxtt2*sinb2)*ddth1
               xxbpl(idt) = xxbpl(idt)+(xxpp1*sinb1-xxpp2*sinb2)*ddth1

               xxxx1= tnspr( 0) * sin0 * cosphb1
     &              + tnspt( 0) * cos0 * cosphb1
     &              - tnspp( 0)        * sinphb1
               xxyy1= tnspr( 0) * sin0 * sinphb1
     &              + tnspt( 0) * cos0 * sinphb1
     &              + tnspp( 0)        * cosphb1
               xxzz1= tnspr( 0) * cos0
     &              - tnspt( 0) * sin0
               xxrr1= xxxx1 * sin0 * cosph0
     &              + xxyy1 * sin0 * sinph0
     &              + xxzz1 * cos0
               xxtt1= xxxx1 * cos0 * cosph0
     &              + xxyy1 * cos0 * sinph0
     &              - xxzz1 * sin0
               xxpp1=-xxxx1       * sinph0
     &              + xxyy1       * cosph0
               xxxx2= tnspr(-1) * sin0 * cosphb2
     &              + tnspt(-1) * cos0 * cosphb2
     &              - tnspp(-1)        * sinphb2
               xxyy2= tnspr(-1) * sin0 * sinphb2
     &              + tnspt(-1) * cos0 * sinphb2
     &              + tnspp(-1)        * cosphb2
               xxzz2= tnspr(-1) * cos0
     &              - tnspt(-1) * sin0
               xxrr2= xxxx2 * sin0 * cosph0
     &              + xxyy2 * sin0 * sinph0
     &              + xxzz2 * cos0
               xxtt2= xxxx2 * cos0 * cosph0
     &              + xxyy2 * cos0 * sinph0
     &              - xxzz2 * sin0
               xxpp2=-xxxx2       * sinph0
     &              + xxyy2       * cosph0
               xxbrl(idt) = xxbrl(idt) + (xxrr1 - xxrr2) * ddph2
               xxbtl(idt) = xxbtl(idt) + (xxtt1 - xxtt2) * ddph2
               xxbpl(idt) = xxbpl(idt) + (xxpp1 - xxpp2) * ddph2

* projection. divB sweeper,
               xxbrl(idt) = xxbrl(idt) - divmag * ur1cc
               xxbtl(idt) = xxbtl(idt) - divmag * ut1cc
               xxbpl(idt) = xxbpl(idt) - divmag * up1cc

!             else ! if not lvolmag
!             endif

** ED, PG
             if (lconene) then
               aa =- urb(1,1,1) * (edb(1,1,1) + pgb(1,1,1)
     &                           +(brb(1,1,1)**2
     &                            +btb(1,1,1)**2
     &                            +bpb(1,1,1)**2)*0.5E+00)
     &             + brb(1,1,1) * (brb(1,1,1) * urb(1,1,1)
     &                            +btb(1,1,1) * utb(1,1,1)
     &                            +bpb(1,1,1) * upb(1,1,1))
               bb =- urb(1,2,1) * (edb(1,2,1) + pgb(1,2,1)
     &                           +(brb(1,2,1)**2
     &                            +btb(1,2,1)**2
     &                            +bpb(1,2,1)**2)*0.5E+00)
     &             + brb(1,2,1) * (brb(1,2,1) * urb(1,2,1)
     &                            +btb(1,2,1) * utb(1,2,1)
     &                            +bpb(1,2,1) * upb(1,2,1))
               dd = (aa + bb) * 0.5E+00 + tvdcr(7,1)
               aa =- urb(0,1,1) * (edb(0,1,1) + pgb(0,1,1)
     &                           +(brb(0,1,1)**2
     &                            +btb(0,1,1)**2
     &                            +bpb(0,1,1)**2)*0.5E+00)
     &             + brb(0,1,1) * (brb(0,1,1) * urb(0,1,1)
     &                            +btb(0,1,1) * utb(0,1,1)
     &                            +bpb(0,1,1) * upb(0,1,1))
               bb =- urb(0,2,1) * (edb(0,2,1) + pgb(0,2,1)
     &                           +(brb(0,2,1)**2
     &                            +btb(0,2,1)**2
     &                            +bpb(0,2,1)**2)*0.5E+00)
     &             + brb(0,2,1) * (brb(0,2,1) * urb(0,2,1)
     &                            +btb(0,2,1) * utb(0,2,1)
     &                            +bpb(0,2,1) * upb(0,2,1))
               ee = (aa + bb) * 0.5E+00 + tvdcr(7,0)
               xx2= (dd - ee) * ddra0 ! dra2*rr(i)**2

               aa =- utb(1,1,2) * (edb(1,1,2) + pgb(1,1,2)
     &                           +(brb(1,1,2)**2
     &                            +btb(1,1,2)**2
     &                            +bpb(1,1,2)**2)*0.5E+00)
     &             + btb(1,1,2) * (brb(1,1,2) * urb(1,1,2)
     &                            +btb(1,1,2) * utb(1,1,2)
     &                            +bpb(1,1,2) * upb(1,1,2))
               bb =- utb(1,2,2) * (edb(1,2,2) + pgb(1,2,2)
     &                           +(brb(1,2,2)**2
     &                            +btb(1,2,2)**2
     &                            +bpb(1,2,2)**2)*0.5E+00)
     &             + btb(1,2,2) * (brb(1,2,2) * urb(1,2,2)
     &                            +btb(1,2,2) * utb(1,2,2)
     &                            +bpb(1,2,2) * upb(1,2,2))
               dd = (aa + bb) * 0.5E+00 + tvdct(7,1)
               aa =- utb(0,1,2) * (edb(0,1,2) + pgb(0,1,2)
     &                           +(brb(0,1,2)**2
     &                            +btb(0,1,2)**2
     &                            +bpb(0,1,2)**2)*0.5E+00)
     &             + btb(0,1,2) * (brb(0,1,2) * urb(0,1,2)
     &                            +btb(0,1,2) * utb(0,1,2)
     &                            +bpb(0,1,2) * upb(0,1,2))
               bb =- utb(0,2,2) * (edb(0,2,2) + pgb(0,2,2)
     &                           +(brb(0,2,2)**2
     &                            +btb(0,2,2)**2
     &                            +bpb(0,2,2)**2)*0.5E+00)
     &             + btb(0,2,2) * (brb(0,2,2) * urb(0,2,2)
     &                            +btb(0,2,2) * utb(0,2,2)
     &                            +bpb(0,2,2) * upb(0,2,2))
               ee = (aa + bb) * 0.5E+00 + tvdct(7,0)
               xx3= (dd * sinb1 - ee * sinb2) * ddth1

               aa =- upb(1,1,3) * (edb(1,1,3) + pgb(1,1,3)
     &                           +(brb(1,1,3)**2
     &                            +btb(1,1,3)**2
     &                            +bpb(1,1,3)**2)*0.5E+00)
     &             + bpb(1,1,3) * (brb(1,1,3) * urb(1,1,3)
     &                            +btb(1,1,3) * utb(1,1,3)
     &                            +bpb(1,1,3) * upb(1,1,3))
               bb =- upb(1,2,3) * (edb(1,2,3) + pgb(1,2,3)
     &                           +(brb(1,2,3)**2
     &                            +btb(1,2,3)**2
     &                            +bpb(1,2,3)**2)*0.5E+00)
     &             + bpb(1,2,3) * (brb(1,2,3) * urb(1,2,3)
     &                            +btb(1,2,3) * utb(1,2,3)
     &                            +bpb(1,2,3) * upb(1,2,3))
               dd = (aa + bb) * 0.5E+00 + tvdcp(7,1)
               aa =- upb(0,1,3) * (edb(0,1,3) + pgb(0,1,3)
     &                           +(brb(0,1,3)**2
     &                            +btb(0,1,3)**2
     &                            +bpb(0,1,3)**2)*0.5E+00)
     &             + bpb(0,1,3) * (brb(0,1,3) * urb(0,1,3)
     &                            +btb(0,1,3) * utb(0,1,3)
     &                            +bpb(0,1,3) * upb(0,1,3))
               bb =- upb(0,2,3) * (edb(0,2,3) + pgb(0,2,3)
     &                           +(brb(0,2,3)**2
     &                            +btb(0,2,3)**2
     &                            +bpb(0,2,3)**2)*0.5E+00)
     &             + bpb(0,2,3) * (brb(0,2,3) * urb(0,2,3)
     &                            +btb(0,2,3) * utb(0,2,3)
     &                            +bpb(0,2,3) * upb(0,2,3))
               ee = (aa + bb) * 0.5E+00 + tvdcp(7,0)
               xx4= (dd - ee) * ddph2

               cc =  rr(i) * omega**2 * sin0**2 !     Centri_r
               dd = gravit !                          Gravitational
               ee =  rr(i) * omega**2 * sin0 * cos0 ! Centri_t
               xx5= ((cc + dd) * ur1cc + ee * ut1cc) * ro1cc

!               if (lmodcone) then
!               else ! usual powell-correction

                 modmagen = - divmag *(br1cc*ur1cc
     &                               + bt1cc*ut1cc
     &                               + bp1cc*up1cc)
!               endif
             endif

             xxedl(idt) = xx2 + xx3 + xx4 + xx5 + modmagen

           enddo ! end of idt - loop ! -----------------------------------------

             xxro = xxrol(2)
             xxmr = xxmrl(2)
             xxmt = xxmtl(2)
             xxmp = xxmpl(2)
             xxbr = xxbrl(2)
             xxbt = xxbtl(2)
             xxbp = xxbpl(2)
             xxed = xxedl(2)

             xxpghd=-1.0E+10 ! should not be used/referred

* update values at grids, j,k,i
           if (i .EQ. 0) then

             rhsall(1,j,k) = xxro * rr2ndmod
             rhsall(2,j,k) = xxmr * rr2ndmod
             rhsall(3,j,k) = xxmt * rr2ndmod
             rhsall(4,j,k) = xxmp * rr2ndmod
             rhsall(5,j,k) = xxbt * rr1stmod
             rhsall(6,j,k) = xxbp * rr1stmod
             rhsall(7,j,k) = xxed * rr2ndmod
             dbrdt(j,k)    = xxbr * rr1stmod

           else

             mrltmp = mr1(j,k,i) + dt * xxmr * drr2 * rr2ndmod
             mtltmp = mt1(j,k,i) + dt * xxmt * drr2 * rr2ndmod
             mpltmp = mp1(j,k,i) + dt * xxmp * drr2 * rr2ndmod
             brltmp = br1(j,k,i) + dt * xxbr * drr1 * rr1stmod
             btltmp = bt1(j,k,i) + dt * xxbt * drr1 * rr1stmod
             bpltmp = bp1(j,k,i) + dt * xxbp * drr1 * rr1stmod
             roltmp = ro1(j,k,i) + dt * xxro * drr2 * rr2ndmod
             edltmp = ed1(j,k,i) + dt * xxed * drr2 * rr2ndmod
               pgltmp =(gammax2 - 1.0E+00)
     &                *(edltmp
     &                -(mrltmp**2 + mtltmp**2 + mpltmp**2)
     &                *0.5E+00/roltmp
     &                -(brltmp**2 + btltmp**2 + bpltmp**2)
     &                *0.5E+00)
             urltmp = mrltmp / roltmp
             utltmp = mtltmp / roltmp
             upltmp = mpltmp / roltmp
* restore to dim.array
             ro9(j,k,i) = roltmp
             pg9(j,k,i) = pgltmp
             ur9(j,k,i) = urltmp
             ut9(j,k,i) = utltmp
             up9(j,k,i) = upltmp
             br9(j,k,i) = brltmp
             bt9(j,k,i) = btltmp
             bp9(j,k,i) = bpltmp

             rhsall(1,j,k) = xxro ! must not be used.
             rhsall(2,j,k) = xxmr
             rhsall(3,j,k) = xxmt
             rhsall(4,j,k) = xxmp
             rhsall(5,j,k) = xxbt
             rhsall(6,j,k) = xxbp
             rhsall(7,j,k) = xxed
           endif

         enddo
         enddo ! end of j,k-loop

* compatibility relation, or, characteristic method at the inner subsonic boundary
         if (i .EQ. 0) then

           idrct2 = 1
           do k = 1, kk
           do j = 1, jj - 1 ! MIND that prog. below is not available just at axis
*
             lplasmod = .false.

             if (lchrsurf) then

* if mass flux is postion-variable
             if (lmflxvar) then
               aa = br2(j,k,i)
               bb = b0
               cc = n0
               dd = v0
               call mflxvsmg(mcrtrn2,vcrtrn2,aa,bb,cc,dd)
             else
               vcrtrn2 = vcrtrn
               mcrtrn2 = mcrtrn
             endif

* some magnetic field map info.
               mgstr = br0(j,k,i)**2 + bt0(j,k,i)**2 + bp0(j,k,i)**2
               mgstr = sqrt(mgstr)
               aa = mgstr**2 / pg0(j,k,i)
               lweakmag = (aa .LT. 1.0E-04)

               brbb  = abs(br2(j,k,i)) / (mgstr + 1.0E-05)
               aa = br2(j,k,i)**2 + bt2(j,k,i)**2 + bp2(j,k,i)**2
               aa = sqrt(aa)
               aa = abs(br2(j,k,i)) / (aa + 1.0E-05)
               cc = 1.0E+00
*               if (br2(j,k,i)*br2(j+2,k  ,i) .LE. 0.0E+00) cc=-1.0E+00 ! broader
               if (br2(j,k,i)*br2(j+1,k+1,i) .LE. 0.0E+00) cc=-1.0E+00
               if (br2(j,k,i)*br2(j+1,k  ,i) .LE. 0.0E+00) cc=-1.0E+00
               if (br2(j,k,i)*br2(j+1,k-1,i) .LE. 0.0E+00) cc=-1.0E+00
               if (br2(j,k,i)*br2(j  ,k-1,i) .LE. 0.0E+00) cc=-1.0E+00
*               if (br2(j,k,i)*br2(j  ,k-2,i) .LE. 0.0E+00) cc=-1.0E+00 ! broader
*               if (br2(j,k,i)*br2(j  ,k+2,i) .LE. 0.0E+00) cc=-1.0E+00 ! broader
               if (br2(j,k,i)*br2(j  ,k+1,i) .LE. 0.0E+00) cc=-1.0E+00
               if (br2(j,k,i)*br2(j-1,k-1,i) .LE. 0.0E+00) cc=-1.0E+00
               if (br2(j,k,i)*br2(j-1,k  ,i) .LE. 0.0E+00) cc=-1.0E+00
               if (br2(j,k,i)*br2(j-1,k+1,i) .LE. 0.0E+00) cc=-1.0E+00
*               if (br2(j,k,i)*br2(j-2,k  ,i) .LE. 0.0E+00) cc=-1.0E+00 ! broader
               lfarmnl = ((aa .GT. brabsb) .AND. (cc .GT. 0.0E+00))
*               lfarmnl = (aa .GT. brabsb)

* try to recover the standard reference state
               if (mr2(j,k,i) .LE. mrzero) then
                 lplasmod = .false.
                 dpgdtmod = 0.0E+00
                 drodtmod = 0.0E+00
               else if (mr2(j,k,i) .LT. mcrtrn2 * 0.999E+00) then
                 lplasmod = .true.
                 dpgdtmod = pg0(j,k,i) - pg2(j,k,i)
                 drodtmod = ro0(j,k,i) - ro2(j,k,i)
               else
!                 if (ichrsurf .EQ. 3) then !  N be fixed at N_0
!                 else if (ichrsurf .EQ. 1) then ! polytrope with gammax2 be fixed
                 if (ichrsurf .EQ. 1) then ! polytrope with gammax2 be fixed
                   lplasmod = .true.
                   drodtmod = 0.0E+00
                   dpgdtmod =  ro2(j,k,i)**gammax2
     &                      *(pg0(j,k,i)/ro0(j,k,i)**gammax2
     &                       -pg2(j,k,i)/ro2(j,k,i)**gammax2)
                 else if (ichrsurf .EQ. 6) then ! polytrope with gamma6 be fixed
                   lplasmod = .true.
                   drodtmod = 0.0E+00
                   dpgdtmod = ro2(j,k,i)**gamma6
     &                      *(pg0(j,k,i)/ro0(j,k,i)**gamma6
     &                       -pg2(j,k,i)/ro2(j,k,i)**gamma6)
!                 else if (ichrsurf .EQ. 4) then
!                 else if (ichrsurf .EQ. 7) then
                 else !                  for ichrsurf is 2, 5
                   lplasmod = .false.
                   dpgdtmod = 0.0E+00
                   drodtmod = 0.0E+00
                 endif
               endif
* amplifying
               drodtmod = drodtmod * ur2(j,k,i) * rr(i)**2 * 1.0E+02
               dpgdtmod = dpgdtmod * ur2(j,k,i) * rr(i)**2 * 1.0E+02
*               drodtmod = drodtmod * rr(i)**2 * 5.0E+02
*               dpgdtmod = dpgdtmod * rr(i)**2 * 5.0E+02

               if ((.NOT. lfarmnl) .OR. (lnorefbot)) then
                 lplasmod = .false.   ! make sure everything freeze
                 dpgdtmod = 0.0E+00
                 drodtmod = 0.0E+00
               endif

* override temperature/pressure  modification if needed.
               if (lmodmap) then
                 lplasmod = .true.
** when T-map is used
                 aa = 1.0E+00 / (timemod * 3600.0E+00 / t0)
                 bb = pg2(j,k,i)
                 cc = ro2(j,k,i) * tmpf0(j,k) / gamma0 ! assume Tc = 1 MK
                 cc =(cc - bb) * aa * rr(i)**2
                 dpgdtmod = cc
                 drodtmod = 0.0E+00 ! assuming ichrsurf is 0 or 6.
** when N-map is used
*                 aa = 1.0E+00 / (timemod * 3600.0E+00 / t0)
*                 bb = ro2(j,k,i)
*                 cc = denf0(j,k) / rhoc
*                 cc =(cc - bb) * aa * rr(i)**2
*                 drodtmod = cc
*                 dpgdtmod = cc * pg2(j,k,i) / ro2(j,k,i) * gammax2
**                 dpgdtmod = cc * pg2(j,k,i) / ro2(j,k,i) * gamma0 ! test
               endif

* make RHS consistent with modified dpg & drho
               if (lplasmod) then
                 dd = drodtmod
                 rhsall(1,j,k) = rhsall(1,j,k) - dd
                 dd = dpgdtmod * dgammax2
     &              + drodtmod
     &              * 0.5E+00*(ur2(j,k,i)**2
     &                        +ut2(j,k,i)**2 + up2(j,k,i)**2)
                 rhsall(7,j,k) = rhsall(7,j,k) - dd
               endif

* make RHS consistent with dmag
!               if (ldmagon) then
!               endif

* first, solving characteristics assuming 0 < Vr < Vc
               roa = ro2(j,k,i) * rr(i)**2
               mra = mr2(j,k,i) * rr(i)**2
               mta = mt2(j,k,i) * rr(i)**2
               mpa = mp2(j,k,i) * rr(i)**2
               bra = br2(j,k,i) * rr(i)
               bta = bt2(j,k,i) * rr(i)
               bpa = bp2(j,k,i) * rr(i)
               eda = ed2(j,k,i) * rr(i)**2
               pga = pg2(j,k,i) * rr(i)**2
               v1 = mra / roa
               v2 = mta / roa
               v3 = mpa / roa
               b1 = bra
               b2 = bta
               b3 = bpa

                 call gtmatl(leigemat,
     &                       roa,pga,v1,v2,v3,b1,b2,b3,gammax2)
                 do m = 1, 7 !  pick up raw of matrix with negative eigen value
                   matl74(m,1) = leigemat(m,3) ! V - V_A
                   matl74(m,2) = leigemat(m,5) ! V - V_F
                   matl74(m,3) = leigemat(m,7) ! V - V_S
                   matl74(m,4) = leigemat(m,1) ! V... ususally >= 0, therefore, not used
                   rhsall2(m)  = rhsall(m,j,k) ! mind this comes from using ww2(j,k,i),,,
                 enddo

                 if (lweakmag) then ! if B << 1

                   n = 2 ! u - V_f = u - a
                   aa = 0.0E+00
                   do m = 1, 7
                     aa = aa + matl74(m,n) * rhsall2(m)
                   enddo
                   ansvec1(1) = aa / matl74(2,n)
                   ansvec1(2) = 0.0E+00
                   ansvec1(3) = 0.0E+00

                   xxro = 0.0E+00
                   xxpg = 0.0E+00
                   xxmr = ansvec1(1) * drr2
                   xxmt = 0.0E+00
                   xxmp = 0.0E+00
                   xxbr = 0.0E+00
                   xxbt = 0.0E+00
                   xxbp = 0.0E+00
                   xxed = 0.0E+00

                 else if (lfarmnl) then

                   do n = 1, 3 !   condition v // B., for v_r, B_theta, B_phi
                     matl33(1,n)= matl74(2,n)
     &                          + matl74(3,n) * b2 / b1
     &                          + matl74(4,n) * b3 / b1
     &                          + matl74(7,n) *(v1
     &                                         +v2 * b2 / b1
     &                                         +v3 * b3 / b1)
                     matl33(2,n)= matl74(3,n) * mra / b1
     &                          + matl74(5,n)
     &                          + matl74(7,n) *(mra / b1 * v2 + b2)
                     matl33(3,n)= matl74(4,n) * mra / b1
     &                          + matl74(6,n)
     &                          + matl74(7,n) *(mra / b1 * v3 + b3)
                   enddo

                   j2 = j
                   k2 = k
                   iloc = 1
                   call gtans33(matl33,matl74,rhsall2,ansvec1,
     &                          j2,k2,iprc2,ncal2,iloc,
     &                          roa,pga,v1,v2,v3,b1,b2,b3,
     &                          ldummy)
                   if (.NOT. ldummy) nabnorm = nabnorm + 1

                   xxro = 0.0E+00
                   xxpg = 0.0E+00
                   xxmr = ansvec1(1) * drr2
                   xxbr = 0.0E+00
                   xxbt = ansvec1(2) * drr1
                   xxbp = ansvec1(3) * drr1
                   xxmt =(mr2(j,k,i)*xxbt + xxmr*bt2(j,k,i))/br2(j,k,i)
                   xxmp =(mr2(j,k,i)*xxbp + xxmr*bp2(j,k,i))/br2(j,k,i)
                   xxed = xxmr * ur2(j,k,i)
     &                  + xxmt * ut2(j,k,i)
     &                  + xxmp * up2(j,k,i)
     &                  + xxbt * bt2(j,k,i)
     &                  + xxbp * bp2(j,k,i)

                 else ! near the neutral line

                   xxro = 0.0E+00
                   xxpg = 0.0E+00
                   xxmr =-1.0E+20 ! enforce Mr < 0
                   xxmt = 0.0E+00
                   xxmp = 0.0E+00
                   xxbr = 0.0E+00
                   xxbt = 0.0E+00
                   xxbp = 0.0E+00
                   xxed = 0.0E+00

                 endif

!               if ((ichrsurf .EQ. 4) .OR. (ichrsurf .EQ. 7)) then
!               else if (ichrsurf .EQ. 2) then
               if (ichrsurf .EQ. 0) then
                 mrdummy = mr1(j,k,i) + xxmr * dt ! calculate provisional Mr
                 if (mrdummy .LT. mrzero) then
                   iclass = 0
                 else
                   iclass = 1
                 endif
                 if (lnorefbot) iclass = 1 ! no Mr-limit
               else
                 mrdummy = mr1(j,k,i) + xxmr * dt ! calculate provisional Mr
                 if (mrdummy .LT. mrzero) then
                   iclass = 0
                 else if (mrdummy .LT. mcrtrn2) then
                   iclass = 1
                 else
                   iclass = 2
                 endif
                 if (lnorefbot) iclass = 1 ! no Mr-limit
               endif

* calculate dW/dt for each cases
               if (lweakmag) then ! if B << 1, ususally only in test mode for Parker.

                   roltmp = ro1(j,k,i)
                   mrltmp = mr1(j,k,i) + dt * xxmr
                   mtltmp = 0.0E+00
                   mpltmp = 0.0E+00
                   brltmp = 0.0E+00
                   btltmp = 0.0E+00
                   bpltmp = 0.0E+00
                   pgltmp = pg1(j,k,i)

               else if (iclass .EQ. 1) then ! |B| >> 0 and ...
* case (1), 0 < v < vc
                 roltmp = ro1(j,k,i) + dt * xxro
                 mrltmp = mr1(j,k,i) + dt * xxmr
                 mtltmp = mt1(j,k,i) + dt * xxmt
                 mpltmp = mp1(j,k,i) + dt * xxmp
                 brltmp = br0(j,k,i)
                 btltmp = bt1(j,k,i) + dt * xxbt
                 bpltmp = bp1(j,k,i) + dt * xxbp

                   pgltmp = pg1(j,k,i) + dt * xxpg

               else if ((iclass .EQ. 0) .OR. (.NOT. lfarmnl)) then ! |B| >> 0 and ...
* case (2), v < 0, stagnant
                 if (lfarmnl) then
                   if (mr1(j,k,i) .LE. mrzero) then
                     aa = 0.0E+00
                   else
                     aa = (mr1(j,k,i) - mrzero)
     &                  / (mr1(j,k,i) - mrdummy)
                   endif
                   if (aa .LT. 0.0E+00) aa = 0.0E+00
                   if (aa .GT. 1.0E+00) aa = 1.0E+00
                   dtlocal = dt * aa
                   roltmp = ro1(j,k,i) + xxro * dtlocal
                   mrltmp = 0.0E+00
                   mtltmp = 0.0E+00
                   mpltmp = 0.0E+00
                   brltmp = br0(j,k,i)
                   btltmp = bt1(j,k,i) + xxbt * dtlocal
                   bpltmp = bp1(j,k,i) + xxbp * dtlocal
                   edltmp = ed1(j,k,i) + xxed * dtlocal
                   pgltmp =(gammax2 - 1.0E+00)
     &                    *(edltmp
     &                    -(brltmp**2 + btltmp**2 + bpltmp**2)
     &                    * 0.5E+00)

                 else ! near the boundary, do nothing except enforcing v=0

                     roltmp = ro1(j,k,i)
                     mrltmp = 0.0E+00
                     mtltmp = 0.0E+00
                     mpltmp = 0.0E+00
                     brltmp = br0(j,k,i)
                     btltmp = bt1(j,k,i)
                     bpltmp = bp1(j,k,i)
                     pgltmp = pg1(j,k,i)
                     edltmp = pgltmp * dgammax2
     &                      +(brltmp**2 + btltmp**2 + bpltmp**2)
     &                      * 0.5E+00

                   dtlocal = 0.0E+00 !  near MNL, do nothing.

                 endif
                 urltmp = 0.0E+00
                 utltmp = 0.0E+00
                 upltmp = 0.0E+00

                 dtlocal = dt - dtlocal

                 roa = roltmp * rr(i)**2
                 mra = 0.0E+00
                 mta = 0.0E+00
                 mpa = 0.0E+00
                 bra = brltmp * rr(i)
                 bta = btltmp * rr(i)
                 bpa = bpltmp * rr(i)
                 eda = edltmp * rr(i)**2
                 pga = pgltmp * rr(i)**2
                 v1 = 0.0E+00
                 v2 = 0.0E+00
                 v3 = 0.0E+00
                 b1 = bra
                 b2 = bta
                 b3 = bpa
                 call gtmatl(leigemat,
     &                       roa,pga,v1,v2,v3,b1,b2,b3,gammax2)
                 do m = 1, 7 !  pick up raw of matrix with negative eigen value
                   matl74(m,1) = leigemat(m,3) ! V - V_A
                   matl74(m,2) = leigemat(m,5) ! V - V_F
                   matl74(m,3) = leigemat(m,7) ! V - V_S
                   matl74(m,4) = leigemat(m,1) ! V... ususally not used
                   rhsall2(m)  = rhsall(m,j,k) ! mind this comes from using ww2(j,k,i),,,
                 enddo
                 aa = (gammax2 - 1.0E+00) / gammab * roa / pga
                 do n = 1, 3  ! B_theta, B_phi, Ene.
                   matl33(1,n)= matl74(5,n)
     &                        - matl74(1,n) * aa * b2
                   matl33(2,n)= matl74(6,n)
     &                        - matl74(1,n) * aa * b3
                   matl33(3,n)= matl74(7,n)
     &                        + matl74(1,n) * aa
                 enddo
*
                 j2 = j
                 k2 = k
                 iloc = 2
                 call gtans33(matl33,matl74,rhsall2,ansvec2,
     &                        j2,k2,iprc2,ncal2,iloc,
     &                        roa,pga,v1,v2,v3,b1,b2,b3,
     &                        ldummy)
                 if (.NOT. ldummy) nabnorm = nabnorm + 1

                 xxbt = ansvec2(1) * drr1
                 xxbp = ansvec2(2) * drr1
                 xxed = ansvec2(3) * drr2
                 xxpg = (gammax2 - 1.0E+00)
     &                * (xxed - btltmp * xxbt - bpltmp * xxbp)
                 xxro = (gammax2 - 1.0E+00) / gammab * roltmp / pgltmp
     &                * (xxed - btltmp * xxbt - bpltmp * xxbp)

                 btltmp = btltmp + xxbt * dtlocal
                 bpltmp = bpltmp + xxbp * dtlocal

                 if (lmodmap) then
*                   roltmp = roltmp + xxro * dtlocal ! when N-map, this must be done always....?
*                   pgltmp = pgltmp + xxpg * dtlocal
                   if (xxro .GT. 0.0E+00) then
                     roltmp = roltmp + xxro * dtlocal
                     pgltmp = pgltmp + xxpg * dtlocal
                   endif
                 else
                   if (imindenp .EQ. 1) then
                     roltmp = roltmp + xxro * dtlocal
                     pgltmp = pgltmp + xxpg * dtlocal
                     if (roltmp .LT. ro0(j,k,i)) roltmp = ro0(j,k,i)
                     if (pgltmp .LT. pg0(j,k,i)) pgltmp = pg0(j,k,i)
                   else if (imindenp .EQ. 2) then
                     roltmp = roltmp + xxro * dtlocal
                     pgltmp = pgltmp + xxpg * dtlocal
*                     if ((.NOT. lfarmnl) .AND.
*     &                   (roltmp .LT. ro0(j,k,i))) then   ! : only near MNL
                     if (roltmp .LT. ro0(j,k,i)) then
                       aa = roltmp
                       roltmp = ro0(j,k,i)
                       pgltmp = pgltmp * (roltmp/aa)**gammab
                     endif
                   else if (imindenp .EQ. 0) then
                     if (xxro .GT. 0.0E+00) then ! actually, dro/dt etc is negative at MNL region.
*                     if ((xxro .GT. 0.0E+00) .AND. (lfarmnl)) then ! thus, diff. is small......
                       roltmp = roltmp + xxro * dtlocal
                       pgltmp = pgltmp + xxpg * dtlocal
                     endif
                   else
                     roltmp = roltmp + xxro * dtlocal
                     pgltmp = pgltmp + xxpg * dtlocal
                   endif
                 endif

               else if (iclass .EQ. 2) then ! |B| >> 0 and ...
* case (3), v > vc
!                 if ((ichrsurf .EQ. 4) .OR. (ichrsurf .EQ. 7)) then
!                 else if (ichrsurf .EQ. 2) then
!                 else
                   aa = mr1(j,k,i)
                   if (aa .GE. mcrtrn2 * 0.999E+00) then
                     aa = 0.0E+00
                   else
                     aa =(aa - mcrtrn2) / (aa - mrdummy)
                   endif
!                 endif

* increment partially
                 if (aa .LT. 0.0E+00) aa = 0.0E+00
                 if (aa .GT. 1.0E+00) aa = 1.0E+00
                 dtlocal = dt * aa
                 roltmp = ro1(j,k,i) + xxro * dtlocal
!                 if ((ichrsurf .EQ. 2) .OR.
!     &               (ichrsurf .EQ. 4) .OR.
!     &               (ichrsurf .EQ. 7)) then
!                   mrltmp = mr1(j,k,i) + xxmr * dtlocal
!                 else
                   mrltmp = mcrtrn2
!                 endif
                 mtltmp = mt1(j,k,i) + xxmt * dtlocal
                 mpltmp = mp1(j,k,i) + xxmp * dtlocal
                 brltmp = br0(j,k,i)
                 btltmp = bt1(j,k,i) + xxbt * dtlocal
                 bpltmp = bp1(j,k,i) + xxbp * dtlocal
                 edltmp = ed1(j,k,i) + xxed * dtlocal
                 pgltmp =(gammax2 - 1.0E+00)
     &                  *(edltmp
     &                  -(mrltmp**2 + mtltmp**2 + mpltmp**2) * 0.5E+00
     &                  / roltmp
     &                  -(brltmp**2 + btltmp**2 + bpltmp**2) * 0.5E+00)
                 urltmp = mrltmp / roltmp
                 utltmp = mtltmp / roltmp
                 upltmp = mpltmp / roltmp

* rest time-increment
                 dtlocal = dt  - dtlocal
                 roa = roltmp * rr(i)**2
                 mra = mrltmp * rr(i)**2
                 mta = mtltmp * rr(i)**2
                 mpa = mpltmp * rr(i)**2
                 bra = brltmp * rr(i)
                 bta = btltmp * rr(i)
                 bpa = bpltmp * rr(i)
                 eda = edltmp * rr(i)**2
                 pga = pgltmp * rr(i)**2
                 v1 = mra / roa
                 v2 = mta / roa
                 v3 = mpa / roa
                 b1 = bra
                 b2 = bta
                 b3 = bpa
!                 if (lchkchar) then
!                 else
                   call gtmatl(leigemat,
     &                         roa,pga,v1,v2,v3,b1,b2,b3,gammax2)
!                 endif
                 do m = 1, 7 !  pick up raw of matrix with negative eigen value
                   matl74(m,1) = leigemat(m,3) ! V - V_A
                   matl74(m,2) = leigemat(m,5) ! V - V_F
                   matl74(m,3) = leigemat(m,7) ! V - V_S
                   matl74(m,4) = leigemat(m,1) ! V... ususally not used
                   rhsall2(m)  = rhsall(m,j,k) ! mind this comes from using ww2(j,k,i),,,
                 enddo

* BC3a : Mr and T(=Pg/N) are fixed.
                 if (ichrsurf .EQ. 1) then

                   do n = 1, 3 !   condition v // B., for rho, B_theta, B_phi
                     matl33(1,n)= matl74(1,n)
     &                          + matl74(7,n)
     &                             *(pga / roa * dgammax2 * gammax2 ! .... if polytrope is kept in
     &                              -(v1**2 + v2**2 + v3**2)*0.5E+00)
                     matl33(2,n)= matl74(3,n) * mra / b1
     &                          + matl74(5,n)
     &                          + matl74(7,n) *(mra / b1 * v2 + b2)
                     matl33(3,n)= matl74(4,n) * mra / b1
     &                          + matl74(6,n)
     &                          + matl74(7,n) *(mra / b1 * v3 + b3)
                   enddo

                   j2 = j
                   k2 = k
                   iloc = 3
                   call gtans33(matl33,matl74,rhsall2,ansvec3,
     &                          j2,k2,iprc2,ncal2,iloc,
     &                          roa,pga,v1,v2,v3,b1,b2,b3,
     &                          ldummy)
                   if (.NOT. ldummy) nabnorm = nabnorm + 1

                   xxro = ansvec3(1) * drr2
                   xxbt = ansvec3(2) * drr1
                   xxbp = ansvec3(3) * drr1
                   xxmt = xxbt * mrltmp / brltmp
                   xxmp = xxbp * mrltmp / brltmp
                   xxpg = xxro * pgltmp / roltmp * gammax2 ! Polytrope
                   cc = roltmp
                   roltmp = roltmp + xxro * dtlocal
                   btltmp = btltmp + xxbt * dtlocal
                   bpltmp = bpltmp + xxbp * dtlocal
                   mtltmp = mtltmp + xxmt * dtlocal
                   mpltmp = mpltmp + xxmp * dtlocal
                   pgltmp = pgltmp + xxpg * dtlocal  ! normal and orthodox
*                   pgltmp = pgltmp *(roltmp / cc)**gammax2 ! Polytrope

* BC3ab : Mr and (P/N^gamma6) are fixed. If gamma6 is 1, then same as case ichrsurf is 1.
                 else if (ichrsurf .EQ. 6) then

                   do n = 1, 3 !   condition v // B., for rho, B_theta, B_phi
                     matl33(1,n)= matl74(1,n)
     &                          + matl74(7,n)
     &                             *(pga / roa * dgammax2 * gamma6
     &                              -(v1**2 + v2**2 + v3**2)*0.5E+00)
                     matl33(2,n)= matl74(3,n) * mra / b1
     &                          + matl74(5,n)
     &                          + matl74(7,n) *(mra / b1 * v2 + b2)
                     matl33(3,n)= matl74(4,n) * mra / b1
     &                          + matl74(6,n)
     &                          + matl74(7,n) *(mra / b1 * v3 + b3)
                   enddo

                   j2 = j
                   k2 = k
                   iloc = 3
                   call gtans33(matl33,matl74,rhsall2,ansvec3,
     &                          j2,k2,iprc2,ncal2,iloc,
     &                          roa,pga,v1,v2,v3,b1,b2,b3,
     &                          ldummy)
                   if (.NOT. ldummy) nabnorm = nabnorm + 1

                   xxro = ansvec3(1) * drr2
                   xxbt = ansvec3(2) * drr1
                   xxbp = ansvec3(3) * drr1
                   xxmt = xxbt * mrltmp / brltmp
                   xxmp = xxbp * mrltmp / brltmp
                   xxpg = xxro * pgltmp / roltmp * gamma6
*                   cc = roltmp !                           save intermediate value
                   roltmp = roltmp + xxro * dtlocal
                   btltmp = btltmp + xxbt * dtlocal
                   bpltmp = bpltmp + xxbp * dtlocal
                   mtltmp = mtltmp + xxmt * dtlocal
                   mpltmp = mpltmp + xxmp * dtlocal
                   pgltmp = pgltmp + xxpg * dtlocal  !        normal and orthodox
*                   pgltmp = pgltmp *(roltmp / cc)**gamma6  ! enforce poly-relation

* BC3b : Mr and N are fixed, thus Vr is also fixed. T and tangential B will vary.
                 else if (ichrsurf .EQ. 3) then

                   do n = 1, 3 !   condition v // B., for Ed, B_theta, B_phi
                     matl33(1,n)= matl74(7,n)
                     matl33(2,n)= matl74(3,n) * mra / b1
     &                          + matl74(5,n)
                     matl33(3,n)= matl74(4,n) * mra / b1
     &                          + matl74(6,n)
                   enddo

                   j2 = j
                   k2 = k
                   iloc = 3
                   call gtans33(matl33,matl74,rhsall2,ansvec3,
     &                          j2,k2,iprc2,ncal2,iloc,
     &                          roa,pga,v1,v2,v3,b1,b2,b3,
     &                          ldummy)
                   if (.NOT. ldummy) nabnorm = nabnorm + 1

                   xxed = ansvec3(1) * drr2
                   xxbt = ansvec3(2) * drr1
                   xxbp = ansvec3(3) * drr1
                   edltmp = edltmp + dtlocal * xxed
                   btltmp = btltmp + dtlocal * xxbt
                   bpltmp = bpltmp + dtlocal * xxbp
                   mrltmp = mcrtrn2
                   cc = abs(brltmp)
                   if (cc .GT. 1.0E-05) then
                     mtltmp = mrltmp * btltmp / brltmp
                     mpltmp = mrltmp * bpltmp / brltmp
                   else
                     mtltmp = 0.0E+00
                     mpltmp = 0.0E+00
                   endif
                   pgltmp =(gammax2 - 1.0E+00)
     &                    *(edltmp
     &                    -(mrltmp**2 + mtltmp**2 + mpltmp**2)
     &                     /roltmp * 0.5E+00
     &                    -(brltmp**2 + btltmp**2 + bpltmp**2)
     &                    * 0.5E+00)

!                 else if (ichrsurf .EQ. 5) then
!                 else if (ichrsurf .EQ. 2) then
!                 else if (ichrsurf .EQ. 4) then

                 else
                   write(*,*) ' No choice ..... stopped at urwofjdnv'
                   stop
                 endif

               else

                 write(*,*) 'something is wrong !!! at vncvnmcvnmc'
                 stop

               endif ! end of if stagnated or not

             else ! if (.NOT. lchrsurf) then ! fixed

               mrltmp = ur1(j,k,i) * ro1(j,k,i)
               mtltmp = ut1(j,k,i) * ro1(j,k,i)
               mpltmp = up1(j,k,i) * ro1(j,k,i)
*               brltmp = br1(j,k,i)
               btltmp = bt1(j,k,i)
               bpltmp = bp1(j,k,i)
               pgltmp = pg1(j,k,i)
               roltmp = ro1(j,k,i)

               dpgdtmod = 0.0E+00
               drodtmod = 0.0E+00

             endif

* given variations
             pgltmp = pgltmp + dpgdtmod * dt * drr2
             roltmp = roltmp + drodtmod * dt * drr2

* fixed or given boundary variables
             brltmp = br0(j,k,i) ! make sure.....
*             roltmp = ro0(j,k,i)
*             pgltmp = pg0(j,k,i)
* aux. variables
             urltmp = mrltmp / roltmp
             utltmp = mtltmp / roltmp
             upltmp = mpltmp / roltmp

* restore to dim.array
             pg9(j,k,i) = pgltmp ! i = 0
             ro9(j,k,i) = roltmp
             ur9(j,k,i) = urltmp
             ut9(j,k,i) = utltmp
             up9(j,k,i) = upltmp
             br9(j,k,i) = brltmp
             bt9(j,k,i) = btltmp
             bp9(j,k,i) = bpltmp

           enddo ! end of k-loop
           enddo ! end of j-loop

* end of compati. at i = 0
         endif ! if (i .EQ. 0)

* end of compati-or-characteristic methods

* some treatment around/at rot.axis
         if (.NOT. (lpoleave)) then
           if (i .GT. 0) then
             do k = 1, kk
*               ro9( 0,k,i) = ro9(   1,k,i)
*               ro9(jj,k,i) = ro9(jj-1,k,i)
*               pg9( 0,k,i) = pg9(   1,k,i)
*               pg9(jj,k,i) = pg9(jj-1,k,i)
               ur9( 0,k,i) = ur9(   1,k,i)
               ur9(jj,k,i) = ur9(jj-1,k,i)
             enddo
           endif
         endif
         if (i .LE. 0) then ! MIND characteristic method is done for j = 1, jj-1
           do k = 1, kk
             ro9( 0,k,i) = ro9(   1,k,i)
             ro9(jj,k,i) = ro9(jj-1,k,i)
             pg9( 0,k,i) = pg9(   1,k,i)
             pg9(jj,k,i) = pg9(jj-1,k,i)
             ur9( 0,k,i) = ur9(   1,k,i)
             ur9(jj,k,i) = ur9(jj-1,k,i)
             ut9( 0,k,i) = ut9(   1,k,i)
             ut9(jj,k,i) = ut9(jj-1,k,i)
             up9( 0,k,i) = up9(   1,k,i)
             up9(jj,k,i) = up9(jj-1,k,i)
             br9( 0,k,i) = br0(   0,k,i)
             br9(jj,k,i) = br0(jj  ,k,i)
             bt9( 0,k,i) = bt9(   1,k,i)
             bt9(jj,k,i) = bt9(jj-1,k,i)
             bp9( 0,k,i) = bp9(   1,k,i)
             bp9(jj,k,i) = bp9(jj-1,k,i)
           enddo
         endif

       enddo ! end of i-loop
!$omp end parallel do

       if (nabnorm .GT. 0) then
         lstop = .true.
         write(*,*) ' Lstop flag on at 1 : ', nabnorm
         return
       endif

* copy new-step values into the global (common) variables

       nabnorm = 0
!$omp parallel do private(i,k,j)
!$omp&            shared(ii,jj,kk)
!$omp&            shared(ro9,pg9,ur9,ut9,up9,br9,bt9,bp9)
!$omp&            shared(ro2,pg2,ur2,ut2,up2,br2,bt2,bp2)
!$omp&            shared(ro1,pg1,gammax)
!$omp&            reduction(+:nabnorm)
       do i = 0, ii + 1 ! <=== MIND : move surface or not
         do k = 1, kk
         do j = 0, jj
** enforce posititivy.. usually not needed.
           if (lenfposi) then
             if (.NOT. (ro9(j,k,i) .GE. 1.0E-10)) ro9(j,k,i) = 1.0E-10
             if (.NOT. (pg9(j,k,i) .GE. 1.0E-10)) pg9(j,k,i) = 1.0E-10
           endif
**
*           if (ipoly .EQ. 0) then
           if (ichrsurf .EQ. 6) then ! for daily-run.
             ro2(j,k,i) = ro9(j,k,i)
             pg2(j,k,i) = pg9(j,k,i)
           else  !                                 MIND ichrsurf should be 1 or 0
!             if (ipoly .EQ. 1) then
               ro2(j,k,i) = ro9(j,k,i)
               pg2(j,k,i) = pg1(j,k,i)
     &                    *(ro9(j,k,i) / ro1(j,k,i))**gammax           ! poly2
!             else
!               ro2(j,k,i) = ro1(j,k,i)
!     &                    *(pg9(j,k,i) / pg1(j,k,i))**(1.0E+00/gammax) ! poly1
!               pg2(j,k,i) = pg9(j,k,i)
!             endif
           endif
           ur2(j,k,i) = ur9(j,k,i)
           ut2(j,k,i) = ut9(j,k,i)
           up2(j,k,i) = up9(j,k,i)
           br2(j,k,i) = br9(j,k,i)
           bt2(j,k,i) = bt9(j,k,i)
           bp2(j,k,i) = bp9(j,k,i)
*           br2(j,k,i) = 0.0E+00 ! HD !!
*           bt2(j,k,i) = 0.0E+00
*           bp2(j,k,i) = 0.0E+00
         enddo
         enddo
       enddo ! end i-loop
!$omp end parallel do

       if (nabnorm .GT. 0) then
         lstop = .true.
         write(*,*) ' Lstop flag on at 2 : ', nabnorm
         return
       endif

       return
       endsubroutine

*
* -----------------------------------------
*
       pure subroutine mflxvsmg(mcrtrn2,vcrtrn2,br,b0,n0,v0)
       implicit none
       intent(out) ::           mcrtrn2,vcrtrn2
       intent(in)  ::                           br,b0,n0,v0
       real    mcrtrn2, vcrtrn2
       real    br,b0,n0,v0
* constant
       real    nsun0
       real    vcrtrn0
       parameter(nsun0 =    2.00E+08) ! in count/cc
       parameter(vcrtrn0 =  5.00E+05) ! cm/s
* local
       real    brcrtrn
       parameter(brcrtrn = 4.0E+00) ! in gauss
       real    mflratio
       parameter(mflratio = nsun0 * vcrtrn0 / brcrtrn) ! 4 gauss gives limit of 5e5 cm/s x 2e8 /cc

       mcrtrn2 = abs(br) * b0 * mflratio / (n0 * v0)
       vcrtrn2 = mcrtrn2 / (nsun0 / n0)

       return
       end
*
* -----------------------------------------
*
       subroutine gtans33(matl33,matl74,rhsall2,ansvec,
     &                    j,k,iprc,ncal,ilocation,
     &                    roa,pga,v1,v2,v3,b1,b2,b3,lok)
       implicit none
* arguments
       intent(out) ::                        ansvec,lok
       intent(in)  ::    matl33,matl74,rhsall2
       intent(in)  ::    j,k,iprc,ncal,ilocation
       intent(in)  ::    roa,pga,v1,v2,v3,b1,b2,b3
       real    rhsall2(7), matl74(7,4),matl33(3,3),ansvec(3)
       integer j, k, iprc, ncal, ilocation
       real    roa, pga, v1, v2, v3, b1, b2, b3
       logical lok
* local
       real    sumrght(4)
       real    mati(3,3), matl74n(7,4), matl33n(3,3)
       integer n, m
       real    aa, detmat

       lok = .true.
* normalize each row
       do n = 1, 3
         aa = matl33(1,n)**2 + matl33(2,n)**2 + matl33(3,n)**2
         aa = sqrt(aa)
         if (aa .GT. 1.0E-05) then ! <=== previous criterion 1.0e-20
           do m = 1, 3
             matl33n(m,n) = matl33(m,n) / aa
           enddo
           do m = 1, 7
             matl74n(m,n) = matl74(m,n) / aa
           enddo
         else
           write(*,'('' a Norm. = 0 at JK Nrow Ncal Iprc '',5i6)')
     &                  j,k,n,ncal,iprc
           write(*,*) ' location = ' , ilocation
           lok = .false.
           write(*,*) ' values of matrixL are : '
           write(*,*) (matl33(m,n),m=1,3)
           write(*,*) ' values of 8 variables are '
           write(*,*) roa, pga
           write(*,*) v1, v2, v3
           write(*,*) b1, b2, b3
           write(*,*) ' substituted by dummy value of 1'
           do m = 1, 3
             matl33n(m,n) = 1.0E+00
           enddo
         endif
       enddo
       detmat = matl33n(1,1)*(matl33n(2,2)*matl33n(3,3)
     &                       -matl33n(2,3)*matl33n(3,2))
     &        - matl33n(1,2)*(matl33n(2,1)*matl33n(3,3)
     &                       -matl33n(2,3)*matl33n(3,1))
     &        + matl33n(1,3)*(matl33n(2,1)*matl33n(3,2)
     &                       -matl33n(2,2)*matl33n(3,1))
       if (abs(detmat) .GT. 1.0E-05) then
         do n = 1, 3
           sumrght(n) = 0.0E+00
           do m = 1, 7
             sumrght(n) = sumrght(n)
     &                  + matl74n(m,n) * rhsall2(m)
           enddo
         enddo

* inverse of matrix MatL (3 x 3)
         mati(1,1) = (matl33n(2,2) * matl33n(3,3)
     &               -matl33n(2,3) * matl33n(3,2)) / detmat
         mati(2,1) =-(matl33n(2,1) * matl33n(3,3)
     &               -matl33n(2,3) * matl33n(3,1)) / detmat
         mati(3,1) = (matl33n(2,1) * matl33n(3,2)
     &               -matl33n(2,2) * matl33n(3,1)) / detmat
         mati(1,2) =-(matl33n(1,2) * matl33n(3,3)
     &               -matl33n(1,3) * matl33n(3,2)) / detmat
         mati(2,2) = (matl33n(1,1) * matl33n(3,3)
     &               -matl33n(1,3) * matl33n(3,1)) / detmat
         mati(3,2) =-(matl33n(1,1) * matl33n(3,2)
     &               -matl33n(1,2) * matl33n(3,1)) / detmat
         mati(1,3) = (matl33n(1,2) * matl33n(2,3)
     &               -matl33n(1,3) * matl33n(2,2)) / detmat
         mati(2,3) =-(matl33n(1,1) * matl33n(2,3)
     &               -matl33n(1,3) * matl33n(2,1)) / detmat
         mati(3,3) = (matl33n(1,1) * matl33n(2,2)
     &               -matl33n(1,2) * matl33n(2,1)) / detmat
* get ans. of 3 linear eqs.
         do n = 1, 3
           ansvec(n) = 0.0E+00
           do m = 1, 3
             ansvec(n) = ansvec(n) + mati(m,n) * sumrght(m)
           enddo
         enddo
       else ! if det = 0...
         lok = .false.
         write(19,*) ' A det = 0 at ',j,k
         do m = 1, 3
           write(19,*) (matl33(n,m),n=1,3)
         enddo
         write(19,*) roa, pga
         write(19,*) v1,v2,v3
         write(19,*) b1,b2,b3
         write(*,*)  ' A det = 0 at ',j,k
         write(*,*) ' location = ' , ilocation
         do m = 1, 3
           write(*,*) (matl33(n,m),n=1,3)
         enddo
         write(*,*) roa, pga
         write(*,*) v1,v2,v3
         write(*,*) b1,b2,b3
         ansvec(1) =-1.0E+20 ! Tekitou
         ansvec(2) =-1.0E+20
         ansvec(3) =-1.0E+20
       endif

       return
       endsubroutine

*
** -------------------------------------------------------------
*
* calculate values at half mesh point, Roe ave, or other
*
** -------------------------------------------------------------
*
       pure subroutine avehalf(gammax, ! idrct,
     &                    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b, ! R
     &                    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c, ! L
     &                    roa,ura,uta,upa,bra,bta,bpa,pga)
       implicit none
* arguments
       intent(in)  ::          gammax ! ,idrct
       intent(in)  ::     ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b
       intent(in)  ::     ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c
       intent(out) ::     roa,ura,uta,upa,bra,bta,bpa,pga
       real    gammax
!       integer idrct
       real    ro1b,mr1b,mt1b,mp1b,br1b,bt1b,bp1b,ed1b
       real    ro1c,mr1c,mt1c,mp1c,br1c,bt1c,bp1c,ed1c
       real    roa, ura, uta, upa, bra, bta, bpa, pga
* local
!       integer iavemenu
!       parameter(iavemenu = 1) ! 1,2,3,4,5,other= straight Roe, simple, mix, modif Roe, average
       real    aa
       real    sqro1, sqro2, mra, mta, mpa, eda

!       if (iavemenu .EQ. 1) then
         sqro1 = sqrt(ro1c)
         sqro2 = sqrt(ro1b)
         roa = sqro1 * sqro2
         aa  = sqro1 + sqro2
         mra =(sqro2*mr1c + sqro1*mr1b) / aa
         mta =(sqro2*mt1c + sqro1*mt1b) / aa
         mpa =(sqro2*mp1c + sqro1*mp1b) / aa
         bra =(sqro2*br1c + sqro1*br1b) / aa
         bta =(sqro2*bt1c + sqro1*bt1b) / aa
         bpa =(sqro2*bp1c + sqro1*bp1b) / aa
*         if (idrct .EQ. 1)  bra = (br1c + br1b) * 0.5E+00 ! needed ?
*         if (idrct .EQ. 2)  bta = (bt1c + bt1b) * 0.5E+00
*         if (idrct .EQ. 3)  bpa = (bp1c + bp1b) * 0.5E+00
         eda =(sqro2*ed1c + sqro1*ed1b) / aa
         pga = (gammax - 1.0E+00) ! orthodox definition of Pg
     &       * (eda - 0.5E+00 * (mra**2 + mta**2 + mpa**2) / roa
     &              - 0.5E+00 * (bra**2 + bta**2 + bpa**2))
*         pga =(sqro2*pg1c + sqro1*pg1b) / aa ! old ... tricky
*         if (pga .LT. 1.0E-05) pga = 1.0E-05 ! old ... criteria
         if (pga .LT. 1.0E-30) pga = 1.0E-30 ! better ?
         ura = mra / roa
         uta = mta / roa
         upa = mpa / roa
!       else if (iavemenu .EQ. 5) then ! test.
!       else if (iavemenu .EQ. 2) then
!       else if (iavemenu .EQ. 3) then
!       else if (iavemenu .EQ. 4) then
!       else
!       endif

       return
       endsubroutine

*
** -------------------------------------------------------------
*
* Interpolations for cell interface
*
** -------------------------------------------------------------
*
*  make 3rd order MUSCL for non-uniform grids..
*
* If not VL diff. limiter,
*   Case kappa = -1 is rather conventional.
*   Case kappa = 1/3 may provide 3rd order accuracy.
*   Usually, kappa = -1, 0, or 1/3 is preferable, while any float between -1 and 1 can be used.
*
*   If omg (compression parameter) = 1, kappa means nothing.
*
** -------------------------------------------------------------
*
       pure function mscl3rd(dyf,dyc,dx1,dx2) ! locally, this try to obtain gradient at i+1/2
       implicit none
       real          mscl3rd
* arguments
       intent(in)  ::        dyf,dyc,dx1,dx2
       real    dyf, dyc   !  dyc is at interface concerned, dyf is the other one.
       real    dx1, dx2
* const.
       logical lweno ! 3rd-WENO
       parameter(lweno = .false.) ! this will override l3rd-muscl
       logical l3rd
       parameter(l3rd = .true.) ! otherwise be of 2nd order, must be off when ilimiter is not 1
       real    kappa
       parameter(kappa = 0.33333E+00) ! -1 =< kappa < 1
       real    omega ! comp. para
       parameter(omega = (3.0E+00 - kappa)/(1.0E+00 - kappa)) ! omega >= 1.0
*       parameter(omega = 3.0E+00) ! safer choice for ilimiter = 1, or minmod().
*       parameter(omega = 1.0E+00) ! if (ilimiter .ne. 1) in limiter()
       real    smll6
       parameter(smll6 = 1.00E-06)
* local
       real   ans, alp
       real   omg1, omg2, ans1, ans2
* function
*       real   limiter
* interface-statement needed when pure subprog is called in pure subprog
       interface
         pure function limiter(a,b)
           real        limiter
           intent(in) ::       a,b
           real   a, b
         endfunction
       endinterface
*
       if (lweno) then
         alp = dx2 / dx1
         ans1 = dyf * alp
         ans2 = dyc ! face-center..
         omg1 = 1.0E+00 / 3.0E+00 * (ans1**2 + smll6)**2
         omg2 = 2.0E+00 / 3.0E+00 * (ans2**2 + smll6)**2
         mscl3rd = (omg1 * ans1 + omg2 * ans2) / (omg1 + omg2)
     &           * 0.5E+00
       else ! usual MUSCL...
         if (dyf * dyc .GT. 0.0E+00) then

           if (l3rd) then
** modified quasi-3rd.
             alp = dx2 / dx1
             ans = limiter(dyf,omega*dyc) * (1.0E+00 - kappa)* alp**2
     &           + limiter(dyf*omega,dyc) * (1.0E+00 + kappa * alp)
             mscl3rd = ans  / (alp + 1.0E+00) * 0.5E+00
** standard quasi-3rd.
*             alp = dx2 / dx1
*             ans = limiter(dyf*alp,omega*dyc) * (1.0E+00 - kappa) ! most standard..
*     &           + limiter(dyf*alp*omega,dyc) * (1.0E+00 + kappa)
*             mscl3rd = ans * 0.25E+00
** no adjustment...
*             ans = limiter(dyf,omega*dyc) * (1.0E+00 - kappa) ! standard..
*     &           + limiter(dyf*omega,dyc) * (1.0E+00 + kappa)
*             mscl3rd = ans * 0.25E+00
**
           else
** 2nd
             alp = dx2 / dx1
             ans = limiter(dyf*alp,dyc)
             mscl3rd = ans * 0.5E+00
           endif

         else
           mscl3rd = 0.0E+00
         endif
       endif

       return
       endfunction

*
** -------------------------------------------------------------
*
       pure function weno3rd(dyf,dyc,dx1,dx2) ! locally, this try to obtain gradient at i+1/2
       implicit none
       real          weno3rd
* arguments
       intent(in)  ::        dyf,dyc,dx1,dx2
       real    dyf, dyc !  dyc is at interface concerned, dyf is the other one.
       real    dx1, dx2
* local
       real   alp
       real   omg1, omg2, ans1, ans2
* const.
       real    smll6
       parameter(smll6 = 1.00E-06)

         alp = dx2 / dx1
         ans1 = dyf * alp
         ans2 = dyc ! face-center..
         omg1 = 1.0E+00 / 3.0E+00 * (ans1**2 + smll6)**2
         omg2 = 2.0E+00 / 3.0E+00 * (ans2**2 + smll6)**2
         weno3rd = (omg1 * ans1 + omg2 * ans2) / (omg1 + omg2)
     &           * 0.5E+00

       return
       endfunction

*
** --------------------------------------------------------------------
*      non-MUSCL 3rd order interpolation, assuming constant grid size.
** --------------------------------------------------------------------
*
       pure function nonmscl3(dy1,dy2,dy3,dx1,dx2,dx3)
       implicit none
       real          nonmscl3
* arguments
       intent(in) ::          dy1,dy2,dy3,dx1,dx2,dx3
       real    dy1,dy2,dy3 ! near-side center far-side
       real    dx1,dx2,dx3
* local
       logical l3rd2 !       MIND this can be different from "l3rd"
       parameter(l3rd2 = .true.) ! otherwise be of 2nd order, must be off when ilimiter is not 1
       real    ans, alp3, alp1

       if (l3rd2) then
         alp1 = dx2 / dx1
         alp3 = dx2 / dx3
         ans = (dy1*alp1 + 6.0E+00 * dy2 - dy3*alp3) / 12.00E+00 ! modif.
*         ans = (dy1 + 6.0E+00 * dy2 - dy3) / 12.00E+00 ! straight !!
       else !                                 2nd, simplified.
         ans = dy2 * 0.5E+00
       endif

       nonmscl3 = ans

       return
       endfunction

*
** -------------------------------------------------------------
*
* flux limiting function(s)
*     ilimiter = 1 : MINMOD
*     ilimiter = 2 ; MC(monotonized central) limiter
*     ilimiter = 3 ; SB(super bee)
*     ilimiter = 4 ; Van Leer
*     ilimiter = 5 ; Van Albada
*     ilimiter = 6 ; Woodward
*     ilimiter = 7 ; Koren
*     ilimiter = 8 ; Osher
*     ilimiter = 9 ; ospre (Optimum Symmetric Polynomial-Ratio Expression)
*     ilimiter =10 ; GPR (generalized polynimial-ratio)
*     otherwise    ; program be terminated
*
* mind that if ilimiter .NE. 1,
*   1) the compression parameter defined in mscl3rd and mscl3rd0 must be 1.00
*   2) l3rd must be false
*   Otherwise, maybe computation will suffer instability.......uuum.
*
* mind that limiter(a,b) = limiter(r) with r = a/b
*
** -------------------------------------------------------------
*
       pure function limiter(a,b)
       implicit none
       real          limiter
* arguments
       intent(in) ::         a,b
       real   a, b
* local
       real   ans, aa, bb ! , cc, dd
       integer ilimiter
       parameter(ilimiter = 1)
       real   smll6
       parameter(smll6 = 1.00E-06)

       if (a * b .LE. 0.0E+00) then
         ans = 0.0E+00
       else
!         if (ilimiter .EQ. 1) then !                 minmod
           aa = abs(a)
           bb = abs(b)
           ans = min(aa,bb)
!         else
!           write(*,*) 'choise of flux limiter wrong ... '
!           stop
!         endif
       endif
       if (a .LT. 0.0E+00) ans = - ans

       limiter = ans

       return
       endfunction

*
** -------------------------------------------------------------
*
       pure subroutine drtp2num(idrct,
     &                   mrd1,mrd2,mrd3,mtd1,mtd2,mtd3,mpd1,mpd2,mpd3,
     &                   brd1,brd2,brd3,btd1,btd2,btd3,bpd1,bpd2,bpd3,
     &                   m1d1,m1d2,m1d3,m2d1,m2d2,m2d3,m3d1,m3d2,m3d3,
     &                   b1d1,b1d2,b1d3,b2d1,b2d2,b2d3,b3d1,b3d2,b3d3)
       implicit none
* arguments
       intent(in)  ::           idrct
       intent(in)  ::    mrd1,mrd2,mrd3,mtd1,mtd2,mtd3,mpd1,mpd2,mpd3
       intent(in)  ::    brd1,brd2,brd3,btd1,btd2,btd3,bpd1,bpd2,bpd3
       intent(out) ::    m1d1,m1d2,m1d3,m2d1,m2d2,m2d3,m3d1,m3d2,m3d3
       intent(out) ::    b1d1,b1d2,b1d3,b2d1,b2d2,b2d3,b3d1,b3d2,b3d3
       integer idrct
       real    mrd1,mrd2,mrd3,mtd1,mtd2,mtd3,mpd1,mpd2,mpd3
       real    brd1,brd2,brd3,btd1,btd2,btd3,bpd1,bpd2,bpd3
       real    m1d1,m1d2,m1d3,m2d1,m2d2,m2d3,m3d1,m3d2,m3d3
       real    b1d1,b1d2,b1d3,b2d1,b2d2,b2d3,b3d1,b3d2,b3d3

       if (idrct .EQ. 1) then
         m1d1 = mrd1
         m1d2 = mrd2
         m1d3 = mrd3
         m2d1 = mtd1
         m2d2 = mtd2
         m2d3 = mtd3
         m3d1 = mpd1
         m3d2 = mpd2
         m3d3 = mpd3
         b1d1 = brd1
         b1d2 = brd2
         b1d3 = brd3
         b2d1 = btd1
         b2d2 = btd2
         b2d3 = btd3
         b3d1 = bpd1
         b3d2 = bpd2
         b3d3 = bpd3
       else if (idrct .EQ. 2) then
         m1d1 = mtd1
         m1d2 = mtd2
         m1d3 = mtd3
         m2d1 = mpd1
         m2d2 = mpd2
         m2d3 = mpd3
         m3d1 = mrd1
         m3d2 = mrd2
         m3d3 = mrd3
         b1d1 = btd1
         b1d2 = btd2
         b1d3 = btd3
         b2d1 = bpd1
         b2d2 = bpd2
         b2d3 = bpd3
         b3d1 = brd1
         b3d2 = brd2
         b3d3 = brd3
       else ! if (idrct .EQ. 3) then
         m1d1 = mpd1
         m1d2 = mpd2
         m1d3 = mpd3
         m2d1 = mrd1
         m2d2 = mrd2
         m2d3 = mrd3
         m3d1 = mtd1
         m3d2 = mtd2
         m3d3 = mtd3
         b1d1 = bpd1
         b1d2 = bpd2
         b1d3 = bpd3
         b2d1 = brd1
         b2d2 = brd2
         b2d3 = brd3
         b3d1 = btd1
         b3d2 = btd2
         b3d3 = btd3
*       else
*         write(*,*) ' Wrong at eiurei'
       endif

       return
       endsubroutine

*
** -------------------------------------------------------------
*
       pure subroutine vnum2rtp(idrct,
     &                     v1,v2,v3,b1,b2,b3,ura,uta,upa,bra,bta,bpa)
       implicit none
* arguments
       intent(in)  ::           idrct
       intent(in)  ::      v1,v2,v3,b1,b2,b3
       intent(out) ::                        ura,uta,upa,bra,bta,bpa
       integer idrct
       real    v1,v2,v3,b1,b2,b3
       real    ura,uta,upa,bra,bta,bpa

       if (idrct .EQ. 1) then
         ura = v1
         uta = v2
         upa = v3
         bra = b1
         bta = b2
         bpa = b3
       else if (idrct .EQ. 2) then
         ura = v3
         uta = v1
         upa = v2
         bra = b3
         bta = b1
         bpa = b2
       else ! if (idrct .EQ. 3) then
         ura = v2
         uta = v3
         upa = v1
         bra = b2
         bta = b3
         bpa = b1
*       else
*         write(*,*) ' Wrong at eiurei'
       endif

       return
       endsubroutine

*
** -------------------------------------------------------------
*
       pure subroutine vrtp2num(idrct,
     &                     ura,uta,upa,bra,bta,bpa,v1,v2,v3,b1,b2,b3)
       implicit none
* arguments
       intent(in)  ::           idrct
       intent(in)  ::      ura,uta,upa,bra,bta,bpa
       intent(out) ::                              v1,v2,v3,b1,b2,b3
       integer idrct
       real    ura,uta,upa,bra,bta,bpa
       real    v1,v2,v3,b1,b2,b3

       if (idrct .EQ. 1) then
         v1 = ura
         v2 = uta
         v3 = upa
         b1 = bra
         b2 = bta
         b3 = bpa
       else if (idrct .EQ. 2) then
         v1 = uta
         v2 = upa
         v3 = ura
         b1 = bta
         b2 = bpa
         b3 = bra
       else ! if (idrct .EQ. 3) then
         v1 = upa
         v2 = ura
         v3 = uta
         b1 = bpa
         b2 = bra
         b3 = bta
*       else
*         write(*,*) ' Wrong at eiurei'
       endif

       return
       endsubroutine

*
** -------------------------------------------------------------
*
* get matrices of eigen-vectors
*
** -------------------------------------------------------------

*
*--------------------------------------------------------------------
*
* get matrix L only,
*
*--------------------------------------------------------------------
*
       pure subroutine gtmatl(leigemat,
     &                        roa,pga,v1,v2,v3,b1,b2,b3,gammax)
       implicit none
* arguments
       intent(out) ::         leigemat
       intent(in)  ::         roa,pga,v1,v2,v3,b1,b2,b3,gammax
       real    leigemat(7,7)
       real    roa, pga, v1, v2, v3, b1, b2, b3, gammax
* local
       integer n, m
       real    acs, alf1, alf2, alf3, alfall, fst, slw, sqro
       real    alft, blk, gammal, alfall2, blk2, acs2
       real    dgammal, dsqro, dgammax, d2acs2ro, d2acs2sqro
       real    droa, dacs2
       real    aa, bb, cc !, rdummy1
       real    sgnb1
       real    beta2, beta3, ss, ff
*
       real    sqrt2
       parameter(sqrt2 = 0.7071067811865475244E+00)

* initialize
       do n = 1, 7
       do m = 1, 7
         leigemat(m,n) = 0.0E+00
       enddo
       enddo

       if (b1 .LT. 0.0E+00) then
         sgnb1 = -1.0E+00
       else
         sgnb1 =  1.0E+00
       endif

       sqro = sqrt(roa)
       droa  = 1.0E+00 / roa
       dsqro = 1.0E+00 / sqro

       acs2  = gammax * pga * droa ! + (2.0E+00 - gammax) * xx
       acs = sqrt(acs2)
       dacs2 = 1.0E+00 / acs2

       alf1 = b1 * dsqro
       alf2 = b2 * dsqro
       alf3 = b3 * dsqro
       alfall2 = alf1**2 + alf2**2 + alf3**2
       alfall  = sqrt(alfall2)

       dgammax = 1.0E+00 / gammax
       gammal = gammax - 1.0E+00
       dgammal = 1.0E+00 / gammal

       d2acs2ro   = 0.5E+00 * dacs2 * droa
       d2acs2sqro = 0.5E+00 * dacs2 * dsqro
       blk2 = v1**2 + v2**2 + v3**2
       blk = sqrt(blk2)

       aa = (acs2 + alfall2)**2
     &    - 4.0E+00 * acs2 * alf1**2
       if (aa .LE. 0.0E+00) then
         aa = 0.0E+00
       else
         aa = sqrt(aa)
       endif
       fst = sqrt(0.5E+00 * (acs2 + alfall2 + aa))
       aa = 0.5E+00 * (acs2 + alfall2 - aa)
       if (aa .LE. 0.0E+00) then
         aa = 0.0E+00
       else
         aa = sqrt(aa)
       endif
       slw = aa

       alft = sqrt(alf2**2 + alf3**2)
       bb = alft / (abs(alf1) + 1.0E-05)
       if (bb .LT. 1.0E-05) then
         beta2 = sqrt2
         beta3 = sqrt2
         if (b2 .LT. 0.0E+00) beta2 = - beta2
         if (b3 .LT. 0.0E+00) beta3 = - beta3
       else
         beta2 = alf2 / alft
         beta3 = alf3 / alft
       endif

       aa = abs(fst**2 - slw**2)
       if (aa .LT. 1.0E-05) then
         ff = sqrt2
         ss = sqrt2
       else
         bb = abs(acs2 - slw**2)
         cc = abs(fst**2 - acs2)
         ff = sqrt(bb / aa)
         ss = sqrt(cc / aa)
       endif

* eigen values and left & right eigen vector(s) for conservatio variables
*       eigenval(1) = v1
       leigemat(1,1) = 1.0E+00
     &                 - gammal * blk2 * 0.5E+00 * dacs2
       leigemat(2,1) =   gammal * v1 * dacs2
       leigemat(3,1) =   gammal * v2 * dacs2
       leigemat(4,1) =   gammal * v3 * dacs2
       leigemat(5,1) =   gammal * b2 * dacs2
       leigemat(6,1) =   gammal * b3 * dacs2
       leigemat(7,1) = - gammal * dacs2

*       eigenval(2) = v1 + abs(alf1)
       leigemat(1,2) = -(beta3 * v2 - beta2 * v3) * droa *0.5E+00*sgnb1
       leigemat(3,2) =  beta3 * 0.5E+00 * droa * sgnb1
       leigemat(4,2) = -beta2 * 0.5E+00 * droa * sgnb1
       leigemat(5,2) = -beta3 * 0.5E+00 * dsqro
       leigemat(6,2) =  beta2 * 0.5E+00 * dsqro
*       eigenval(3) = v1 - abs(alf1)
       leigemat(1,3) =  (beta3 * v2 - beta2 * v3) * droa *0.5E+00*sgnb1
       leigemat(3,3) = -beta3 * 0.5E+00 * droa * sgnb1
       leigemat(4,3) =  beta2 * 0.5E+00 * droa * sgnb1
       leigemat(5,3) = -beta3 * 0.5E+00 * dsqro
       leigemat(6,3) =  beta2 * 0.5E+00 * dsqro

*       eigenval(4) = v1 + fst
       leigemat(1,4) =(ff * (gammal * blk2 * 0.5E+00 - fst * v1)
     &                +ss * slw * sgnb1 * (beta2 * v2 + beta3 * v3))
     &               * d2acs2ro
       leigemat(2,4) = ff *( fst - gammal * v1)
     &               * d2acs2ro
       leigemat(3,4) = (-ss * slw * beta2 * sgnb1 - ff * gammal * v2)
     &               * d2acs2ro
       leigemat(4,4) = (-ss * slw * beta3 * sgnb1 - ff * gammal * v3)
     &               * d2acs2ro
       leigemat(5,4) = (ss * acs * beta2 - ff * gammal * alf2)
     &               * d2acs2sqro
       leigemat(6,4) = (ss * acs * beta3 - ff * gammal * alf3)
     &               * d2acs2sqro
       leigemat(7,4) = gammal * ff * d2acs2ro
*       eigenval(5) = v1 - fst
       leigemat(1,5) =(ff * (gammal * blk2 * 0.5E+00 + fst * v1)
     &                -ss * slw * sgnb1 * (beta2 * v2 + beta3 * v3))
     &               * d2acs2ro
       leigemat(2,5) = ff *(-fst - gammal * v1)
     &               * d2acs2ro
       leigemat(3,5) = ( ss * slw * beta2 * sgnb1 - ff * gammal * v2)
     &               * d2acs2ro
       leigemat(4,5) = ( ss * slw * beta3 * sgnb1 - ff * gammal * v3)
     &               * d2acs2ro
       leigemat(5,5) = (ss * acs * beta2 - ff * gammal * alf2)
     &               * d2acs2sqro
       leigemat(6,5) = (ss * acs * beta3 - ff * gammal * alf3)
     &               * d2acs2sqro
       leigemat(7,5) = gammal * ff * d2acs2ro

*       eigenval(6) = v1 + slw
       leigemat(1,6) =(ss * (gammal * blk2 * 0.5E+00 - slw * v1)
     &              -ff * fst * sgnb1 * (beta2 * v2 + beta3 * v3))
     &               * d2acs2ro
       leigemat(2,6) = ss *( slw - gammal * v1)
     &               * d2acs2ro
       leigemat(3,6) = (ff * fst * beta2 * sgnb1 - ss * gammal * v2)
     &               * d2acs2ro
       leigemat(4,6) = (ff * fst * beta3 * sgnb1 - ss * gammal * v3)
     &               * d2acs2ro
       leigemat(5,6) =-(ff * acs * beta2 + ss * gammal * alf2)
     &               * d2acs2sqro
       leigemat(6,6) =-(ff * acs * beta3 + ss * gammal * alf3)
     &               * d2acs2sqro
       leigemat(7,6) = gammal * ss * d2acs2ro
*       eigenval(7) = v1 - slw
       leigemat(1,7) =(ss * (gammal * blk2 * 0.5E+00 + slw * v1)
     &                +ff * fst * sgnb1 * (beta2 * v2 + beta3 * v3))
     &               * d2acs2ro
       leigemat(2,7) = ss *(-slw - gammal * v1)
     &               * d2acs2ro
       leigemat(3,7) = (-ff * fst * beta2 * sgnb1 - ss * gammal * v2)
     &               * d2acs2ro
       leigemat(4,7) = (-ff * fst * beta3 * sgnb1 - ss * gammal * v3)
     &               * d2acs2ro
       leigemat(5,7) =-(ff * acs * beta2 + ss * gammal * alf2)
     &               * d2acs2sqro
       leigemat(6,7) =-(ff * acs * beta3 + ss * gammal * alf3)
     &               * d2acs2sqro
       leigemat(7,7) = gammal * ss * d2acs2ro

* make continurous eigen vector....
       if ((beta2 .LT. 0.0E+00) .OR.
     &     (abs(beta2) .LT. 1.0E-05) .AND. (beta3 .LT. 0.0E+00)) then
         aa = acs2 - alfall2
         if (aa .GT. 0.0E+00) then
           do n = 1, 7
             leigemat(n,6) = - leigemat(n,6)
             leigemat(n,7) = - leigemat(n,7)
           enddo
         else
           do n = 1, 7
             leigemat(n,4) = - leigemat(n,4)
             leigemat(n,5) = - leigemat(n,5)
           enddo
         endif
       endif

       return
       endsubroutine

*
*--------------------------------------------------------------------
*
* END of This file
*
*--------------------------------------------------------------------
*
