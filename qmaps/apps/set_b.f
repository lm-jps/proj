c#######################################################################
      subroutine set_b (nr,nt,np,br,bt,bp,r,t,p,b)
c
c-----------------------------------------------------------------------
c
c ****** Set up magnetic field vector structure B based on the
c ****** magnetic field components passed in.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use vars
      use params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nr,nt,np
      real(r_typ), dimension(nr,nt,np), target :: br,bt,bp
      real(r_typ), dimension(nr), target :: r
      real(r_typ), dimension(nt), target :: t
      real(r_typ), dimension(np), target :: p
      type(vec) :: b
c
c-----------------------------------------------------------------------
c
c ****** Load the magnetic field structure components.
c
      b%r%ndim=3
      b%r%dims(1)=nr
      b%r%dims(2)=nt
      b%r%dims(3)=np
      b%r%scale=.true.
      b%r%hdf32=.false.
      b%r%scales(1)%f=>r
      b%r%scales(2)%f=>t
      b%r%scales(3)%f=>p
      b%r%f=>br
c
      b%t%ndim=3
      b%t%dims(1)=nr
      b%t%dims(2)=nt
      b%t%dims(3)=np
      b%t%scale=.true.
      b%t%hdf32=.false.
      b%t%scales(1)%f=>r
      b%t%scales(2)%f=>t
      b%t%scales(3)%f=>p
      b%t%f=>bt
c
      b%p%ndim=3
      b%p%dims(1)=nr
      b%p%dims(2)=nt
      b%p%dims(3)=np
      b%p%scale=.true.
      b%p%hdf32=.false.
      b%p%scales(1)%f=>r
      b%p%scales(2)%f=>t
      b%p%scales(3)%f=>p
      b%p%f=>bp
c
c ****** Set the type of magnetic field files.
c
      call set_btype (b)
c
c ****** Build the inverse interpolation tables.
c
      call build_inverse_tables (b%r,b%inv(1))
      call build_inverse_tables (b%t,b%inv(2))
      call build_inverse_tables (b%p,b%inv(3))
c
c ****** If cubic spline interpolation was requested, get the
c ****** spline coefficients.
c
      if (cubic) then
        b%cubic=.true.
        if (verbose) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//
     &                ' for Br ...'
        end if
        call compute_spline_3d (b%r%dims(1),b%r%dims(2),b%r%dims(3),
     &                          b%r%scales(1)%f,
     &                          b%r%scales(2)%f,
     &                          b%r%scales(3)%f,
     &                          b%r%f,b%spl%r)
        if (verbose) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//
     &                ' for Bt ...'
        end if
        call compute_spline_3d (b%t%dims(1),b%t%dims(2),b%t%dims(3),
     &                          b%t%scales(1)%f,
     &                          b%t%scales(2)%f,
     &                          b%t%scales(3)%f,
     &                          b%t%f,b%spl%t)
        if (verbose) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//
     &                ' for Bp ...'
        end if
        call compute_spline_3d (b%p%dims(1),b%p%dims(2),b%p%dims(3),
     &                          b%p%scales(1)%f,
     &                          b%p%scales(2)%f,
     &                          b%p%scales(3)%f,
     &                          b%p%f,b%spl%p)
      else
        b%cubic=.false.
      end if
c
      return
      end
c#######################################################################
      subroutine set_btype (b)
c
c-----------------------------------------------------------------------
c
c ****** Store the primary (r,t,p) scales in structure B.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
c
c-----------------------------------------------------------------------
c
c ****** Assume that all components of B are on the same mesh.
c
      b%nrs=b%r%dims(1)
      b%nts=b%r%dims(2)
      b%nps=b%r%dims(3)
c
      b%rs=>b%r%scales(1)%f
      b%ts=>b%r%scales(2)%f
      b%ps=>b%r%scales(3)%f
c
      b%lim0(1)=b%rs(1)
      b%lim1(1)=b%rs(b%nrs)
      b%lim0(2)=b%ts(1)
      b%lim1(2)=b%ts(b%nts)
      b%lim0(3)=b%ps(1)
      b%lim1(3)=b%ps(b%nps)
c
c ****** Build the inverse interpolation tables for the
c ****** main mesh.
c
      b%rs_invtab%n=b%nrs
      allocate (b%rs_invtab%f(b%rs_invtab%n))
      call getinv (b%rs,b%nrs,b%rs_invtab)
c
      b%ts_invtab%n=b%nts
      allocate (b%ts_invtab%f(b%ts_invtab%n))
      call getinv (b%ts,b%nts,b%ts_invtab)
c
      b%ps_invtab%n=b%nps
      allocate (b%ps_invtab%f(b%ps_invtab%n))
      call getinv (b%ps,b%nps,b%ps_invtab)
c
c ****** Compute the mesh cell dimensions on the main mesh.
c
c ****** These are used in setting the field line integration
c ****** step size.
c
      allocate (b%drs(b%nrs))
      allocate (b%dts(b%nts))
      allocate (b%dps(b%nps))
      allocate (b%sts(b%nts))
c
      call get_dx (b%nrs,b%rs,b%drs)
      call get_dx (b%nts,b%ts,b%dts)
      call get_dx (b%nps,b%ps,b%dps)
      b%sts=sin(b%ts)
      b%sts(    1)=max(b%sts(    1),sin(b%dts(    1)))
      b%sts(b%nts)=max(b%sts(b%nts),sin(b%dts(b%nts)))
c
      return
      end
c#######################################################################
      subroutine set_slice_coordinates (r_q,t_q,p_q)
c
c-----------------------------------------------------------------------
c
c ****** Set up the coordinates that define the slice in
c ****** the 3D volume.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(npss,ntss,1), target :: r_q,t_q,p_q
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      logical, external :: same_structure_sds
c
c-----------------------------------------------------------------------
c
c ****** Set the coordinate names.
c
      if (slice_coords_are_xyz) then
        slice_coord_name=(/'x','y','z'/)
      else
        slice_coord_name=(/'r','t','p'/)
      end if
c
c ****** Set the coordinates of the slice.
c
      slice_c1%ndim=2
      slice_c1%dims(1)=npss
      slice_c1%dims(2)=ntss
      slice_c1%dims(3)=1
      slice_c1%scale=.false.
      slice_c1%hdf32=.false.
      slice_c1%f=>r_q
c
      slice_c2%ndim=2
      slice_c2%dims(1)=npss
      slice_c2%dims(2)=ntss
      slice_c2%dims(3)=1
      slice_c2%scale=.false.
      slice_c2%hdf32=.false.
      slice_c2%f=>t_q
c
      slice_c3%ndim=2
      slice_c3%dims(1)=npss
      slice_c3%dims(2)=ntss
      slice_c3%dims(3)=1
      slice_c3%scale=.false.
      slice_c3%hdf32=.false.
      slice_c3%f=>p_q
c
      return
      end
