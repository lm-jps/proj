      subroutine mapfl_wrapper(
     1 bp, bt, br, p, t, r, np, nt, nr,
     2 p_q, t_q, r_q, np_q, nt_q, r_ch,
     3 do_chmap, vb, do_cubic,
     4 qmap, chmap)

!==============================================================!

c      use number_types
      use ident
      use params
      use types
      use files
      use mesh
      use field
      use vars
      use field_line_params
      use step_size_stats
      use debug

c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nr,nt,np,np_q,nt_q,do_chmap
      real(r_typ), dimension(nr,nt,np), target :: br,bt,bp
      real(r_typ), dimension(nr), target :: r
      real(r_typ), dimension(nt), target :: t
      real(r_typ), dimension(np), target :: p
      real(r_typ) :: r_ch
      real(r_typ), dimension(np_q,nt_q), target :: p_q,t_q,r_q
      real(r_typ), dimension(np_q,nt_q), target :: qmap,chmap
      integer :: vb, do_cubic

c-----------------------------------------------------------------------

      integer :: ierr
      logical :: trace_fwd,trace_bwd
      logical :: compute_q_on_slice
      logical :: compute_ch_map
      real(r_typ) :: ch_map_r
      integer :: rtest,ttest,ptest
      character(len=80) :: fmt
c
      ttest = 363 + 1
      ptest = 325 + 1
      fmt = "(3F20.16)"
c-----------------------------------------------------------------------

      if (vb .eq. 0) then
        verbose = .false.
      else
        verbose = .true.
      end if
      debug_level = 0_r_typ
      domain_r_min = 0.0_r_typ
      domain_r_max = 1000000.0_r_typ
      if (do_cubic .eq. 0) then
        cubic = .false.
      else
        cubic = .true.
      end if
      ds%variable = .true.
      ds%over_rc = 0.001_r_typ
      set_ds_automatically = .false.
      ds%min = 0.0002_r_typ
      ds%max = 0.002_r_typ
      dsmult = 1._r_typ
      ds%limit_by_local_mesh = .true.
      ds%local_mesh_factor = 1.
      ds%lmax = 100._r_typ
      trace_fwd = .false.
      trace_bwd = .false.
      new_t_mesh = .true.
      mesh_file_t = ' '
      ntss = nt_q
      t0 = 0._r_typ
      t1 = 0._r_typ
      new_p_mesh = .true.
      mesh_file_p = ' '
      npss = np_q
      p0 = 0._r_typ
      p1 = 0._r_typ
      compute_q_on_slice = .true.
      slice_coords_are_xyz = .false.
      q_increment_h = 0.0001_r_typ
      ch_map_r = r_ch
      if (do_chmap .eq. 0) then
        compute_ch_map = .false.
      else
        compute_ch_map = .true.
      end if
c
      if (verbose) then
        write (*,*)
        write (*,*) '### ',cname,' Version ',cvers,' of ',cdate,'.'
      end if
c
c ****** Some test printing.
      write(*,'(A,F19.16)') 'ph=', p_q(ptest,ttest)
      write(*,'(A,F19.16)') 'th=', t_q(ptest,ttest)
      write(*,'(A,F19.16)') 'r=', r_q(ptest,ttest)
c
c ****** Set the field line integration parameters.
c
      ds%max_increase_factor=max_increase_factor
      ds%max_decrease_factor=max_decrease_factor
      ds%predictor_min_clip_fraction=predictor_min_clip_fraction
      ds%short_fl_min_points=short_fl_min_points
      ds%short_fl_max_tries=short_fl_max_tries
      ds%short_fl_shrink_factor=short_fl_shrink_factor

c
c ****** Read the magnetic field.
c
      call set_b (nr,nt,np,br,bt,bp,r,t,p,b)

c-----------------------------------------------------------------------

      if (.false.) then
      write(*,*) br(1,1,1), br(2,1,1), br(3,1,1)
      write(*,*) br(1,1,1), br(1,2,1), br(1,3,1)
      write(*,*) br(1,1,1), br(1,1,2), br(1,1,3)
      write(*,*) b%r%f(1,1,1), b%r%f(2,1,1), b%r%f(3,1,1)
      write(*,*) b%r%f(1,1,1), b%r%f(1,2,1), b%r%f(1,3,1)
      write(*,*) b%r%f(1,1,1), b%r%f(1,1,2), b%r%f(1,1,3)
      write(*,*) '--------------'
      write(*,*) p_q(1,1), p_q(1,2), p_q(1,3)
      write(*,*) p_q(1,1), p_q(2,1), p_q(3,1)
      write(*,*) t_q(1,1), t_q(1,2), t_q(1,3)
      write(*,*) t_q(1,1), t_q(2,1), t_q(3,1)
      write(*,*) r_q(1,1), r_q(1,2), r_q(1,3)
      write(*,*) r_q(1,1), r_q(2,1), r_q(3,1)
      endif
      if (.true.) then
      write(*,fmt) b%r%f(1,1,1), b%r%f(2,1,1), b%r%f(3,1,1)
      write(*,fmt) b%r%f(1,1,1), b%r%f(1,2,1), b%r%f(1,3,1)
      write(*,fmt) b%r%f(1,1,1), b%r%f(1,1,2), b%r%f(1,1,3)
      endif

c-----------------------------------------------------------------------

c
c ****** Set the radial domain limits to those specified.
c
      b%lim0(1)=max(b%lim0(1),domain_r_min)
      b%lim1(1)=min(b%lim1(1),domain_r_max)
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Domain limits:'
        write (*,*) 'Lower boundary value: ',b%lim0(1)
        write (*,*) 'Upper boundary value: ',b%lim1(1)
      end if
c
c ****** Make the new r, t, and p meshes.
c
      call make_new_meshes (b)
c
c ****** Set the default step size.
c
      call set_ds (b)
c
c ****** Set the flag to gather step size statistics.
c
      gather_stats=verbose
c
c ****** Trace the field lines forward, if requested.
c
      if (trace_fwd) call map_forward
c
c ****** Trace the field lines backward, if requested.
c
      if (trace_bwd) call map_backward
c
c ****** Determine Q on a slice, if requested.
c
      if (compute_q_on_slice) then
        call set_slice_coordinates(r_q,t_q,p_q)
c ****** the returned qmap is already rotated to C convention
        call get_q_on_slice(qmap)
c        no need to deallocate mesh as the pointer is passed in
c        call deallocate_slice_coordinates
      end if
c
c ****** Compute a coronal hole map, if requested.
c
      if (compute_ch_map) then
c ****** the returned chmap is already rotated to C convention
        call get_ch_map (ch_map_r,chmap)
      end if
c
      if (verbose) then
        stat_ds_avg=0.
        if (stat_n.ne.0) stat_ds_avg=stat_ds_sum/stat_n
        write (*,*)
        write (*,*) '### Field line integration step size statistics:'
        write (*,*) 'Number of field line segments = ',stat_n
        write (*,*) 'Minimum step size used = ',stat_ds_min
        write (*,*) 'Maximum step size used = ',stat_ds_max
        write (*,*) 'Average step size used = ',stat_ds_avg
      end if

c-----------------------------------------------------------------------

      end subroutine mapfl_wrapper
