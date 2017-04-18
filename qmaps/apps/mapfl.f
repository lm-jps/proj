c
c-----------------------------------------------------------------------
c
c ****** Get the field line mapping for MAS code runs, computing
c ****** the structural quantity Q.
c
c-----------------------------------------------------------------------
c
c ****** Updates and bug fixes:
c
c        08/30/2005, ZM, Version 1.06:
c
c         - Converted MAPFL v1.05 into a tool.
c         - Improved the tracing accuracy of short field lines
c           (e.g., those near the neutral line).
c         - Improved the field line integrator.
c         - Added the ability to perform a 3D mapping of the
c           field lines.
c
c        06/16/2006, ZM, Version 1.07:
c
c         - Fixed the effect of rondoff errors in generating the
c           meshes for the calculation which caused initial points
c           to lie outside the computation domain.
c
c        06/21/2006, ZM, Version 1.08:
c
c         - Allowed the ability of specifying planes for the
c           field line mapping region and the 3D mapping.
c           This extends the capability of the tool to get
c           Q on 2D slices.
c
c        07/07/2006, ZM, Version 1.09:
c
c         - Added the ability to do new POT3D files,
c           new-old POT3D files, and old POT3D files.
c
c        08/22/2006, ZM, Version 1.10:
c
c         - Corrected a bug in the computation of Q for field
c           lines that wrap periodically at the join between
c           phi=0 and phi=2*pi.
c         - Corrected the computation of Q on open field lines
c           to take into account the difference in the radii
c           of the initial and final field line footpoints.
c         - Added the ability to use cubic spline interpolation.
c           This increases the storage requirements substantially,
c           and makes the computation significantly slower.
c
c        09/17/2006, ZM, Version 1.11:
c
c         - Reverted to the Cartesian field line integrator.
c           The accuracy of the spherical integrator was being
c           cast into doubt, though it was not proven to be
c           bad.
c         - Improved the accuracy with which the final point
c           (i.e., the end point at the r boundaries) was
c           being clipped to the boundary.
c         - Changed the computation of Q to be staggered with
c           respect to the mapping quantities.
c         - Fixed the backward mapping routine to compute Q.
c
c        09/29/2006, ZM, Version 1.12:
c
c         - Improved the field line integrator to use a variable
c           step size.  It is now possible to select either a uniform
c           step size for the field line integration, as before,
c           or to use a variable step size based on the local radius
c           of curvature of the field line, and the local mesh size
c           of the B files along the field line trace.
c         - This has changed the format of the input file.
c
c        01/28/2007, ZM, Version 1.13:
c
c         - Added the ability to compute the mapping on 2D slices
c           through the 3D volume.
c         - This has changed the format of the input file.
c         - Cleaned up some formatting of the output.
c
c        02/19/2008, ZM, Version 1.14:
c
c         - Changed the default behavior that terminated the code
c           when more than 100 bad field line traces were found.
c           By default this is now disabled, but can be put back
c           by setting the variable MAX_BAD_FIELDLINES.
c
c        04/06/2009, ZM, Version 1.15:
c
c         - Added the ability to use an analytic function to define
c           the magnetic field.
c         - Corrected bugs involving how the expansion factor and Q
c           were computed for the backward trace.
c         - Performed a cosmetic clean-up of the code.
c
c        04/21/2009, ZM, Version 1.16:
c
c         - Added the ability to use Slava's new formulation
c           to compute Q directly on slices within the domain
c           by tracing a bundle of field lines from points in
c           the domain forward and backward to the boundaries.
c         - Cleaned up the way mapping along "slices" in the 3D
c           volume is implemented.  These "slices" can now be
c           lines (i.e., 1D files), 2D slices, or 3D volumes.
c           These are all defined by reading in rectilinear
c           files (1D, 2D, or 3D) that contain the (r,t,p) or
c           (x,y,z) starting coordinates for the mapping.
c         - Added a check to make sure that the input file
c           has the correct number of lines.
c
c        02/18/2010, ZM, Version 1.17:
c
c         - Made the ability to trace from a slice more flexible.
c           It is now possible to map from a slice along the
c           forward and backward directions separately.  It is
c           also possible to specify this direction to be either
c           along B or along the direction of increasing radius.
c         - Added the ability to directly compute coronal hole
c           maps at a given radius.  These are coded by magnetic
c           field polarity.
c
c        04/12/2010, ZM, Version 1.18:
c
c         - Added the ability of specifying the r, t, and p meshes
c           to be used for the calculation.  These can be specified
c           as 1D HDF files or as uniform meshes.
c         - Allowed the phi coordinate to be outside the [0,2*pi]
c           range in the input file.  It is now properly wrapped
c           into the main interval during the calculation.
c
c        04/26/2010, ZM, Version 1.19:
c
c         - Added the ability to compute 3D coronal hole maps.
c           These are useful to compute 3D open/closed field regions
c           for use, perhaps, in developing heating masks.
c         - Removed the writing of warning messages about field
c           lines that did not reach the boundaries in the coronal
c           hole map traces.  This was not really necessary since
c           such field lines are already flagged (with values of
c           "-2") in the coronal hole maps.
c
c        02/07/2011, ZM, Version 1.20:
c
c         - Added a multi-threading capability using OpenMP.
c           This version runs in parallel.  This required a
c           restructuring of the code to improve the modularity,
c           to make the code thread-safe, and to improve the
c           parallel performance.
c         - It is necessary to use "thread-safe" programming
c           in the parallel sections from now on.
c
c        05/18/2011, ZM, Version 1.21:
c
c         - Corrected a bug in the periodic wrapping of the phi
c           coordinate when it was outside the range [0,2*pi].
c           It was not being wrapped to the main periodic interval
c           in the GETB routine.
c
c        08/17/2012, ZM, Version 1.22:
c
c         - Fixed a minor bug that was discovered by compiling
c           with GFORTRAN.
c
c        04/29/2013, ZM, Version 1.23:
c
c         - Changed the interpretation of the value of DEBUG_LEVEL
c           for debugging output.  This was done to optimize
c           Slava's output of field line traces.  The functionality
c           of the program was not otherwise modified.
c
c        05/13/2013, ZM, Version 1.24:
c
c         - Added the ability to write field line traces to HDF
c           files.  This option can be selected when tracing from
c           a "slice".  It is not intended to be used when doing
c           very high resolution mapping computations, since it
c           can produce a large amount of data.  To get the HDF
c           files, set DEBUG_LEVEL.ge.2.  The field lines traces
c           for the forward and/or backward trace will be written
c           to HDF files named:
c
c             field_line_traces_forward_rtp.hdf
c             field_line_traces_forward_xyz.hdf
c             field_line_traces_backward_rtp.hdf
c             field_line_traces_backward_xyz.hdf
c
c        09/25/2015, ZM, Version 1.25:
c
c         - Changed the way short field lines are treated.
c           Because of strange behavior noticed by Slava in a
c           certain case, I changed the way the step size for
c           "short" field lines was controlled.  It turned out
c           that a "short" field line could become a very long
c           field line when retraced!  Previously, these field
c           lines were retraced with the miniumum step size, which
c           led to a very long field line with lots of points.
c           I relaxed this requirement on the step size once the
c           number of points in the retraced field line exceeded
c           the minimum number of points.
c
c        01/13/2016, ZM, Version 1.25_SU:
c
c         - Simplified version for use in calculating a limited set
c           of quantities in the SDO data pipeline.
c
c        02/09/2016, ZM, Version 1.26_SU:
c
c         - Fixed the check that adds a phi point to the B arrays
c           for cases when all the B components are on the same
c           mesh.  The code now checks if the point at phi=2*pi is
c           already present, and only adds a point if it is not.
c         - Fixed the storage of field line traces in the field line
c           buffers.  Previously, when "short" field lines were
c           detected, the buffers were not reinitialized when these
c           short field lines were retraced with a smaller step size,
c           leading to incorrect saved field line traces.
c         - Fixed the computation of Q on slices for field line
c           footpoints that approach the poles in routine GETQ.
c           The previous differencing was not accurate for field
c           line footpoints near the poles.  The new scheme switches
c           to Cartesian basis vectors in a small neighborhood of
c           the poles.
c
c-----------------------------------------------------------------------
c
c#######################################################################
c
c ****** This tool uses modules from ZM's tools library.
c
c#######################################################################
      module ident
c
      character(*), parameter :: cname='MAPFL'
      character(*), parameter :: cvers='1.26_SU'
      character(*), parameter :: cdate='02/09/2016'
c
      end module
c#######################################################################
      module debug
c
      implicit none
c
c ****** Debugging level.
c
      integer :: debug_level=0
c
      end module
c#######################################################################
      module constants
c
c-----------------------------------------------------------------------
c ****** Constants.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
      real(r_typ), parameter :: halfpi=.5_r_typ*pi
      real(r_typ), parameter :: twopi=2._r_typ*pi
c
      end module
c#######################################################################
      module types
c
c-----------------------------------------------------------------------
c ****** Definition of data structures.
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use invint_def
      use spline_def
c
      implicit none
c
c ****** Maximum number of dimensions.
c
      integer, parameter, private :: ndim_max=3
c
c ****** Inverse interpolation table structure for a vector field.
c
      type :: vtab
        type(itab), dimension(ndim_max) :: c
      end type
c
c ****** Vector spline structure.
c
      type :: vspl3d
        type(spl3d) :: r
        type(spl3d) :: t
        type(spl3d) :: p
      end type
c
c ****** Magnetic field vector structure.
c
      type :: vec
        type(sds) :: r
        type(sds) :: t
        type(sds) :: p
        integer :: nrs
        integer :: nts
        integer :: nps
        real(r_typ), dimension(:), pointer :: rs
        real(r_typ), dimension(:), pointer :: ts
        real(r_typ), dimension(:), pointer :: ps
        real(r_typ) :: lim0(ndim_max)
        real(r_typ) :: lim1(ndim_max)
        type(vtab), dimension(ndim_max) :: inv
        real(r_typ), dimension(:), pointer :: drs
        real(r_typ), dimension(:), pointer :: dts
        real(r_typ), dimension(:), pointer :: dps
        real(r_typ), dimension(:), pointer :: sts
        type(itab) :: rs_invtab
        type(itab) :: ts_invtab
        type(itab) :: ps_invtab
        logical :: cubic
        type(vspl3d) :: spl
      end type
c
c ****** Data structure to hold file names of a vector component.
c
      type :: vfile
        character(512) :: r
        character(512) :: t
        character(512) :: p
      end type
c
c ****** Initial size for the field line trace buffer.
c
      integer, parameter, private :: fl_buffer_size=1000
c
c ****** Trajectory structure definition.
c
      type :: traj
        integer :: ndim
        integer :: initial_size=fl_buffer_size
        integer :: size
        integer :: npts
        type(rp1d), dimension(:), pointer :: x
      end type
c
c ****** Dual representation Cartesian and spherical position vector.
c
      type :: csvec
        real(r_typ), dimension(3) :: c
        real(r_typ), dimension(3) :: s
      end type
c
c ****** "Inside domain" structure.
c
      type :: inout
        logical :: domain
        logical :: r0
        logical :: r1
        logical :: t0
        logical :: t1
        logical :: p0
        logical :: p1
        logical :: r
        logical :: t
        logical :: p
      end type
c
c ****** Field line integration parameters.
c
      type :: flparam
        logical :: variable
        integer :: direction
        logical :: direction_is_along_b
        real(r_typ) :: min
        real(r_typ) :: max
        real(r_typ) :: over_rc
        real(r_typ) :: lmax
        logical :: limit_by_local_mesh
        real(r_typ) :: local_mesh_factor
        real(r_typ) :: max_increase_factor
        real(r_typ) :: max_decrease_factor
        integer :: short_fl_min_points
        integer :: short_fl_max_tries
        real(r_typ) :: short_fl_shrink_factor
        real(r_typ) :: predictor_min_clip_fraction
      end type
c
      end module
c#######################################################################
      module mesh
c
c-----------------------------------------------------------------------
c ****** Meshes.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ****** Secondary (t,p) meshes.
c
      integer :: ntss
      integer :: npss
c
      real(r_typ), dimension(:), pointer :: tss
      real(r_typ), dimension(:), pointer :: pss
c
      end module
c#######################################################################
      module field
c
c-----------------------------------------------------------------------
c ****** Magnetic field storage.
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
      implicit none
c
c ****** Structure that holds the magnetic field.
c
      type(vec) :: b
c
      end module
c#######################################################################
      module vars
c
c-----------------------------------------------------------------------
c ****** Input variables, switches, etc.
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use types
c
      implicit none
c
c ****** Spherical geometry domain limits.
c
      real(r_typ) :: domain_r_min=1._r_typ
      real(r_typ) :: domain_r_max=30._r_typ
c
c ****** New output mesh flags.
c
      logical :: new_t_mesh=.false.
      logical :: new_p_mesh=.false.
c
c ****** New output mesh limits.
c
      real(r_typ) :: r0=0.
      real(r_typ) :: r1=0.
      real(r_typ) :: t0=0.
      real(r_typ) :: t1=0.
      real(r_typ) :: p0=0.
      real(r_typ) :: p1=0.
c
c ****** Flag to use tri-cubic interpolation (when .TRUE.) or
c ****** simple linear interpolation (when .FALSE.) to
c ****** interpolate B between mesh points.
c
      logical :: cubic
c
c ****** Field line integration.
c
      type(flparam) :: ds
      logical :: set_ds_automatically
      real(r_typ) :: dsmult
c
c ****** Parameters for the slice mapping.
c
      logical :: slice_coords_are_xyz
c
c ****** Names of the slice coodrinates.
c
      character, dimension(3) :: slice_coord_name
c
c ****** Structures that hold the slice coordinates.
c
      type (sds) :: slice_c1,slice_c2,slice_c3
c
c ****** Increment to compute Q directly on a slice.
c
      real(r_typ) :: q_increment_h
c
c ****** Flag to use 32-bit HDF output files.
c
      logical :: hdf32=.true.
c
      end module
c#######################################################################
      module diags
c
c-----------------------------------------------------------------------
c ****** Variables that control diagnostic output.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ****** Number of iterations between prints of diagnostics
c ****** during execution.
c
      integer :: diagnostic_interval=1000
c
c ****** File sequence interval for diagnostic files.
c
      integer :: diag_seq=0
c
      end module
c#######################################################################
      module files
c
c-----------------------------------------------------------------------
c ****** File names.
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
      implicit none
c
      type(vfile) :: bfile
c
      character(512) :: rffile,tffile,pffile,qffile
      character(512) :: rbfile,tbfile,pbfile,qbfile
c
c ****** File names for the t and p meshes.
c
      character(512) :: mesh_file_t
      character(512) :: mesh_file_p
c
      type(vfile) :: slice_input_file
c
      character(512) :: slice_q_output_file
c
c ****** File name for the output coronal hole map.
c
      character(512) :: ch_map_output_file
c
      end module
c#######################################################################
      module field_line_params
c
c-----------------------------------------------------------------------
c ****** Parameters that control the field line integration.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Parameters that control variable step-size tracing.
c-----------------------------------------------------------------------
c
c ****** These factors control how much the field line integration
c ****** step size can change from one step to another for the
c ****** case when a variable step size is being used.
c
c ****** MAX_INCREASE_FACTOR should be greater than 1, and
c ****** MAX_DECREASE_FACTOR should be less than 1.

      real(r_typ), parameter :: max_increase_factor=1.5_r_typ
      real(r_typ), parameter :: max_decrease_factor=.1_r_typ
c
c-----------------------------------------------------------------------
c ****** Parameters that control tracing of short field lines.
c-----------------------------------------------------------------------
c
c ****** If a field line trace has a smaller number of points
c ****** than SHORT_FL_MIN_POINTS, the integration step size
c ****** is decreased by the factor SHORT_FL_SHRINK_FACTOR,
c ****** and it is retraced, up to a maximum number of tries
c ****** equal to SHORT_FL_MAX_TRIES.
c
c ****** SHORT_FL_SHRINK_FACTOR should be less than 1.
c
      integer, parameter :: short_fl_min_points=10
      integer, parameter :: short_fl_max_tries=5
      real(r_typ), parameter :: short_fl_shrink_factor=.1_r_typ
c
c-----------------------------------------------------------------------
c ****** Parameters that control clipping to boundaries.
c-----------------------------------------------------------------------
c
c ****** The factor PREDICTOR_MIN_CLIP_FRACTION determines when to
c ****** clip a trace to the radial boundary in the predictor.
c ****** When the normalized distance to the r boundary (as a
c ****** fraction of the current step size) is less than
c ****** PREDICTOR_MIN_CLIP_FRACTION, the field line is clipped
c ****** to the boundary in the predictor without doing a
c ****** corrector step.  This number should be between 0 and 1.
c
      real(r_typ), parameter :: predictor_min_clip_fraction=.1_r_typ
c
c-----------------------------------------------------------------------
c ****** Maximum number of "bad" field line traces after
c ****** which to terminate.  Set to -1 to disable the termination.
c-----------------------------------------------------------------------
c
      integer, parameter :: max_bad_fieldlines=-1
c
      end module
c#######################################################################
      module step_size_stats
c
c-----------------------------------------------------------------------
c ****** Variable step size statistics.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      logical :: gather_stats
c
      integer(8) :: stat_n=0
      real(r_typ) :: stat_ds_sum=0._r_typ
      real(r_typ) :: stat_ds_avg=0._r_typ
      real(r_typ) :: stat_ds_min=huge(0._r_typ)
      real(r_typ) :: stat_ds_max=0._r_typ
c
      end module
c#######################################################################
      module openmp_vars
c
c-----------------------------------------------------------------------
c ****** Variables to control OpenMP parallelization.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Number of iterations to do in each thread.
c
      integer :: iterations_per_thread=500
c
      end module
c#######################################################################
      module params
c
c-----------------------------------------------------------------------
c ****** Parameters.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      character(512) :: infile
      logical :: verbose
c
      end module
c#######################################################################
      module interp_interface
      interface
        subroutine interp (n,x,xv,i,ip1,alpha,tab)
        use number_types
        use invint_def
        use locate_interval_interface
        integer :: n
        real(r_typ), dimension(n) :: x
        real(r_typ) :: xv
        integer :: i
        integer :: ip1
        real(r_typ) :: alpha
        type(itab), optional :: tab
        end
      end interface
      end module
c#######################################################################
      module tracefl_interface
      interface
        subroutine tracefl (b,ds,s0,s1,bs0,bs1,s,
     &                      traced_to_r_boundary,xt)
        use number_types
        use types
        type(vec) :: b
        type(flparam) :: ds
        real(r_typ), dimension(3) :: s0,s1
        real(r_typ), dimension(3) :: bs0,bs1
        real(r_typ) :: s
        logical :: traced_to_r_boundary
        type(traj), optional :: xt
        end
      end interface
      end module
c#######################################################################
      program MAPFL
c
c-----------------------------------------------------------------------
c
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
c ****** The number of lines required in the input file.
c
      integer, parameter :: input_file_lines_required=88
c
c-----------------------------------------------------------------------
c
      integer :: ierr
      integer :: nlines
      real(r_typ) :: ch_map_r
      logical :: trace_fwd,trace_bwd
      logical :: compute_q_on_slice
      logical :: compute_ch_map
      logical :: compute_ch_map_3d
c
c-----------------------------------------------------------------------
c
c ****** Set the parameters.
c
      call set_parameters
c
c ****** Check that the input file has the correct number of lines
c ****** to reduce the chance of using an invalid input file.
c
      call ffopen (1,infile,'r',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in MAPFL:'
        write (*,*) '### The input file does not exist'//
     &              ' or cannot be read.'
        write (*,*) 'File name: ',trim(infile)
        call exit (1)
      end if
c
      nlines=0
      do while (.true.)
        read (1,*,iostat=ierr)
        if (ierr.lt.0) exit
        nlines=nlines+1
      enddo
c
      close (1)
c
      if (nlines.ne.input_file_lines_required) then
        write (*,*)
        write (*,*) '### ERROR in MAPFL:'
        write (*,*) '### The input file does not have the'//
     &              ' correct number of lines:'
        write (*,*) 'Number of lines required = ',
     &              input_file_lines_required
        write (*,*) 'Number of lines present  = ',nlines
        call exit (1)
      end if
c
c ****** Read the input file.
c
      call ffopen (1,trim(infile),'r',ierr)
c
      if (ierr.ne.0) call exit (1)
c
      read (1,*)
      read (1,*) debug_level
      read (1,*)
      read (1,*) domain_r_min
      read (1,*)
      read (1,*) domain_r_max
      read (1,*)
      read (1,*) bfile%r
      read (1,*)
      read (1,*) bfile%t
      read (1,*)
      read (1,*) bfile%p
      read (1,*)
      read (1,*) cubic
      read (1,*)
      read (1,*) ds%variable
      read (1,*)
      read (1,*) ds%over_rc
      read (1,*)
      read (1,*) set_ds_automatically
      read (1,*)
      read (1,*) ds%min
      read (1,*)
      read (1,*) ds%max
      read (1,*)
      read (1,*) dsmult
      read (1,*)
      read (1,*) ds%limit_by_local_mesh
      read (1,*)
      read (1,*) ds%local_mesh_factor
      read (1,*)
      read (1,*) ds%lmax
      read (1,*)
      read (1,*) trace_fwd
      read (1,*)
      read (1,*) trace_bwd
      read (1,*)
      read (1,*) rffile
      read (1,*)
      read (1,*) tffile
      read (1,*)
      read (1,*) pffile
      read (1,*)
      read (1,*) qffile
      read (1,*)
      read (1,*) rbfile
      read (1,*)
      read (1,*) tbfile
      read (1,*)
      read (1,*) pbfile
      read (1,*)
      read (1,*) qbfile
      read (1,*)
      read (1,*) new_t_mesh
      read (1,*)
      read (1,*) mesh_file_t
      read (1,*)
      read (1,*) ntss
      read (1,*)
      read (1,*) t0,t1
      read (1,*)
      read (1,*) new_p_mesh
      read (1,*)
      read (1,*) mesh_file_p
      read (1,*)
      read (1,*) npss
      read (1,*)
      read (1,*) p0,p1
      read (1,*)
      read (1,*) compute_q_on_slice
      read (1,*)
      read (1,*) slice_coords_are_xyz
      read (1,*)
      read (1,*) q_increment_h
      read (1,*)
      read (1,*) slice_input_file%r
      read (1,*)
      read (1,*) slice_input_file%t
      read (1,*)
      read (1,*) slice_input_file%p
      read (1,*)
      read (1,*) slice_q_output_file
      read (1,*)
      read (1,*) compute_ch_map
      read (1,*)
      read (1,*) ch_map_r
      read (1,*)
      read (1,*) ch_map_output_file
c
      close (1)
c
      if (verbose) then
        write (*,*)
        write (*,*) '### ',cname,' Version ',cvers,' of ',cdate,'.'
      end if
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
      call readb (bfile,b)
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
        call read_slice_coordinates
        call get_q_on_slice
        call deallocate_slice_coordinates
      end if
c
c ****** Compute a coronal hole map, if requested.
c
      if (compute_ch_map) then
        call get_ch_map (ch_map_r)
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
c
      call exit (0)
c
      end
c#######################################################################
      subroutine readb (bfile,b)
c
c-----------------------------------------------------------------------
c
c ****** Read the magnetic field from the files specified by
c ****** BFILE into the magnetic field vector structure B.
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
      type(vfile) :: bfile
      type(vec) :: b
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Read the magnetic field components.
c
c ****** Br.
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Reading data file: ',trim(bfile%r)
      end if
c
      call rdhdf (bfile%r,b%r,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Br.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%r)
        call exit (1)
      end if
c
      if (b%r%ndim.ne.3.or..not.b%r%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Br file.'
        write (*,*) 'File name: ',trim(bfile%r)
        call exit (1)
      end if
c
c ****** Bt.
c
      if (verbose) then
        write (*,*) 'Reading data file: ',trim(bfile%t)
      end if
c
      call rdhdf (bfile%t,b%t,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Bt.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%t)
        call exit (1)
      end if
c
      if (b%t%ndim.ne.3.or..not.b%t%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Bt file.'
        write (*,*) 'File name: ',trim(bfile%t)
        call exit (1)
      end if
c
c ****** Bp.
c
      if (verbose) then
        write (*,*) 'Reading data file: ',trim(bfile%p)
      end if
c
      call rdhdf (bfile%p,b%p,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Bp.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%p)
        call exit (1)
      end if
c
      if (b%p%ndim.ne.3.or..not.b%p%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Bp file.'
        write (*,*) 'File name: ',trim(bfile%p)
        call exit (1)
      end if
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
c ****** Determine the primary (r,t,p) scales and the mesh limits
c ****** from the type of magnetic field in structure B, and
c ****** store them in structure B.
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
c ****** Tolerance for checking the extent of phi scales.
c
      real(r_typ), parameter :: eps=2.e-6_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: z
      integer :: n1,n2,n3
      logical :: add_phi_point
c
c-----------------------------------------------------------------------
c
c ****** Check the type of magnetic field files read in.
c
c ****** Check for new MAS code files.
c
      if (b%r%dims(1).eq.b%t%dims(1)+1.and.
     &    b%r%dims(1).eq.b%p%dims(1)+1.and.
     &    b%r%dims(2).eq.b%t%dims(2)-1.and.
     &    b%r%dims(2).eq.b%p%dims(2)  .and.
     &    b%r%dims(3).eq.b%t%dims(3)  .and.
     &    b%r%dims(3).eq.b%p%dims(3)-1) then
c
        b%nrs=b%r%dims(1)-1
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
c
        b%rs=>b%t%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
        add_phi_point=.false.
c
c ****** Check for old MAS code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)+1.and.
     &         b%r%dims(1).eq.b%p%dims(1)+1.and.
     &         b%r%dims(2).eq.b%t%dims(2)-1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)  ) then
c
        b%nrs=b%r%dims(1)-1
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
c
        b%rs=>b%t%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
        add_phi_point=.true.
c
c ****** Check for new POT3D code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)+1) then
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)-1
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%p%scales(3)%f
c
        add_phi_point=.false.
c
c ****** Check for new-old POT3D code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)-1) then
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
        add_phi_point=.false.
c
c ****** Check for old POT3D code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)  ) then
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
        add_phi_point=.true.
c
c ****** Check for non-staggered files.
c
      else if (b%r%dims(1).eq.b%t%dims(1).and.
     &         b%r%dims(1).eq.b%p%dims(1).and.
     &         b%r%dims(2).eq.b%t%dims(2).and.
     &         b%r%dims(2).eq.b%p%dims(2).and.
     &         b%r%dims(3).eq.b%t%dims(3).and.
     &         b%r%dims(3).eq.b%p%dims(3)) then
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
c ****** Do not add a phi point if the phi interval already includes
c ****** the whole interval.
c
        if (abs((b%ps(b%nps)-b%ps(1))-twopi).lt.eps) then
          add_phi_point=.false.
        else
          add_phi_point=.true.
        end if
c
      else
c
c ****** Invalid file type.
c
        write (*,*)
        write (*,*) '### ERROR in SET_BTYPE:'
        write (*,*) '### Unrecognized magnetic field file staggering.'
        call exit (1)
c
      end if
c
c ****** If appropriate, add a point in the phi dimension to
c ****** take care of periodic wrap-around.
c
      if (add_phi_point) then
c
        n1=b%r%dims(1)
        n2=b%r%dims(2)
        n3=b%r%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%r%f(:,:,:)
        f(:,:,n3+1)=b%r%f(:,:,1)
        z(1:n3)=b%r%scales(3)%f(:)
        z(n3+1)=b%r%scales(3)%f(1)+twopi
        deallocate (b%r%f)
        deallocate (b%r%scales(3)%f)
        b%r%dims(3)=n3+1
        b%r%f=>f
        b%r%scales(3)%f=>z
c
        n1=b%t%dims(1)
        n2=b%t%dims(2)
        n3=b%t%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%t%f(:,:,:)
        f(:,:,n3+1)=b%t%f(:,:,1)
        z(1:n3)=b%t%scales(3)%f(:)
        z(n3+1)=b%t%scales(3)%f(1)+twopi
        deallocate (b%t%f)
        deallocate (b%t%scales(3)%f)
        b%t%dims(3)=n3+1
        b%t%f=>f
        b%t%scales(3)%f=>z
c
        n1=b%p%dims(1)
        n2=b%p%dims(2)
        n3=b%p%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%p%f(:,:,:)
        f(:,:,n3+1)=b%p%f(:,:,1)
        z(1:n3)=b%p%scales(3)%f(:)
        z(n3+1)=b%p%scales(3)%f(1)+twopi
        deallocate (b%p%f)
        deallocate (b%p%scales(3)%f)
        b%p%dims(3)=n3+1
        b%p%f=>f
        b%p%scales(3)%f=>z
c
        b%nps=b%r%dims(3)
        b%ps=>b%r%scales(3)%f
c
      end if
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
      subroutine get_dx (n,x,dx)
c
c-----------------------------------------------------------------------
c
c ****** Get the cell size DX(N) from the 1D mesh in X(N).
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x,dx
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      if (n.le.1) then
        dx(1)=0.
      else if (n.eq.2) then
        dx(1)=x(2)-x(1)
        dx(2)=dx(1)
      else
        do i=2,n-1
          dx(i)=half*(x(i+1)-x(i-1))
        enddo
        dx(1)=dx(2)
        dx(n)=dx(n-1)
      end if
c
      return
      end
c#######################################################################
      subroutine make_new_meshes (b)
c
c-----------------------------------------------------------------------
c
c ****** Make new t and p meshes, if requested, or link them
c ****** to the meshes in the B files.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use params
      use files
      use vars
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
      type(sds) :: s
      integer :: i,ierr
      real(r_typ) :: d
c
c-----------------------------------------------------------------------
c
c ****** Make the t mesh.
c
      if (new_t_mesh) then
c
c ****** Check if the mesh is to be read from a 1D HDF file.
c
        if (mesh_file_t.ne.' ') then
c
          if (verbose) then
            write (*,*)
            write (*,*) '### Reading the t mesh from file: ',
     &                  trim(mesh_file_t)
          end if
c
          call rdhdf (mesh_file_t,s,ierr)
c
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the t mesh'//
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_t)
            call exit (1)
          end if
c
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the t mesh'//
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_t)
            call exit (1)
          end if
c
          ntss=s%dims(1)
          allocate (tss(ntss))
          tss=s%f(:,1,1)
c
          call deallocate_sds (s)
c
          if (verbose) then
            write (*,*)
            write (*,*) '### Mesh read in for the t mesh:'
            write (*,*) 'Number of points = ',ntss
            write (*,*) 'Lower limit = ',tss(1)
            write (*,*) 'Upper limit = ',tss(ntss)
          end if
c
        else
c
c ****** Generate a uniform mesh.
c
          if (ntss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//
     &                  ' for the uniform t mesh.'
            write (*,*) 'Number of points specified = ',ntss
            call exit (1)
          end if
c
          allocate (tss(ntss))
c
          if (t0.eq.0..and.t1.eq.0.) then
            t0=b%lim0(2)
            t1=b%lim1(2)
          end if
c
          if (verbose) then
            write (*,*)
            write (*,*) '### Generating a uniform t mesh:'
            write (*,*) 'Number of points = ',ntss
            write (*,*) 'Lower limit = ',t0
            write (*,*) 'Upper limit = ',t1
          end if
c
          if (ntss.ne.1) then
            d=(t1-t0)/(ntss-1)
          else
            d=0.
          end if
c
          do i=1,ntss
            tss(i)=t0+(i-1)*d
          enddo
          tss(1)=t0
          if (ntss.gt.1) then
            tss(ntss)=t1
          end if
c
        end if
c
      else
c
c ****** Use the same mesh as the primary B field mesh.
c
        ntss=b%nts
        tss=>b%ts
c
      end if
c
c ****** Make the p mesh.
c
      if (new_p_mesh) then
c
c ****** Check if the mesh is to be read from a 1D HDF file.
c
        if (mesh_file_p.ne.' ') then
c
          if (verbose) then
            write (*,*)
            write (*,*) '### Reading the p mesh from file: ',
     &                  trim(mesh_file_p)
          end if
c
          call rdhdf (mesh_file_p,s,ierr)
c
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the p mesh'//
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_p)
            call exit (1)
          end if
c
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the p mesh'//
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_p)
            call exit (1)
          end if
c
          npss=s%dims(1)
          allocate (pss(npss))
          pss=s%f(:,1,1)
c
          call deallocate_sds (s)
c
          if (verbose) then
            write (*,*)
            write (*,*) '### Mesh read in for the p mesh:'
            write (*,*) 'Number of points = ',npss
            write (*,*) 'Lower limit = ',pss(1)
            write (*,*) 'Upper limit = ',pss(npss)
          end if
c
        else
c
c ****** Generate a uniform mesh.
c
          if (npss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//
     &                  ' for the uniform p mesh.'
            write (*,*) 'Number of points specified = ',npss
            call exit (1)
          end if
c
          allocate (pss(npss))
c
          if (p0.eq.0..and.p1.eq.0.) then
            p0=b%lim0(3)
            p1=b%lim1(3)
          end if
c
          if (verbose) then
            write (*,*)
            write (*,*) '### Generating a uniform p mesh:'
            write (*,*) 'Number of points = ',npss
            write (*,*) 'Lower limit = ',p0
            write (*,*) 'Upper limit = ',p1
          end if
c
          if (npss.ne.1) then
            d=(p1-p0)/(npss-1)
          else
            d=0.
          end if
c
          do i=1,npss
            pss(i)=p0+(i-1)*d
          enddo
          pss(1)=p0
          if (npss.gt.1) then
            pss(npss)=p1
          end if
c
        end if
c
      else
c
c ****** Use the same mesh as the primary B field mesh.
c
        npss=b%nps
        pss=>b%ps
c
      end if
c
      return
      end
c#######################################################################
      subroutine set_ds (b)
c
c-----------------------------------------------------------------------
c
c ****** Set the field line integration step size.
c
c-----------------------------------------------------------------------
c
c ****** If SET_DS_AUTOMATICALLY=.T., the miniumum step size is set
c ****** to the minimum of the cell dimensions from the magnetic
c ****** field files, and the maximum step size is set to the
c ****** maximum of the cell dimensions.  Otherwise, the values read
c ****** in for DS%MIN and DS%MAX are used.
c
c ****** After being set in the above way, DS%MIN and DS%MAX are
c ****** multiplied by the factor DSMULT.  Thus, DSMULT provides a
c ****** quick way to change the integration step size.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use constants
      use vars
      use params
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
      integer :: i,j,k
      real(r_typ) :: dr,dt,dp
      real(r_typ) :: drmin,dtmin,dpmin
      real(r_typ) :: drmax,dtmax,dpmax
c
c-----------------------------------------------------------------------
c
      if (set_ds_automatically) then
c
        drmin=abs(b%lim1(1)-b%lim0(1))
        drmax=0.
        do i=1,b%nrs-1
          dr=abs(b%rs(i+1)-b%rs(i))
          drmin=min(drmin,dr)
          drmax=max(drmax,dr)
        enddo
c
        dtmin=pi
        dtmax=0.
        do j=1,b%nts-1
          dt=abs(b%ts(j+1)-b%ts(j))
          dtmin=min(dtmin,dt)
          dtmax=max(dtmax,dt)
        enddo
c
        dpmin=twopi
        dpmax=0.
        do k=1,b%nps-1
          dp=abs(b%ps(k+1)-b%ps(k))
          dpmin=min(dpmin,dp)
          dpmax=max(dpmax,dp)
        enddo
c
        ds%min=min(drmin,b%lim0(1)*dtmin,b%lim0(1)*dtmin*dpmin)
        ds%max=max(drmax,b%lim1(1)*dtmax,b%lim1(1)*dpmax)
c
      end if
c
      if (dsmult.le.0.) then
        write (*,*)
        write (*,*) '### ERROR in SET_DS:'
        write (*,*) '### DSMULT must be positive.'
        write (*,*) 'DSMULT= ',dsmult
        call exit (1)
      end if
c
      ds%over_rc=ds%over_rc*dsmult
      ds%min=ds%min*dsmult
      ds%max=ds%max*dsmult
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Field line integration parameters:'
        if (ds%variable) then
          write (*,*)
          write (*,*) '### Integration step size control:'//
     &                ' variable step size'
          write (*,*)
          write (*,*) '### Step size parameters:'
          write (*,*) 'DS%OVER_RC = ',ds%over_rc
          write (*,*) 'DS%MIN = ',ds%min
          write (*,*) 'DS%MAX = ',ds%max
          write (*,*)
          if (ds%limit_by_local_mesh) then
            write (*,*) '### Step size limited by the local'//
     &                  ' B mesh: yes'
            write (*,*) 'DS%LOCAL_MESH_FACTOR = ',
     &                  ds%local_mesh_factor
          else
            write (*,*) '### Step size limited by the local'//
     &                  ' B mesh: no'
          end if
        else
          write (*,*)
          write (*,*) '### Integration step size control:'//
     &                ' uniform step size'
          write (*,*)
          write (*,*) '### Step size parameters:'
          write (*,*) 'DS = ',ds%min
        end if
      end if
c
      return
      end
c#######################################################################
      subroutine map_forward
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines outward from r=R0.
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
      use field_line_params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
c
c-----------------------------------------------------------------------
c
c ****** Storage for the mapping.
c
      real(r_typ), dimension(ntss,npss) :: rfl,tfl,pfl,efl,kfl
c
      real(r_typ), dimension(:), allocatable :: tssh,pssh
      real(r_typ), dimension(:,:), allocatable :: qfl
c
c-----------------------------------------------------------------------
c
      integer :: ierr,j,k,nbad
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      real(r_typ) :: dtdt_m,dtdt_p
      real(r_typ) :: dtdp_m,dtdp_p
      real(r_typ) :: dpdt_m,dpdt_p
      real(r_typ) :: dpdp_m,dpdp_p
      real(r_typ) :: dtdt,dtdp,dpdt,dpdp
      real(r_typ) :: dt,dp,aa,bb,cc,dd,stm,stp,tmav,efav
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: modulo_twopi
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines, starting from each (T,P) cell at r=R0,
c ****** until the field line hits r=R1, or goes back to r=R0,
c ****** or exhausts the field line length allowed.
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a forward mapping from R0:'
      end if
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      nbad=0
c
      ds%direction_is_along_b=.false.
      ds%direction=1
c
      n_total=ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
          if (verbose) then
c$omp critical
            n_completed=n_completed+1
            nc=n_completed
c$omp end critical
          end if
c
          xfl0(1)=b%lim0(1)
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
c
          call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1, and set
c ****** the expansion factor.
c
          if (ttb) then
            rfl(j,k)=xfl1(1)
            tfl(j,k)=xfl1(2)
            pfl(j,k)=xfl1(3)
            if (bs1(1).ne.0.) then
              efl(j,k)=abs((bs0(1)*xfl0(1)**2)/(bs1(1)*xfl1(1)**2))
            else
              efl(j,k)=0.
            end if
            if (bs1(1).ne.0.) then
              kfl(j,k)=log10(max(abs(bs0(1)/bs1(1)),tiny(bs0(1))))
            else
              kfl(j,k)=-50._r_typ
            end if
          else
c$omp critical
            nbad=nbad+1
            write (*,*)
            write (*,*) '### WARNING from MAP_FORWARD:'
            write (*,*) '### A field line did not reach R0 or R1.'
            write (*,*) 'Initial theta = ',xfl0(2)
            write (*,*) 'Initial phi   = ',xfl0(3)
            write (*,*) 'Final field line radius = ',xfl1(1)
c$omp end critical
            rfl(j,k)=-1._r_typ
            tfl(j,k)=-1._r_typ
            pfl(j,k)=-1._r_typ
            efl(j,k)=0.
            kfl(j,k)=-50._r_typ
          end if
c
          if (max_bad_fieldlines.gt.0) then
            if (nbad.gt.max_bad_fieldlines) then
c$omp critical
              write (*,*)
              write (*,*) '### ERROR in MAP_FORWARD:'
              write (*,*) '### Too many field lines did not reach'//
     &                    ' R0 or R1.'
              write (*,*) 'Number of bad traces = ',max_bad_fieldlines
              call exit (1)
c$omp end critical
            end if
          end if
c
c ****** Write progress diagnostics if requested. 
c
          if (verbose) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
c
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the mapping.
c
      wrote_cr=.false.
c
      if (rffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//
     &                'r to file: ',
     &                trim(rffile)
        end if
        call wrhdf_2d (rffile,.true.,ntss,npss,rfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
c
      if (tffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//
     &                't to file: ',
     &                trim(tffile)
        end if
        call wrhdf_2d (tffile,.true.,ntss,npss,tfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
c
      if (pffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//
     &                'p to file: ',
     &                trim(pffile)
        end if
        call wrhdf_2d (pffile,.true.,ntss,npss,pfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
c
c ****** Compute Q (if requested).
c
      if (qffile.eq.' ') return
c
c ****** This can only be done if NTSS and NPSS exceed 1.
c
      if (ntss.le.1.or.npss.le.1) then
        write (*,*)
        write (*,*) '### WARNING from MAP_FORWARD:'
        write (*,*) '### Could not compute the Q factor.'
        write (*,*) '### To compute Q, NTSS and NPSS'//
     &              ' must be greater than 1.'
        return
      end if
c
      allocate (tssh(ntss-1))
      allocate (pssh(npss-1))
      allocate (qfl(ntss-1,npss-1))
c
c ****** Define the half-meshes (on which Q is computed).
c
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
c
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
c
      do k=1,npss-1
        dp=pss(k+1)-pss(k)
        do j=1,ntss-1
          dt=tss(j+1)-tss(j)
          if (efl(j  ,k  ).eq.0..or.efl(j+1,k  ).eq.0..or.
     &        efl(j  ,k+1).eq.0..or.efl(j+1,k+1).eq.0.) then
            qfl(j,k)=0.
          else
            efav=quarter*(efl(j,k  )+efl(j+1,k  )+
     &                    efl(j,k+1)+efl(j+1,k+1))
            if (efav.ne.0.) then
              tmav=quarter*(tfl(j,k  )+tfl(j+1,k  )+
     &                      tfl(j,k+1)+tfl(j+1,k+1))
              stm=sin(tmav)
              stp=sin(tssh(j))
              dtdt_m=(tfl(j+1,k  )-tfl(j  ,k  ))/dt
              dtdt_p=(tfl(j+1,k+1)-tfl(j  ,k+1))/dt
              dtdp_m=(tfl(j  ,k+1)-tfl(j  ,k  ))/dp
              dtdp_p=(tfl(j+1,k+1)-tfl(j+1,k  ))/dp
              dpdt_m=modulo_twopi(pfl(j+1,k  )-pfl(j  ,k  ))/dt
              dpdt_p=modulo_twopi(pfl(j+1,k+1)-pfl(j  ,k+1))/dt
              dpdp_m=modulo_twopi(pfl(j  ,k+1)-pfl(j  ,k  ))/dp
              dpdp_p=modulo_twopi(pfl(j+1,k+1)-pfl(j+1,k  ))/dp
              dtdt=half*(dtdt_m+dtdt_p)
              dtdp=half*(dtdp_m+dtdp_p)
              dpdt=half*(dpdt_m+dpdt_p)
              dpdp=half*(dpdp_m+dpdp_p)
              aa=stm*dpdp/stp
              bb=stm*dpdt
              cc=dtdp/stp
              dd=dtdt
              qfl(j,k)=(aa**2+bb**2+cc**2+dd**2)/efav
            else
              qfl(j,k)=0.
            end if
          end if
        enddo
      enddo
c
      if (verbose) then
        if (.not.wrote_cr) then
          write (*,*)
          wrote_cr=.true.
        end if
        write (*,*) 'Writing the forward mapping '//
     &              'Q to file: ',
     &              trim(qffile)
      end if
      call wrhdf_2d (qffile,.true.,ntss-1,npss-1,qfl,tssh,pssh,
     &               hdf32,ierr)
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in MAP_FORWARD:'
        write (*,*) '### Could not write the Q factor file.'
        call exit (1)
      end if
c
      deallocate (tssh)
      deallocate (pssh)
      deallocate (qfl)
c
      return
      end
c#######################################################################
      subroutine map_backward
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines inward from r=R1.
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
      use field_line_params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
c
c-----------------------------------------------------------------------
c
c ****** Storage for the mapping.
c
      real(r_typ), dimension(ntss,npss) :: rfl,tfl,pfl,efl,kfl
c
      real(r_typ), dimension(:), allocatable :: tssh,pssh
      real(r_typ), dimension(:,:), allocatable :: qfl
c
c-----------------------------------------------------------------------
c
      integer :: ierr,j,k,nbad
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      real(r_typ) :: dtdt_m,dtdt_p
      real(r_typ) :: dtdp_m,dtdp_p
      real(r_typ) :: dpdt_m,dpdt_p
      real(r_typ) :: dpdp_m,dpdp_p
      real(r_typ) :: dtdt,dtdp,dpdt,dpdp
      real(r_typ) :: dt,dp,aa,bb,cc,dd,stm,stp,tmav,efav
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: modulo_twopi
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines, starting from each (T,P) cell at r=R1,
c ****** until the field line hits r=R0, or goes back to r=R1,
c ****** or exhausts the number of segments allowed.
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a backward mapping from R1:'
      end if
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      nbad=0
c
      ds%direction_is_along_b=.false.
      ds%direction=-1
c
      n_total=ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
          if (verbose) then
c$omp critical
            n_completed=n_completed+1
            nc=n_completed
c$omp end critical
          end if
c
          xfl0(1)=b%lim1(1)
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
c
          call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
c
          if (ttb) then
            rfl(j,k)=xfl1(1)
            tfl(j,k)=xfl1(2)
            pfl(j,k)=xfl1(3)
            if (bs0(1).ne.0.) then
              efl(j,k)=abs((bs1(1)*xfl1(1)**2)/(bs0(1)*xfl0(1)**2))
            else
              efl(j,k)=0.
            end if
            if (bs1(1).ne.0.) then
              kfl(j,k)=log10(max(abs(bs0(1)/bs1(1)),tiny(bs0(1))))
            else
              kfl(j,k)=-50._r_typ
            end if
          else
c$omp critical
            nbad=nbad+1
            write (*,*)
            write (*,*) '### WARNING from MAP_BACKWARD:'
            write (*,*) '### A field line did not reach R0 or R1.'
            write (*,*) 'Initial theta = ',xfl0(2)
            write (*,*) 'Initial phi   = ',xfl0(3)
            write (*,*) 'Final field line radius = ',xfl1(1)
c$omp end critical
            rfl(j,k)=-1._r_typ
            tfl(j,k)=-1._r_typ
            pfl(j,k)=-1._r_typ
            efl(j,k)=0.
            kfl(j,k)=-50._r_typ
          end if
c
          if (max_bad_fieldlines.gt.0) then
            if (nbad.gt.max_bad_fieldlines) then
c$omp critical
              write (*,*)
              write (*,*) '### ERROR in MAP_BACKWARD:'
              write (*,*) '### Too many field lines did not reach'//
     &                    ' R0 or R1.'
              write (*,*) 'Number of bad traces = ',max_bad_fieldlines
              call exit (1)
c$omp end critical
            end if
          end if
c
c ****** Write progress diagnostics if requested. 
c
          if (verbose) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
c
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the mapping.
c
      wrote_cr=.false.
c
      if (rbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//
     &                'r to file: ',
     &                trim(rbfile)
        end if
        call wrhdf_2d (rbfile,.true.,ntss,npss,rfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
c
      if (tbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//
     &                't to file: ',
     &                trim(tbfile)
        end if
        call wrhdf_2d (tbfile,.true.,ntss,npss,tfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
c
      if (pbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//
     &                'p to file: ',
     &                trim(pbfile)
        end if
        call wrhdf_2d (pbfile,.true.,ntss,npss,pfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
c
c ****** Compute Q (if requested).
c
      if (qbfile.eq.' ') return
c
c ****** This can only be done if NTSS and NPSS exceed 1.
c
      if (ntss.le.1.or.npss.le.1) then
        write (*,*)
        write (*,*) '### WARNING from MAP_BACKWARD:'
        write (*,*) '### Could not compute the Q factor.'
        write (*,*) '### To compute Q, NTSS and NPSS'//
     &              ' must be greater than 1.'
        return
      end if
c
      allocate (tssh(ntss-1))
      allocate (pssh(npss-1))
      allocate (qfl(ntss-1,npss-1))
c
c ****** Define the half-meshes (on which Q is computed).
c
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
c
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
c
      do k=1,npss-1
        dp=pss(k+1)-pss(k)
        do j=1,ntss-1
          dt=tss(j+1)-tss(j)
          if (efl(j  ,k  ).eq.0..or.efl(j+1,k  ).eq.0..or.
     &        efl(j  ,k+1).eq.0..or.efl(j+1,k+1).eq.0.) then
            qfl(j,k)=0.
          else
            efav=quarter*(efl(j,k  )+efl(j+1,k  )+
     &                    efl(j,k+1)+efl(j+1,k+1))
            if (efav.ne.0.) then
              tmav=quarter*(tfl(j,k  )+tfl(j+1,k  )+
     &                      tfl(j,k+1)+tfl(j+1,k+1))
              stm=sin(tmav)
              stp=sin(tssh(j))
              dtdt_m=(tfl(j+1,k  )-tfl(j  ,k  ))/dt
              dtdt_p=(tfl(j+1,k+1)-tfl(j  ,k+1))/dt
              dtdp_m=(tfl(j  ,k+1)-tfl(j  ,k  ))/dp
              dtdp_p=(tfl(j+1,k+1)-tfl(j+1,k  ))/dp
              dpdt_m=modulo_twopi(pfl(j+1,k  )-pfl(j  ,k  ))/dt
              dpdt_p=modulo_twopi(pfl(j+1,k+1)-pfl(j  ,k+1))/dt
              dpdp_m=modulo_twopi(pfl(j  ,k+1)-pfl(j  ,k  ))/dp
              dpdp_p=modulo_twopi(pfl(j+1,k+1)-pfl(j+1,k  ))/dp
              dtdt=half*(dtdt_m+dtdt_p)
              dtdp=half*(dtdp_m+dtdp_p)
              dpdt=half*(dpdt_m+dpdt_p)
              dpdp=half*(dpdp_m+dpdp_p)
              aa=stm*dpdp/stp
              bb=stm*dpdt
              cc=dtdp/stp
              dd=dtdt
              qfl(j,k)=(aa**2+bb**2+cc**2+dd**2)*efav
            else
              qfl(j,k)=0.
            end if
          end if
        enddo
      enddo
c
      if (verbose) then
        if (.not.wrote_cr) then
          write (*,*)
          wrote_cr=.true.
        end if
        write (*,*) 'Writing the backward mapping '//
     &              'Q to file: ',
     &              trim(qbfile)
      end if
      call wrhdf_2d (qbfile,.true.,ntss-1,npss-1,qfl,tssh,pssh,
     &               hdf32,ierr)
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in MAP_BACKWARD:'
        write (*,*) '### Could not write the Q factor file.'
        call exit (1)
      end if
c
      deallocate (tssh)
      deallocate (pssh)
      deallocate (qfl)
c
      return
      end
c#######################################################################
      function modulo_twopi (x)
c
c-----------------------------------------------------------------------
c
c ****** Return the smallest value of X, modulo 2*pi.
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x
      real(r_typ) :: modulo_twopi
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xm,xp,x_min
c
c-----------------------------------------------------------------------
c
      xm=abs(x-twopi)
      xp=abs(x+twopi)
      x_min=min(xm,abs(x),xp)
      if (xm.eq.x_min) then
        modulo_twopi=x-twopi
      else if (xp.eq.x_min) then
        modulo_twopi=x+twopi
      else
        modulo_twopi=x
      end if
c
      return
      end
c#######################################################################
      subroutine read_slice_coordinates
c
c-----------------------------------------------------------------------
c
c ****** Read the coordinates that define the slice in the 3D volume.
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
c ****** Read the coordinates of the slice.
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Reading the coordinates of the slice ...'
      end if
c
c ****** Read the x/r coordinate file.
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Reading the '//slice_coord_name(1)//
     &              ' coordinate from file: ',
     &              trim(slice_input_file%r)
      end if
c
      call rdhdf (slice_input_file%r,slice_c1,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(1)//
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%r)
        call exit (1)
      end if
c
c ****** Read the y/t coordinate file.
c
      if (verbose) then
        write (*,*) 'Reading the '//slice_coord_name(2)//
     &              ' coordinate from file: ',
     &              trim(slice_input_file%t)
      end if
c
      call rdhdf (slice_input_file%t,slice_c2,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(2)//
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%t)
        call exit (1)
      end if
c
c ****** Check that the y/t coordinate has the same structure as the
c ****** x/r coordinate.
c
      if (.not.same_structure_sds(slice_c1,slice_c2)) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### The data sets for coordinates '//
     &              slice_coord_name(1)//' and '//
     &              slice_coord_name(2)//' do not have'//
     &              ' the same structure.'
        write (*,*) 'Coordinate '//slice_coord_name(1)//
     &              ' file name: ',trim(slice_input_file%r)
        write (*,*) 'Coordinate '//slice_coord_name(2)//
     &              ' file name: ',trim(slice_input_file%t)
        call exit (1)
      end if
c
c ****** Read the z/p coordinate file.
c
      if (verbose) then
        write (*,*) 'Reading the '//slice_coord_name(3)//
     &              ' coordinate from file: ',
     &              trim(slice_input_file%p)
      end if
c
      call rdhdf (slice_input_file%p,slice_c3,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(3)//
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%p)
        call exit (1)
      end if
c
c ****** Check that the z/p coordinate has the same structure as the
c ****** x/r coordinate.
c
      if (.not.same_structure_sds(slice_c1,slice_c3)) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### The data sets for coordinates '//
     &              slice_coord_name(1)//' and '//
     &              slice_coord_name(3)//' do not have'//
     &              ' the same structure.'
        write (*,*) 'Coordinate '//slice_coord_name(1)//
     &              ' file name: ',trim(slice_input_file%r)
        write (*,*) 'Coordinate '//slice_coord_name(3)//
     &              ' file name: ',trim(slice_input_file%p)
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine deallocate_slice_coordinates
c
c-----------------------------------------------------------------------
c
      use sds_def
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Read the coordinates of the slice.
c
      call deallocate_sds (slice_c1)
      call deallocate_sds (slice_c2)
      call deallocate_sds (slice_c3)
c
      return
      end
c#######################################################################
      function same_structure_sds (s1,s2)
c
c-----------------------------------------------------------------------
c
c ****** Check if the two data sets S1 and S2 have the same
c ****** structure.  If they do, return .TRUE; otherwise, return
c ****** .FALSE.
c
c-----------------------------------------------------------------------
c
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s1,s2
      logical :: same_structure_sds
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      same_structure_sds=.false.
c
      if (s1%ndim.ne.s2%ndim) return
c
      if (s1%scale.neqv.s2%scale) return
c
      do i=1,s1%ndim
        if (s1%dims(i).ne.s2%dims(i)) return
      enddo
c
      same_structure_sds=.true.
c
      return
      end
c#######################################################################
      subroutine get_ch_map (rv)
c
c-----------------------------------------------------------------------
c
c ****** Compute a coronal hole map at radius r=RV.
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
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
c
c-----------------------------------------------------------------------
c
c ****** Storage for the coronal hole map.
c
      real(r_typ), dimension(npss,ntss) :: ch
c
c-----------------------------------------------------------------------
c
      integer :: ierr,j,k
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: f_trace_reached_boundary
      logical :: b_trace_reached_boundary
      logical :: f_trace_on_r0,f_trace_on_r1
      logical :: b_trace_on_r0,b_trace_on_r1
      logical :: f_br_positive
      logical :: b_br_positive
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a coronal hole map at r = ',rv
      end if
c
c ****** Check that the radius specified is valid.
c
      if (rv.lt.b%lim0(1).or.rv.gt.b%lim1(1)) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### Invalid radius specified.'
        write (*,*) '### The radius is outside the domain limits:'
        write (*,*) 'Lower radial domain limit = ',b%lim0(1)
        write (*,*) 'Upper radial domain limit = ',b%lim1(1)
        write (*,*) 'Specified radius          = ',rv
        call exit (1)
      end if
c
c ****** Check that the coronal hole map output file name is not
c ****** blank, since this does not make sense.
c
      if (ch_map_output_file.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### A coronal hole map was requested, yet'//
     &              ' the output file name is blank.'
        call exit (1)
      end if
c
c ****** Set the tracing direction to be along the direction
c ****** of the magnetic field.
c
      ds%direction_is_along_b=.true.
c
      ds_f=ds
      ds_f%direction=1
c
      ds_b=ds
      ds_b%direction=-1
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      n_total=ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(f_trace_reached_boundary,f_br_positive)
c$omp& private(f_trace_on_r0,f_trace_on_r1)
c$omp& private(b_trace_reached_boundary,b_br_positive)
c$omp& private(b_trace_on_r0,b_trace_on_r1)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
          if (verbose) then
c$omp critical
            n_completed=n_completed+1
            nc=n_completed
c$omp end critical
          end if
c
          xfl0(1)=rv
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
c
c ****** Trace a field line in the forward direction along B.
c
          call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1.
c
          if (ttb) then
            f_trace_reached_boundary=.true.
            f_trace_on_r0=xfl1(1).eq.b%lim0(1)
            f_trace_on_r1=xfl1(1).eq.b%lim1(1)
            f_br_positive=bs1(1).ge.0.
          else
            f_trace_reached_boundary=.false.
          end if
c
c ****** Trace a field line in the backward direction along B.
c
          call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1.
c
          if (ttb) then
            b_trace_reached_boundary=.true.
            b_trace_on_r0=xfl1(1).eq.b%lim0(1)
            b_trace_on_r1=xfl1(1).eq.b%lim1(1)
            b_br_positive=bs1(1).ge.0.
          else
            b_trace_reached_boundary=.false.
          end if
c
c ****** Set the coronal hole map value.
c
c ****** Note that the following values are set in the output
c ****** coronal hole map:
c ******
c ******    -1: open field line with negative polarity
c ******     1: open field line with positive polarity
c ******     0: closed field line with both footpoints
c ******        on r=R0
c ******     2: closed field line with both footpoints
c ******        on r=R1
c ******    -2: field line that does not reach either
c ******        the r=R0 boundary or the r=R1 boundary
c
          if (f_trace_reached_boundary.and.
     &        b_trace_reached_boundary) then
            if (f_trace_on_r0.and.b_trace_on_r1) then
              if (f_br_positive) then
                ch(k,j)=one
              else
                ch(k,j)=-one
              end if
            else if (f_trace_on_r1.and.b_trace_on_r0) then
              if (b_br_positive) then
                ch(k,j)=one
              else
                ch(k,j)=-one
              end if
            else if (f_trace_on_r0.and.b_trace_on_r0) then
              ch(k,j)=0.
            else if (f_trace_on_r1.and.b_trace_on_r1) then
              ch(k,j)=two
            else
              ch(k,j)=-two
            end if
          else
            ch(k,j)=-two
          end if
c
c ****** Write progress diagnostics if requested. 
c
          if (verbose) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
c
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the coronal hole map.
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Writing the coronal hole map to file: ',
     &              trim(ch_map_output_file)
      end if
c
      call wrhdf_2d (ch_map_output_file,.true.,
     &               npss,ntss,ch,pss,tss,
     &               hdf32,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### Could not write the coronal hole map.'
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine get_q_on_slice
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines from points on a slice in the 3D volume,
c ****** getting Q at each point.
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
      use diags
      use openmp_vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Storage for the mapping.
c
      real(r_typ), dimension(:,:,:), allocatable, target :: qfl
c
c-----------------------------------------------------------------------
c
      type (sds) :: out
      integer :: ierr,i,j,k,n1,n2,n3
      real(r_typ), dimension(3) :: s0
      real(r_typ), dimension(3) :: c
      real(r_typ) :: h
      real(r_typ) :: q
      logical :: gotq
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing Q on a slice:'
      end if
c
      n1=slice_c1%dims(1)
      n2=slice_c1%dims(2)
      n3=slice_c1%dims(3)
c
c ****** Allocate the storage for the output Q array.
c
      allocate (qfl(n1,n2,n3))
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Getting Q ...'
        write (*,*)
      end if
c
c ****** Calculate Q at each point on the slice.
c
      ds%direction_is_along_b=.true.
      ds%direction=1
c
      n_total=n1*n2*n3
      n_completed=0
c
c$omp parallel do
c$omp& private(i,j,k,c,s0,q,gotq)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(3)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,n3
        do j=1,n2
          do i=1,n1
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
            if (verbose) then
c$omp critical
              n_completed=n_completed+1
              nc=n_completed
c$omp end critical
            end if
c
            if (slice_coords_are_xyz) then
              c=(/slice_c1%f(i,j,k),
     &            slice_c2%f(i,j,k),
     &            slice_c3%f(i,j,k)/)
              call c2s (c,s0)
            else
              s0=(/slice_c1%f(i,j,k),
     &             slice_c2%f(i,j,k),
     &             slice_c3%f(i,j,k)/)
            end if
c
            call getq (ds,s0,q_increment_h,q,gotq)
c
            if (gotq) then
              qfl(i,j,k)=q
            else
              qfl(i,j,k)=0.
            end if
c
c ****** Write progress diagnostics if requested. 
c
            if (verbose) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
c
          enddo
        enddo
      enddo
c$omp end parallel do
c
      if (verbose) then
        write (*,*)
      end if
c
c ****** Write the Q slice.
c
      if (slice_q_output_file.ne.' ') then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>qfl
        if (verbose) then
          write (*,*) 'Writing Q in the slice to file: ',
     &                trim(slice_q_output_file)
        end if
        call wrhdf (slice_q_output_file,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in GET_Q_ON_SLICE:'
          write (*,*) '### Could not write Q in the slice'//
     &                ' to file: ',trim(slice_q_output_file)
          call exit (1)
        end if
      end if
c
      deallocate (qfl)
c
      return
      end
c#######################################################################
      subroutine getq (ds,s0,h,q,valid)
c
c-----------------------------------------------------------------------
c
c ****** Obtain Q at the point given by the spherical coordinates in
c ****** vector S0 by tracing 5 field lines forwards and backwards
c ****** from S0 to the boundaries.
c
c ****** If Q was successfully obtained, return the Q value in
c ****** variable Q, and VALID=.T.; otherwise, return VALID=.F..
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use field
      use debug
      use diags
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(flparam) :: ds
      real(r_typ), dimension(3) :: s0
      real(r_typ) :: h
      real(r_typ) :: q
      logical :: valid
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
c ****** Set the limit for how to evaluate Q derivatives near the
c ****** poles.
c
      real(r_typ), parameter :: st_pole_limit_max=5.e-3_r_typ
c
c ****** When sin(theta) of the central field line is smaller than
c ****** ST_POLE_LIMIT_MAX, Cartesian basis vectors are used to
c ****** compute derivatives of the field line mapping; otherwise,
c ****** spherical basis vectors are used.
c
c-----------------------------------------------------------------------
c
c ****** Flag to write warning messages.
c
      logical, parameter :: write_warning_messages=.true.
c
c-----------------------------------------------------------------------
c
      logical :: ttb
      type(flparam) :: ds_local
      real(r_typ), dimension(3) :: x0,x00
      real(r_typ), dimension(3) :: s00
      real(r_typ), dimension(3) :: bs0,bc0
      real(r_typ), dimension(3) :: bs1_f,bs1_b
      real(r_typ) :: s
      real(r_typ) :: bmag,br0_f,br0_b
      real(r_typ), dimension(3) :: bhat,bhat_abs
      real(r_typ) :: r_f,r_b,t_f,t_b
      real(r_typ) :: st_f,st_b
      integer, dimension(1) :: index_min
      real(r_typ), dimension(3) :: e,e1,e2
      real(r_typ), dimension(3) :: s1_f,s1_b
      real(r_typ), dimension(3) :: s1_1p_f,s1_1m_f
      real(r_typ), dimension(3) :: s1_1p_b,s1_1m_b
      real(r_typ), dimension(3) :: s1_2p_f,s1_2m_f
      real(r_typ), dimension(3) :: s1_2p_b,s1_2m_b
      real(r_typ) :: a_f,b_f,c_f,d_f
      real(r_typ) :: a_b,b_b,c_b,d_b
      real(r_typ) :: nsq
      type(csvec) :: xcs
      type(inout) :: outside
      logical :: outside_1p,outside_1m
      logical :: outside_2p,outside_2m
      real(r_typ) :: dx_1p,dx_1m,dx_2p,dx_2m
      real(r_typ) :: x_1p_f,x_1m_f,y_1p_f,y_1m_f
      real(r_typ) :: x_2p_f,x_2m_f,y_2p_f,y_2m_f
      real(r_typ) :: x_1p_b,x_1m_b,y_1p_b,y_1m_b
      real(r_typ) :: x_2p_b,x_2m_b,y_2p_b,y_2m_b
      character(6) :: ch_seq
      logical :: save_trace_points
c
c ****** Field line trace storage buffer.
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: modulo_twopi
      logical, external :: outside_domain
c
c-----------------------------------------------------------------------
c
c ****** Set the flag to save field line traces.
c
      if (debug_level.ge.2) then
        save_trace_points=.true.
        call allocate_trajectory_buffer (xt)
c$omp critical
        diag_seq=diag_seq+1
        write (ch_seq,'(i6.6)') diag_seq
c$omp end critical
      else
        save_trace_points=.false.
      end if
c
c ****** Initialize Q and the validity flag.

      q=0.
      valid=.false.
c
c ****** Initilaize the local DS from that supplied in the
c ****** argument list.
c
      ds_local=ds
c
c ****** Set the flag to interpret the tracing direction as the
c ****** direction along the magnetic field line.
c
      ds_local%direction_is_along_b=.true.
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### COMMENT from GETQ:'
        write (*,*) '### Diagnostics for computation of Q:'
        write (*,*) 's0=',s0
      end if
c
c ****** Check that the requested position is not outside the domain.
c ****** If it is, return without calculating a valid Q.
c
      xcs%s=s0
      call sph_to_cart (xcs)
c
      if (outside_domain(b,xcs,outside)) return
c
c-----------------------------------------------------------------------
c ****** Trace the central field line.
c-----------------------------------------------------------------------
c
      s00=s0
c
      ds_local%direction=1
      if (save_trace_points) then
        call tracefl (b,ds_local,s00,s1_f,bs0,bs1_f,s,ttb,xt)
      else
        call tracefl (b,ds_local,s00,s1_f,bs0,bs1_f,s,ttb)
      end if
c
      if (debug_level.ge.2) then
c$omp critical
        call write_trace ('fl_00_f_'//ch_seq//'.dat',xt)
c$omp end critical
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line, forward trace:'
        write (*,*) 's0=',s00
        write (*,*) 's1=',s1_f
        write (*,*) 'b0=',bs0
        write (*,*) 'b1=',bs1_f
        write (*,*) 'ttb=',ttb
        write (*,*) 's=',s
      end if
c
      if (.not.ttb) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The central field line did not reach'//
     &                ' the domain boundaries'
          write (*,*) '### during the forward trace.'
          write (*,*) 'Central launch point: ',s0
        end if
        return
      end if
c
      ds_local%direction=-1
      if (save_trace_points) then
        call tracefl (b,ds_local,s00,s1_b,bs0,bs1_b,s,ttb,xt)
      else
        call tracefl (b,ds_local,s00,s1_b,bs0,bs1_b,s,ttb)
      end if
c
      if (debug_level.ge.2) then
c$omp critical
        call write_trace ('fl_00_b_'//ch_seq//'.dat',xt)
c$omp end critical
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line, backward trace:'
        write (*,*) 's0=',s00
        write (*,*) 's1=',s1_b
        write (*,*) 'b0=',bs0
        write (*,*) 'b1=',bs1_b
        write (*,*) 'ttb=',ttb
        write (*,*) 's=',s
      end if
c
      if (.not.ttb) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The central field line did not reach'//
     &                ' the domain boundaries'
          write (*,*) '### during the backward trace.'
          write (*,*) 'Central launch point: ',s0
        end if
        return
      end if
c
c ****** Save the central field line endpoints.
c
      r_f=s1_f(1)
      r_b=s1_b(1)
      t_f=s1_f(2)
      t_b=s1_b(2)
      st_f=sin(t_f)
      st_b=sin(t_b)
c
c ****** Get the radial component of the magnetic field at the
c ****** central field line endpoints.
c
      br0_f=bs1_f(1)
      br0_b=bs1_b(1)
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line:'
        write (*,*) 'r_f=',r_f
        write (*,*) 'r_b=',r_b
        write (*,*) 't_f=',t_f
        write (*,*) 't_b=',t_b
        write (*,*) 'br0_f=',br0_f
        write (*,*) 'br0_b=',br0_b
      end if
c
c ****** Generate the basis vectors that define the plane of
c ****** the Q computation at the starting location
c ****** (in Cartesian coordinates).
c
c ****** Initialize the position vector at the starting point
c ****** in Cartesian coordinates in X0.
c
      call s2c (s0,x0)
c
c ****** Get the Cartesian magnetic field vector at S0.
c
      call sv_to_cv (s0,bs0,bc0)
c
      bmag=sqrt(bc0(1)**2+bc0(2)**2+bc0(3)**2)
c
      if (debug_level.ge.2) then
        write (*,*) 'bs0=',bs0
        write (*,*) 'bc0=',bc0
        write (*,*) 'bmag=',bmag
      end if
c
c ****** If we hit a null point (B=0), exit with an error.
c
      if (bmag.eq.0.) then
        if (debug_level.ge.2) then
          write (*,*)
          write (*,*) '### WARNING in GETQ:'
          write (*,*) 'B = 0 at the launch point.'
          write (*,*) 'Exiting with an error ...'
        end if
        return
      end if
c
      bhat=bc0/bmag
c
      if (debug_level.ge.2) then
        write (*,*) 'bhat=',bhat
      end if
c
c ****** Select the unit vector that is most perpendicular to B.
c
      bhat_abs=abs(bhat)
      index_min=minloc(bhat_abs)
c
      if (debug_level.ge.2) then
        write (*,*) 'index_min=',index_min
      end if
c
      e=0.
      e(index_min(1))=one
c
      if (debug_level.ge.2) then
        write (*,*) 'e=',e
      end if
c
c ****** The triplet E1, E2, and BHAT form an orthogonal basis.
c
      call normalized_cross_product (e,bhat,e1)
      call normalized_cross_product (e1,bhat,e2)
c
      if (debug_level.ge.2) then
        write (*,*) 'e1=',e1
        write (*,*) 'e2=',e2
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0+H*E1.
c-----------------------------------------------------------------------
c
      x00=x0+h*e1
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_1p=.true.
        dx_1p=0.
        s1_1p_f=s1_f
        s1_1p_b=s1_b
      else
        outside_1p=.false.
        dx_1p=h
      end if
c
      if (.not.outside_1p) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1p_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1p_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_p0_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e1 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1p_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1p_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1p_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_p0_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e1 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1p_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0-H*E1.
c-----------------------------------------------------------------------
c
      x00=x0-h*e1
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_1m=.true.
        dx_1m=0.
        s1_1m_f=s1_f
        s1_1m_b=s1_b
      else
        outside_1m=.false.
        dx_1m=h
      end if
c
c ****** If both the plus and minus perturbed launch points
c ****** are outside the domain, return without calculating a
c ****** valid Q.  Write a warning if this happens.
c
      if (outside_1p.and.outside_1m) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The x0+h*e1 and x0-h*e1 launch points'//
     &                ' are both outside the domain.'
          write (*,*) 'Central launch point: ',s0
          x00=x0+h*e1
          call c2s (x00,s00)
          write (*,*) 'Launch point x0+h*e1: ',s00
          x00=x0-h*e1
          call c2s (x00,s00)
          write (*,*) 'Launch point x0-h*e1: ',s00
        end if
        return
      end if
c
      if (.not.outside_1m) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1m_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1m_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_m0_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e1 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1m_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1m_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1m_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_m0_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e1 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1m_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0+H*E2.
c-----------------------------------------------------------------------
c
      x00=x0+h*e2
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_2p=.true.
        dx_2p=0.
        s1_2p_f=s1_f
        s1_2p_b=s1_b
      else
        outside_2p=.false.
        dx_2p=h
      end if
c
      if (.not.outside_2p) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2p_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2p_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0p_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e2 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2p_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2p_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2p_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0p_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e2 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2p_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0-H*E2.
c-----------------------------------------------------------------------
c
      x00=x0-h*e2
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_2m=.true.
        dx_2m=0.
        s1_2m_f=s1_f
        s1_2m_b=s1_b
      else
        outside_2m=.false.
        dx_2m=h
      end if
c
c ****** If both the plus and minus perturbed launch points
c ****** are outside the domain, return without calculating a
c ****** valid Q.  This should never happen.
c
      if (outside_2p.and.outside_2m) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The x0+h*e2 and x0-h*e2 launch points'//
     &                ' are both outside the domain.'
          write (*,*) 'Central launch point: ',s0
          x00=x0+h*e2
          call c2s (x00,s00)
          write (*,*) 'Launch point x0+h*e2: ',s00
          x00=x0-h*e2
          call c2s (x00,s00)
          write (*,*) 'Launch point x0-h*e2: ',s00
        end if
        return
      end if
c
      if (.not.outside_2m) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2m_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2m_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0m_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e2 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2m_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: Spherical s0 = ',s0
            write (*,*) 'Launch point: Spherical x0-h*e2 = ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2m_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2m_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0m_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e2 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2m_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: Spherical s0 = ',s0
            write (*,*) 'Launch point: Spherical x0-h*e2 = ',s00
          end if
          return
        end if
c
      end if
c
c ****** Get the value of the Jacobian matrix coefficients.
c
      if (st_f.lt.st_pole_limit_max) then
        x_1p_f=r_f*sin(s1_1p_f(2))*cos(s1_1p_f(3))
        x_1m_f=r_f*sin(s1_1m_f(2))*cos(s1_1m_f(3))
        y_1p_f=r_f*sin(s1_1p_f(2))*sin(s1_1p_f(3))
        y_1m_f=r_f*sin(s1_1m_f(2))*sin(s1_1m_f(3))
        x_2p_f=r_f*sin(s1_2p_f(2))*cos(s1_2p_f(3))
        x_2m_f=r_f*sin(s1_2m_f(2))*cos(s1_2m_f(3))
        y_2p_f=r_f*sin(s1_2p_f(2))*sin(s1_2p_f(3))
        y_2m_f=r_f*sin(s1_2m_f(2))*sin(s1_2m_f(3))
        a_f=(x_1p_f-x_1m_f)/(dx_1p+dx_1m)
        b_f=(x_2p_f-x_2m_f)/(dx_2p+dx_2m)
        c_f=(y_1p_f-y_1m_f)/(dx_1p+dx_1m)
        d_f=(y_2p_f-y_2m_f)/(dx_2p+dx_2m)
      else
        a_f=modulo_twopi(s1_1p_f(3)-s1_1m_f(3))*r_f*st_f/(dx_1p+dx_1m)
        b_f=modulo_twopi(s1_2p_f(3)-s1_2m_f(3))*r_f*st_f/(dx_2p+dx_2m)
        c_f=(s1_1p_f(2)-s1_1m_f(2))*r_f/(dx_1p+dx_1m)
        d_f=(s1_2p_f(2)-s1_2m_f(2))*r_f/(dx_2p+dx_2m)
      end if
c
      if (st_b.lt.st_pole_limit_max) then
        x_1p_b=r_b*sin(s1_1p_b(2))*cos(s1_1p_b(3))
        x_1m_b=r_b*sin(s1_1m_b(2))*cos(s1_1m_b(3))
        y_1p_b=r_b*sin(s1_1p_b(2))*sin(s1_1p_b(3))
        y_1m_b=r_b*sin(s1_1m_b(2))*sin(s1_1m_b(3))
        x_2p_b=r_b*sin(s1_2p_b(2))*cos(s1_2p_b(3))
        x_2m_b=r_b*sin(s1_2m_b(2))*cos(s1_2m_b(3))
        y_2p_b=r_b*sin(s1_2p_b(2))*sin(s1_2p_b(3))
        y_2m_b=r_b*sin(s1_2m_b(2))*sin(s1_2m_b(3))
        a_b=(x_1p_b-x_1m_b)/(dx_1p+dx_1m)
        b_b=(x_2p_b-x_2m_b)/(dx_2p+dx_2m)
        c_b=(y_1p_b-y_1m_b)/(dx_1p+dx_1m)
        d_b=(y_2p_b-y_2m_b)/(dx_2p+dx_2m)
      else
        a_b=modulo_twopi(s1_1p_b(3)-s1_1m_b(3))*r_b*st_b/(dx_1p+dx_1m)
        b_b=modulo_twopi(s1_2p_b(3)-s1_2m_b(3))*r_b*st_b/(dx_2p+dx_2m)
        c_b=(s1_1p_b(2)-s1_1m_b(2))*r_b/(dx_1p+dx_1m)
        d_b=(s1_2p_b(2)-s1_2m_b(2))*r_b/(dx_2p+dx_2m)
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Jacobian matrix coefficients:'
        write (*,*) 'a_f=',a_f
        write (*,*) 'b_f=',b_f
        write (*,*) 'c_f=',c_f
        write (*,*) 'd_f=',d_f
        write (*,*) 'a_b=',a_b
        write (*,*) 'b_b=',b_b
        write (*,*) 'c_b=',c_b
        write (*,*) 'd_b=',d_b
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Div B conservation check:'
        write (*,*) 'abs(Det(D_f))*abs(br0_f)/bmag = ',
     &              abs(a_f*d_f-b_f*c_f)*abs(br0_f)/bmag
        write (*,*) 'abs(Det(D_b))*abs(br0_b)/bmag = ',
     &              abs(a_b*d_b-b_b*c_b)*abs(br0_b)/bmag
      end if
c
c ****** Get the value of Q.
c
      nsq= (a_f*d_b-b_f*c_b)**2
     &    +(a_b*b_f-a_f*b_b)**2
     &    +(d_b*c_f-d_f*c_b)**2
     &    +(a_b*d_f-b_b*c_f)**2
c
      q=nsq*abs(br0_f*br0_b)/bmag**2
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Final Q value:'
        write (*,*) 'q=',q
      end if
c
      valid=.true.
c
c ****** Deallocate the storage for the field line trace buffer
c ****** if it was used.
c
      if (save_trace_points) then
        call deallocate_trajectory_buffer (xt)
      end if
c
      return
      end
c#######################################################################
      subroutine normalized_cross_product (a,b,c)
c
c-----------------------------------------------------------------------
c
c ****** Return the unit vector C = (A x B)/|A x B|.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: a,b,c
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: cnorm
c
c-----------------------------------------------------------------------
c
c ****** Set C to the cross-product of A and B.
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
c
c ****** Normalize C to unit length.
c
      cnorm=sqrt(c(1)**2+c(2)**2+c(3)**2)
c    
      if (cnorm.ne.0.) then
        c=c/cnorm
      end if
c
      return
      end
c#######################################################################
      subroutine write_trace (fname,xt)
c
c-----------------------------------------------------------------------
c
c ****** Write the field line trace in structure XT to the
c ****** text file named FNAME.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      character, parameter :: TAB=achar(9)
c
c-----------------------------------------------------------------------
c
      integer :: i,ierr
c
c-----------------------------------------------------------------------
c
      call ffopen (1,fname,'rw',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_TRACE:'
        write (*,*) '### Could not create a field line output'//
     &              ' file.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
c
c ****** Write the output coordinates.
c
      write (1,'(5a)') 'r',TAB,'t',TAB,'p'
      do i=1,xt%npts
        write (1,'(3(1pe23.16,a))') xt%x(1)%f(i),TAB,
     &                              xt%x(2)%f(i),TAB,
     &                              xt%x(3)%f(i)
      enddo
c
      close (1)
c
      return
      end
c#######################################################################
      subroutine tracefl (b,ds,s0,s1,bs0,bs1,s,
     &                    traced_to_r_boundary,xt)
c
c-----------------------------------------------------------------------
c
c ****** Trace a magnetic field line.
c
c-----------------------------------------------------------------------
c
c ****** The 3D magnetic field is specified by structure B.
c ****** The structure DS has the field line integration parameters,
c ****** and S0 contains the spherical coordinates of the launch
c ****** point.
c
c ****** TRACED_TO_R_BOUNDARY=.T. is returned if the field line
c ****** was traced all the way to a radial domain boundary, in
c ****** which case S1 contains the spherical coordinates of the
c ****** final location, and S has the traced field line length.
c
c ****** Otherwise, TRACED_TO_R_BOUNDARY=.F. is returned, and
c ****** S and S1 do not necessarily have valid values.
c ****** This also occurs if |B|=0 is encountered during the trace.
c
c ****** The magnetic field vectors in spherical coordinates
c ****** at the starting and ending footpoints, respectively,
c ****** are returned in BS0 and BS1.
c
c ****** If the field line trace is needed, pass in the optional
c ****** field line trace buffer XT.  The buffer XT needs to be
c ****** allocated prior to being passed in to this routine.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use step_size_stats
      use debug
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(flparam) :: ds
      real(r_typ), dimension(3) :: s0,s1
      real(r_typ), dimension(3) :: bs0,bs1
      real(r_typ) :: s
      logical :: traced_to_r_boundary
      type(traj), optional :: xt
c
      intent(in) :: b,ds,s0
      intent(out) :: s1,bs0,bs1,s,traced_to_r_boundary
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      logical :: store_trace
      real(r_typ) :: ds0,dss,dsss,frac,dsmult_corrector
      logical :: done_tracing,first,nullb
      integer :: idir0,i,n,ntry
      type (csvec) :: x,xp,xo,bv,bhat1,bhat2
      type (inout) :: outside
      integer :: ierr
      real(r_typ) :: arc_length_remaining,ds_b
      type(flparam) :: current_ds
c
      integer :: local_stat_n
      real(r_typ) :: local_stat_ds_sum,
     &               local_stat_ds_min,
     &               local_stat_ds_max
c
c-----------------------------------------------------------------------
c
      logical, external :: outside_domain
c
c-----------------------------------------------------------------------
c
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) '### COMMENT from TRACEFL:'
        write (*,*) '### Starting a new trace:'
        write (*,*) 'S0 = ',S0
      end if
c
c ****** Initializations.
c
      store_trace=present(xt)
c
      current_ds=ds
      ntry=1
      first=.true.
      s1=s0
c
      if (gather_stats) then
        local_stat_n=0
        local_stat_ds_sum=0.
        local_stat_ds_min=huge(one)
        local_stat_ds_max=0.
      end if
c
c ****** If storage for the field line trace was requested,
c ****** initialize the field line trace buffer.
c
      if (store_trace) xt%npts=0
c
  100 continue
c
c-----------------------------------------------------------------------
c ****** Trace a field line starting at launch point S0.
c-----------------------------------------------------------------------
c
      done_tracing=.false.
c
c ****** If things go wrong, return TRACED_TO_R_BOUNDARY=.F..
c
      traced_to_r_boundary=.false.
c
      x%s=s0
      call sph_to_cart (x)
      s=0.
      n=1
      if (store_trace) call add_trajectory_point (xt,x%s)
c
c ****** Set the starting step size.
c
      ds0=current_ds%min
c
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) 'Start of trace:'
        write (*,*) 'X = ',x%s
      end if
c
c ****** If the initial point is outside the domain, stop tracing.
c
      if (outside_domain(b,x,outside)) go to 500
c
  200 continue
c
c-----------------------------------------------------------------------
c ****** Trace to the next step.
c-----------------------------------------------------------------------
c
c ****** Set the tracing direction.
c
      call getb (b,x,bv)
c
c ****** On the first time in, set the tracing direction in
c ****** variable IDIR0 based on that specified by DS%DIRECTION
c ****** and DS%DIRECTION_IS_ALONG_B.  Store the magnetic field at
c ****** the initial location (in spherical coordinates) in BS0.
c
c ****** When DS%DIRECTION_IS_ALONG_B=.F., DS%DIRECTION=1 traces
c ****** in the direction of increasing r, and DS%DIRECTION=-1
c ****** traces in the opposite direction.
c
c ****** When DS%DIRECTION_IS_ALONG_B=.T., DS%DIRECTION=1 traces
c ****** in the direction of the magnetic field vector, and
c ****** DS%DIRECTION=-1 traces in the opposite direction.
c
      if (first) then
        first=.false.
        bs0=bv%s
        if (ds%direction.gt.0) then
          idir0=1
        else
          idir0=-1
        end if
        if (.not.ds%direction_is_along_b) then
          if (bv%s(1).lt.0.) idir0=-idir0
        end if
      end if
c
      call normalize_v (bv,nullb)
      if (nullb) then
        write (*,*)
        write (*,*) '### WARNING from TRACEFL:'
        write (*,*) '### The trace encountered a null point (B = 0).'
        write (*,*) '### This occurred at the start of the trace.'
        write (*,*) '### Abandoning the trace ...'
        write (*,*) 'Location (r,t,p) = ',x%s
        go to 500
      end if
c
      if (n.gt.1) then
        ds_b=abs(dsss)
        bhat1=bhat2
      end if
      bhat2=bv
c
c ****** If a variable step size is being used, set the step
c ****** size for the next step.
c
c ****** If this is a repeat trace (of a short field line,
c ****** NTRY.gt.1), then leave the step size at the minimum
c ****** value until the number of points exceeds the minimum
c ****** allowed number.
c
      if (ds%variable.and.n.gt.1) then
        if (.not.(ntry.gt.1.and.n.le.ds%short_fl_min_points)) then
          call get_ds (b,x,bhat1,bhat2,ds_b,current_ds,ds0)
        end if
      end if
c
c ****** Check to see if this trace segment ends the trace
c ****** (i.e., exceeds the arc length specified, DS%LMAX).
c
      arc_length_remaining=ds%lmax-s
c
      if (ds0.ge.arc_length_remaining) then
        dss=arc_length_remaining
        done_tracing=.true.
      else
        dss=ds0
      end if
c
      if (dss.le.0.) go to 500
c
c ****** Gather step size statistics.
c
      if (gather_stats) then
        local_stat_n=local_stat_n+1
        local_stat_ds_sum=local_stat_ds_sum+abs(ds0)
        local_stat_ds_min=min(local_stat_ds_min,abs(ds0))
        local_stat_ds_max=max(local_stat_ds_max,abs(ds0))
      end if
c
c-----------------------------------------------------------------------
c ****** Predictor.
c-----------------------------------------------------------------------
c
c ****** Advance for a half-step to achieve second-order accuracy.
c
      dsmult_corrector=one
      xo=x
      xp=x
c
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) 'Predictor:'
        write (*,*) 'N = ',n
        write (*,*) 'BV = ',bv%s
      end if
c
      dsss=half*idir0*dss
      call advance (xp,bv,dsss)
c
      if (debug_level.ge.4) then
        write (*,*) 'DSSS = ',dsss
        write (*,*) 'XP = ',xp%s
      end if
c
c ****** Check if the field line has exited the domain.
c
      if (outside_domain(b,xp,outside)) then
c
        if (outside%t.or.outside%p) then
c
c ****** The field line has crossed the theta or phi boundaries;
c ****** stop tracing.
c
          go to 500
c
        else if (outside%r) then
c
c ****** The field line has crossed the r boundary.
c
          if (debug_level.ge.4) then
            write (*,*) '### Predictor: Outside r domain:'
          end if
c
c ****** Get the clip fraction, FRAC, that clips the segment
c ****** to the r boundary.
c
          call get_r_clip_fraction (b,xo,xp,outside%r0,frac,ierr)
c
          if (debug_level.ge.4) then
            write (*,*) 'After GET_R_CLIP_FRACTION (predictor):'
            write (*,*) 'IERR = ',ierr
            write (*,*) 'FRAC = ',frac
          end if
c
c ****** If there is an error in getting the clip fraction,
c ****** stop tracing.  This should never happen: write
c ****** detailed debugging information.
c
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ANOMALY in TRACEFL:'
            write (*,*) '### Predictor, 1st step:'
            write (*,*) '### Could not get the clip fraction.'
            write (*,*)
            write (*,*) '### Debugging info:'
            write (*,*) 'S0 = ',s0
            write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',
     &                  ds%direction_is_along_b
            write (*,*) 'DS%DIRECTION = ',ds%direction
            write (*,*) 'DS%MIN = ',current_ds%min
            write (*,*) 'DS%MAX = ',current_ds%max
            write (*,*) 'DSSS = ',dsss
            write (*,*) 'BV = ',bv
            write (*,*) 'XO = ',xo
            write (*,*) 'XP = ',xp
            go to 500
          end if
c
          if (frac.eq.0..and.n.eq.1) then
c
c ****** The initial point is exactly on the boundary, and is
c ****** being traced out of the boundary.  We are done.
c
            traced_to_r_boundary=.true.
            done_tracing=.true.
            go to 500
          end if
c
          if (frac.le.ds%predictor_min_clip_fraction) then
c
c ****** The starting point is close to the boundary.  Clip the
c ****** final point to the radial boundary and stop tracing.
c
            call clip_to_r_boundary (b,xo,xp,outside%r0,frac)
c
            if (debug_level.ge.4) then
              write (*,*) 'After CLIP_TO_R_BOUNDARY (predictor):'
              write (*,*) 'XP = ',xp%s
            end if
c
            x=xp
            dsss=frac*dsss
            traced_to_r_boundary=.true.
            done_tracing=.true.
            go to 300
c
          end if
c
c ****** If the starting point is not close enough to the radial
c ****** boundary, predict again half way to the radial boundary
c ****** to achieve second-order accuracy in the corrector.
c
          dsss=half*frac*dsss
          xp=x
          call advance (xp,bv,dsss)
c
          if (debug_level.ge.4) then
            write (*,*) 'After 2nd predictor:'
            write (*,*) 'XP = ',xp%s
          end if
c
c ****** Check if the predicted point has exited the domain.
c ****** This should never happen: write detailed debugging
c ****** information.
c
          if (outside_domain(b,xp,outside)) then
            write (*,*)
            write (*,*) '### ANOMALY in TRACEFL:'
            write (*,*) '### Predictor, 2nd step:'
            write (*,*) '### Point is outside the boundary.'
            write (*,*)
            write (*,*) '### Debugging info:'
            write (*,*) 'S0 = ',s0
            write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',
     &                  ds%direction_is_along_b
            write (*,*) 'DS%DIRECTION = ',ds%direction
            write (*,*) 'DS%MIN = ',current_ds%min
            write (*,*) 'DS%MAX = ',current_ds%max
            write (*,*) 'DSSS = ',dsss
            write (*,*) 'BV = ',bv
            write (*,*) 'XO = ',xo
            write (*,*) 'XP = ',xp
            go to 500
          end if
c
c ****** Reduce the step size for the corrector.  This minimizes the
c ****** chance of "going through the radial boundary and back into
c ****** the domain again", which can happen when the step size
c ****** is too big.
c
          dsmult_corrector=half*frac
c
        end if
c
      end if
c
c-----------------------------------------------------------------------
c ****** Corrector.
c-----------------------------------------------------------------------
c
      call getb (b,xp,bv)
      call normalize_v (bv,nullb)
      if (nullb) then
        write (*,*)
        write (*,*) '### WARNING from TRACEFL:'
        write (*,*) '### The trace encountered a null point (B = 0).'
        write (*,*) '### This occurred during the corrector.'
        write (*,*) '### Abandoning the trace ...'
        write (*,*) 'Location (r,t,p) = ',xp%s
        go to 500
      end if
c
      dsss=idir0*dsmult_corrector*dss
      call advance (x,bv,dsss)
c
      if (debug_level.ge.4) then
        write (*,*) 'After corrector advance:'
        write (*,*) 'BV = ',bv%s
        write (*,*) 'DSSS = ',dsss
        write (*,*) 'X = ',x%s
      end if
c
      if (outside_domain(b,x,outside)) then
        if (outside%t.or.outside%p) then
          go to 500
        else if (outside%r) then
c
          if (debug_level.ge.4) then
            write (*,*) '### Corrector: Outside r domain:'
          end if
c
          call get_r_clip_fraction (b,xo,x,outside%r0,frac,ierr)
c
c ****** If there is an error in getting the clip fraction,
c ****** stop tracing.  This should never happen: write
c ****** detailed debugging information.
c
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ANOMALY in TRACEFL:'
            write (*,*) '### Corrector step:'
            write (*,*) '### Could not get the clip fraction.'
            write (*,*)
            write (*,*) '### Debugging info:'
            write (*,*) 'S0 = ',s0
            write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',
     &                  ds%direction_is_along_b
            write (*,*) 'DS%DIRECTION = ',ds%direction
            write (*,*) 'DS%MIN = ',current_ds%min
            write (*,*) 'DS%MAX = ',current_ds%max
            write (*,*) 'DSSS = ',dsss
            write (*,*) 'BV = ',bv
            write (*,*) 'XO = ',xo
            write (*,*) 'X  = ',x
            go to 500
          end if
c
          if (debug_level.ge.4) then
            write (*,*) 'After GET_R_CLIP_FRACTION (corrector):'
            write (*,*) 'IERR = ',ierr
            write (*,*) 'FRAC = ',frac
          end if
c
          done_tracing=.true.
          traced_to_r_boundary=.true.
c
          if (frac.eq.0.) then
c
c ****** This should only happen when the previous point was exactly
c ****** on the boundary, and is being traced out of the domain.
c ****** In this case, do not store this point: we are done.
c
            x=xo
c
            if (debug_level.ge.4) then
              write (*,*) 'The trace went from the boundary to'//
     &                    ' the outside (corrector):'
              write (*,*) 'X = ',x%s
            end if
c
            go to 400
c
          else
c
            call clip_to_r_boundary (b,xo,x,outside%r0,frac)
            dsss=frac*dsss
c
            if (debug_level.ge.4) then
              write (*,*) 'After CLIP_TO_R_BOUNDARY (corrector):'
              write (*,*) 'X = ',x%s
            end if
c
          end if
        end if
      end if
c
  300 continue
c
c ****** Increment the number of points and the arc length.
c
      n=n+1
      s=s+abs(dsss)
c
c ****** Add the current position to the field line buffer
c ****** if requested.
c
      if (store_trace) call add_trajectory_point (xt,x%s)
c
      if (.not.done_tracing) go to 200
c
  400 continue
c
c ****** Finished tracing the field line.
c
c ****** Check that the number of points in the field line trace
c ****** exceeds the minimum allowed, DS%SHORT_FL_MIN_POINTS.
c ****** If not, reduce the step size and retrace the field line,
c ****** up to a maximum of DS%SHORT_FL_MAX_TRIES times.
c
      if (n.lt.ds%short_fl_min_points) then
        if (ntry.le.ds%short_fl_max_tries) then
c
          if (debug_level.ge.4) then
            write (*,*) 'Short field line:'
            write (*,*) 'N = ',n
            write (*,*) 'CURRENT_DS%OVER_RC = ',current_ds%over_rc
            write (*,*) 'CURRENT_DS%MIN = ',current_ds%min
            write (*,*) 'CURRENT_DS%MAX = ',current_ds%max
          end if
c
          ntry=ntry+1
          current_ds%over_rc=current_ds%over_rc*
     &                       ds%short_fl_shrink_factor
          current_ds%min=current_ds%min*
     &                   ds%short_fl_shrink_factor
          current_ds%max=current_ds%max*
     &                   ds%short_fl_shrink_factor
c
          if (store_trace) then
            call deallocate_trajectory_buffer (xt)
            call allocate_trajectory_buffer (xt)
            xt%npts=0
          end if
c
          go to 100
c
        else
          write (*,*)
          write (*,*) '### WARNING from TRACEFL:'
          write (*,*) 'Short field line after max tries:'
          write (*,*) 'Number of points in field line = ',n
          write (*,*) 'Number of tries = ',ntry
          write (*,*) 'CURRENT_DS%OVER_RC = ',current_ds%over_rc
          write (*,*) 'CURRENT_DS%MIN = ',current_ds%min
          write (*,*) 'CURRENT_DS%MAX = ',current_ds%max
          write (*,*) 'S0 = ',s0
          write (*,*) 'S1 = ',x%s
        end if
      end if
c
  500 continue
c
c ****** Update the step size statistics.
c
      if (gather_stats) then
c$omp critical
        stat_n=stat_n+local_stat_n
        stat_ds_sum=stat_ds_sum+local_stat_ds_sum
        stat_ds_min=min(stat_ds_min,local_stat_ds_min)
        stat_ds_max=max(stat_ds_max,local_stat_ds_max)
c$omp end critical
      end if
c
c ****** Store the final location in S1.
c
      s1=x%s
c
c ****** If the field line was not traced to the radial boundary,
c ****** do not attempt to return the magnetic field vector
c ****** at the field line endpoint.
c
      if (.not.traced_to_r_boundary) then
        bs1=0.
        return
      end if
c
c ****** Store the magnetic field vector (in spherical coordinates)
c ****** at the field line endpoint in BS1.
c
      call getb (b,x,bv)
      bs1=bv%s
c
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) '### COMMENT from TRACEFL:'
        write (*,*) '### About to exit:'
        write (*,*) 'N = ',n
        write (*,*) 'S = ',s
        write (*,*) 'S1 = ',s1
        write (*,*) 'BS0 = ',bs0
        write (*,*) 'BS1 = ',bs1
      end if
c
      return
      end
c#######################################################################
      subroutine get_ds (b,x,v0,v1,ds_v,ds,deltas)
c
c-----------------------------------------------------------------------
c
c ****** Set the integration step size based on the radius of
c ****** curvature of the field line, as estimated from the unit
c ****** vectors V0 and V1, and the local mesh cell size of the
c ****** magnetic field B at X (if requested).
c
c ****** The vectors V0 and V1 are assumed to have been evaluated
c ****** along the field line a distance DS_V apart.
c
c ****** On input, DELTAS should have the present value of the
c ****** step size.  On return, DELTAS is overwritten by the new
c ****** step size.
c
c ****** It is assumed that DS_V and DELTAS are positive.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x
      type(csvec) :: v0,v1
      real(r_typ) :: ds_v
      type(flparam) :: ds
      real(r_typ) :: deltas
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: dv,factor,ds_mesh
c
c-----------------------------------------------------------------------
c
c ****** If DELTAS is zero, set it to the minimum value.
c
      if (deltas.eq.0.) then
        deltas=ds%min
        return
      end if
c
c ****** First, set the step size based on the local radius
c ****** of curvature of the field line.
c
c ****** Compute the factor by which DELTAS must be multiplied to
c ****** achieve the specified ratio of step size to radius
c ****** of curvature, DS%OVER_RC.
c
      dv=sqrt( (v0%c(1)-v1%c(1))**2
     &        +(v0%c(2)-v1%c(2))**2
     &        +(v0%c(3)-v1%c(3))**2)
c
      if (dv.ne.0.) then
        factor=ds_v*ds%over_rc/(dv*deltas)
      else
        factor=huge(dv)
      end if
c
c ****** Limit the change in DELTAS to the maximum permitted
c ****** change per step.
c
      factor=min(factor,ds%max_increase_factor)
      factor=max(factor,ds%max_decrease_factor)
c
c ****** Set the new step size.
c
      deltas=factor*deltas
c
c ****** Next, if requested, limit DELTAS by the local mesh
c ****** cell size of B at X.
c
c ****** Note that this only needs to be done if DELTAS is bigger
c ****** than DS%MIN, since only then is there a possibility of
c ****** reducing DELTAS further.
c
      if (deltas.gt.ds%min.and.ds%limit_by_local_mesh) then
c
c ****** Get the local mesh cell size at X.
c
        call get_local_mesh_size (b,x%s,ds_mesh)
c
c ****** Set DELTAS so that it does not exceed DS%LOCAL_MESH_FACTOR
c ****** times the local mesh size.
c
        deltas=min(deltas,ds%local_mesh_factor*ds_mesh)
c
      end if
c
c ****** Limit DELTAS by the maximum and mimimum allowed values.
c
      deltas=max(deltas,ds%min)
      deltas=min(deltas,ds%max)
c
      return
      end
c#######################################################################
      subroutine get_local_mesh_size (b,s,ds)
c
c-----------------------------------------------------------------------
c
c ****** Get the local mesh size DS from the magnetic field in
c ****** structure B at the spherical position S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use interp_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      real(r_typ), dimension(3) :: s
      real(r_typ) :: ds
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ar,at,ap,drv,dtv,dpv,stv
c
c-----------------------------------------------------------------------
c
c ****** Get the local size of the main mesh of B at the
c ****** specified point.
c
      call interp (b%nrs,b%rs,s(1),i,ip1,ar,b%rs_invtab)
      drv=(one-ar)*b%drs(i)+ar*b%drs(ip1)
c
      call interp (b%nts,b%ts,s(2),j,jp1,at,b%ts_invtab)
      dtv=(one-at)*b%dts(j)+at*b%dts(jp1)
      stv=(one-at)*b%sts(j)+at*b%sts(jp1)
c
      call interp (b%nps,b%ps,s(3),k,kp1,ap,b%ps_invtab)
      dpv=(one-ap)*b%dps(k)+ap*b%dps(kp1)
c
      ds=min(drv,s(1)*dtv,s(1)*stv*dpv)
c
      return
      end
c#######################################################################
      subroutine normalize_v (v,null)
c
c-----------------------------------------------------------------------
c
c ****** Normalize the vector V to return a unit vector along V.
c ****** If V has zero length, then set NULL=.T. and leave V
c ****** unchanged; otherwise, set NULL=.F..
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: v
      logical :: null
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: vmag
c
c-----------------------------------------------------------------------
c
c ****** Use the spherical representation to compute the norm.
c
      vmag=sqrt(v%s(1)**2+v%s(2)**2+v%s(3)**2)
c
      if (vmag.eq.0.) then
        null=.true.
      else
        null=.false.
        v%c=v%c/vmag
        v%s=v%s/vmag
      end if
c
      return
      end
c#######################################################################
      subroutine advance (x,v,ds)
c
c-----------------------------------------------------------------------
c
c ****** Advance the position vector X by the step DS using the
c ****** velocity vector V.
c
c ****** This routine updates both Cartesian and spherical
c ****** represenations of X in the dual represenations position
c ****** vector X.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: x,v
      real(r_typ) :: ds
c
c-----------------------------------------------------------------------
c
c ****** Advance the Cartesian position.
c
      x%c=x%c+ds*v%c
c
c ****** Transform the new Cartesian position to spherical
c ****** coordinates.
c
      call cart_to_sph (x)
c
      return
      end
c#######################################################################
      subroutine get_r_clip_fraction (b,x0,x1,outside_r0,frac,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Get the fraction FRAC (between 0 and 1) that expresses
c ****** the normalized distance between X0 and X1
c ****** corresponding to the location of the radial boundary.
c
c ****** The boundary location is obtained from the 3D magnetic
c ****** field in structure B.
c
c ****** It is assumed that X0 and X1 lie on different sides
c ****** of the r boundary specified by flag OUTSIDE_R0.
c
c ****** FRAC can be used to clip the position to the
c ****** radial boundary.
c
c ****** For a normal return, IERR=0 is returned when it was possible
c ****** to estimate FRAC; in the case of an inconsistency,
c ****** IERR=1 is returned and FRAC is invalid.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x0,x1
      logical :: outside_r0
      real(r_typ) :: frac
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rb,dr0,dr1,dssq,x0dotb
      real(r_typ) :: eps,term1,term2,ratio,disc
      integer :: ipm
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      if (outside_r0) then
c
c ****** The segment crossed the inner radial boundary.
c
        rb=b%lim0(1)
        ipm=-1
c
      else
c
c ****** The segment crossed the outer radial boundary.
c
        rb=b%lim1(1)
        ipm=1
c
      end if
c
c ****** Check that X0 and X1 are consistent with a radial
c ****** boundary crossing.
c
      dr0=rb-x0%s(1)
      dr1=rb-x1%s(1)
c
      if (dr0.eq.0..and.dr1.eq.0.) then
        ierr=1
        return
      end if
c
      if (dr0*dr1.gt.0.) then
        ierr=1
        return
      end if
c
c ****** Treat the special case when X0 is exactly on the boundary.
c
      if (dr0.eq.0.) then
        frac=0.
        return
      end if
c
      dssq= (x1%c(1)-x0%c(1))**2
     &     +(x1%c(2)-x0%c(2))**2
     &     +(x1%c(3)-x0%c(3))**2
c
      if (dssq.le.0.) then
        ierr=1
        return
      end if
c
      x0dotb= x0%c(1)*(x1%c(1)-x0%c(1))
     &       +x0%c(2)*(x1%c(2)-x0%c(2))
     &       +x0%c(3)*(x1%c(3)-x0%c(3))
c
c ****** Get the square root of the discriminant, expanding small
c ****** arguments for accuracy.
c
      eps=10*sqrt(spacing(one))
c
      term1=x0dotb**2
      term2=dssq*dr0*(two*x0%s(1)+dr0)
      if (term1+term2.lt.0.) then
        ierr=1
        return
      end if
      if (term1.ne.0.) then
        ratio=term2/term1
        if (abs(ratio).lt.eps) then
          disc=sqrt(term1)*(one+half*ratio)
        else
          disc=sqrt(term1+term2)
        end if
      else
        disc=sqrt(term1+term2)
      end if
c
      frac=(-x0dotb+ipm*disc)/dssq
c
      return
      end
c#######################################################################
      subroutine clip_to_r_boundary (b,x0,x1,outside_r0,frac)
c
c-----------------------------------------------------------------------
c
c ****** Reset the position X1 to the normalized distance FRAC
c ****** between X0 and X1.
c
c ****** The boundary location is obtained from the 3D magnetic
c ****** field in structure B.
c
c ****** The flag OUTSIDE_R0=.T. indicates that the lower
c ****** radial boundary was crossed in going from X0 to X1;
c ****** otherwise, the upper radial boundary was crossed.
c
c ****** In conjunction with GET_R_CLIP_FRACTION this routine
c ****** can be used to clip points that cross the radial boundaries
c ****** to the boundary.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x0,x1
      logical :: outside_r0
      real(r_typ) :: frac
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rval
c
c-----------------------------------------------------------------------
c
c ****** Reset the Cartesian position in X1.
c
      x1%c=(one-frac)*x0%c+frac*x1%c
c
c ****** Transform the new Cartesian position to spherical
c ****** coordinates.
c
      call cart_to_sph (x1)
c
c ****** Set the radius to the appropriate radial boundary
c ****** value (exactly).  This is done to take care of roundoff.
c
      if (outside_r0) then
        rval=b%lim0(1)
      else
        rval=b%lim1(1)
      end if
c
      x1%s(1)=rval
c
c ****** Transform the new spherical position to Cartesian
c ****** coordinates.
c
      call sph_to_cart (x1)
c
      return
      end
c#######################################################################
      function outside_domain (b,x,outside)
c
c-----------------------------------------------------------------------
c
c ****** If the spherical position in X lies outside the limits
c ****** of the domain, return a function result of .T.;
c ****** otherwise, return .F..
c
c ****** The detailed in/out location for each coordinate
c ****** is set in structure OUTSIDE.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x
      type(inout) :: outside
      logical :: outside_domain
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: pv
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Get the value of phi in the main interval [0,2*pi].
c
      pv=fold(zero,twopi,x%s(3))
c
      outside%r0=x%s(1).lt.b%lim0(1)
      outside%t0=x%s(2).lt.b%lim0(2)
      outside%p0=pv.lt.b%lim0(3)
      outside%r1=x%s(1).gt.b%lim1(1)
      outside%t1=x%s(2).gt.b%lim1(2)
      outside%p1=pv.gt.b%lim1(3)
c
      outside%r=outside%r0.or.outside%r1
      outside%t=outside%t0.or.outside%t1
      outside%p=outside%p0.or.outside%p1
c
      outside%domain=outside%r.or.outside%t.or.outside%p
c
      outside_domain=outside%domain
c
      return
      end
c#######################################################################
      function fold (x0,x1,x)
c
c-----------------------------------------------------------------------
c
c ****** "Fold" X into the periodic interval [X0,X1].
c
c ****** On return, X is such that X0.le.X.lt.X1.
c
c-----------------------------------------------------------------------
c
c ****** It is assumed that X0 does not equal X1, as is physically
c ****** necessary.  If X0 and X1 are equal, the routine just
c ****** returns with FOLD=X.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: fold
      real(r_typ) :: x0,x1,x
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xl
c
c-----------------------------------------------------------------------
c
      fold=x
c
      if (x0.eq.x1) return
c
      xl=x1-x0
c
      fold=mod(x-x0,xl)+x0
c
      if (fold.lt.x0) fold=fold+xl
      if (fold.ge.x1) fold=fold-xl
c
      return
      end
c#######################################################################
      subroutine cart_to_sph (x)
c
c-----------------------------------------------------------------------
c
c ****** Update the dual representation position vector X so that
c ****** the Cartesian position corresponds to the spherical position.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: x
c
c-----------------------------------------------------------------------
c
      call c2s (x%c,x%s)
c
      return
      end
c#######################################################################
      subroutine sph_to_cart (x)
c
c-----------------------------------------------------------------------
c
c ****** Update the dual representation position vector X so that
c ****** the spherical position corresponds to the Cartesian position.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: x
c
c-----------------------------------------------------------------------
c
      call s2c (x%s,x%c)
c
      return
      end
c#######################################################################
      subroutine c2s (x,s)
c
c-----------------------------------------------------------------------
c
c ****** Convert the vector X = (x,y,z) from Cartesian coordinates
c ****** to spherical coordinates S = (r,t,p).
c
c ****** This routine returns T and P in radians, in the
c ****** following range:
c
c          0. .le. t .le. pi
c          0. .le. p .lt. 2.*pi
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: x
      real(r_typ), dimension(3) :: s
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r,t,p
c
c-----------------------------------------------------------------------
c
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
c
      if (r.eq.0.) then
        t=0.
      else
        t=acos(x(3)/r)
      end if
c
      if (x(1).eq.0.) then
        if (x(2).ge.0.) then
          p= halfpi
        else
          p=-halfpi
        end if
      else
        p=atan2(x(2),x(1))
      end if
      if (p.lt.0.) p=p+twopi
c
      s(1)=r
      s(2)=t
      s(3)=p
c
      return
      end
c#######################################################################
      subroutine s2c (s,x)
c
c-----------------------------------------------------------------------
c
c ****** Convert the vector S = (r,t,p) from spherical coordinates
c ****** to Cartesian coordinates X = (x,y,z).
c
c ****** This routine assumes that T and P are in radians.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: s
      real(r_typ), dimension(3) :: x
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rst
c
c-----------------------------------------------------------------------
c
      rst=s(1)*sin(s(2))
      x(1)=rst*cos(s(3))
      x(2)=rst*sin(s(3))
      x(3)=s(1)*cos(s(2))
c
      return
      end
c#######################################################################
      subroutine cv_to_sv (s,cv,sv)
c
c-----------------------------------------------------------------------
c
c ****** Convert the Cartesian vector CV to spherical vector SV
c ****** at spherical position S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: s,cv,sv
      intent(in) :: s,cv
      intent(out) :: sv
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: st,ct,sp,cp
c
c-----------------------------------------------------------------------
c
      st=sin(s(2))
      ct=cos(s(2))
      sp=sin(s(3))
      cp=cos(s(3))
c
      sv(1)= cv(1)*st*cp+cv(2)*st*sp+cv(3)*ct
      sv(2)= cv(1)*ct*cp+cv(2)*ct*sp-cv(3)*st
      sv(3)=-cv(1)*sp   +cv(2)*cp
c
      return
      end
c#######################################################################
      subroutine sv_to_cv (s,sv,cv)
c
c-----------------------------------------------------------------------
c
c ****** Convert the spherical vector SV to Cartesian vector CV
c ****** at spherical position S.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: s,sv,cv
      intent(in) :: s,sv
      intent(out) :: cv
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: st,ct,sp,cp
c
c-----------------------------------------------------------------------
c
      st=sin(s(2))
      ct=cos(s(2))
      sp=sin(s(3))
      cp=cos(s(3))
c
      cv(1)= sv(1)*st*cp+sv(2)*ct*cp-sv(3)*sp
      cv(2)= sv(1)*st*sp+sv(2)*ct*sp+sv(3)*cp
      cv(3)= sv(1)*ct   -sv(2)*st
c
      return
      end
c#######################################################################
      subroutine getb (b,x,bv)
c
c-----------------------------------------------------------------------
c
c ****** Get the interpolated magnetic field vector BV at the
c ****** dual representation position vector X from the
c ****** magnetic field vector in structure B.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use evaluate_spline_3d_interface
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x,bv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: pv
      real(r_typ) :: br,bt,bp,st,ct,sp,cp
      real(r_typ), dimension(3) :: xs
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Get the value of phi in the main interval [0,2*pi].
c
      pv=fold(zero,twopi,x%s(3))
c
c ****** Get Br, Bt, and Bp.
c
      if (b%cubic) then
c
c ****** Cubic spline interpolation.
c
        br=evaluate_spline_3d(b%spl%r,x%s(1),x%s(2),pv,
     &                        b%inv(1)%c(1),
     &                        b%inv(1)%c(2),
     &                        b%inv(1)%c(3))
c
        bt=evaluate_spline_3d(b%spl%t,x%s(1),x%s(2),pv,
     &                        b%inv(2)%c(1),
     &                        b%inv(2)%c(2),
     &                        b%inv(2)%c(3))
c
        bp=evaluate_spline_3d(b%spl%p,x%s(1),x%s(2),pv,
     &                        b%inv(3)%c(1),
     &                        b%inv(3)%c(2),
     &                        b%inv(3)%c(3))
c
      else
c
c ****** Linear interpolation.
c
        call interp_3d (b%r%dims(1),b%r%dims(2),b%r%dims(3),
     &                  b%r%scales(1)%f,
     &                  b%r%scales(2)%f,
     &                  b%r%scales(3)%f,
     &                  b%inv(1),
     &                  b%r%f,x%s(1),x%s(2),pv,br)
c
        call interp_3d (b%t%dims(1),b%t%dims(2),b%t%dims(3),
     &                  b%t%scales(1)%f,
     &                  b%t%scales(2)%f,
     &                  b%t%scales(3)%f,
     &                  b%inv(2),
     &                  b%t%f,x%s(1),x%s(2),pv,bt)
c
        call interp_3d (b%p%dims(1),b%p%dims(2),b%p%dims(3),
     &                  b%p%scales(1)%f,
     &                  b%p%scales(2)%f,
     &                  b%p%scales(3)%f,
     &                  b%inv(3),
     &                  b%p%f,x%s(1),x%s(2),pv,bp)
c
      end if
c
      bv%s(1)=br
      bv%s(2)=bt
      bv%s(3)=bp
c
  100 continue
c
c ****** Transform the spherical components of the magnetic field
c ****** to the Cartesian components.
c
      st=sin(x%s(2))
      ct=cos(x%s(2))
      sp=sin(pv)
      cp=cos(pv)
      bv%c(1)=bv%s(1)*st*cp+bv%s(2)*ct*cp-bv%s(3)*sp
      bv%c(2)=bv%s(1)*st*sp+bv%s(2)*ct*sp+bv%s(3)*cp
      bv%c(3)=bv%s(1)*ct   -bv%s(2)*st
c
      return
      end
c#######################################################################
      subroutine interp_3d (nx,ny,nz,x,y,z,inv,f,xv,yv,zv,fv)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the value of the 3D field FV at (XV,YV,ZV) from
c ****** array F(NX,NY,NZ), defined on the mesh X(NX) x Y(NY) x Z(NZ).
c ****** The structure INV holds the inverse interpolation tables.
c
c ****** Note that if the point (XV,YV,ZV) is outside the bounds of
c ****** the X x Y x Z mesh, FV=0. is returned.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use interp_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx,ny,nz
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      type(vtab) :: inv
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ) :: xv,yv,zv,fv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ax,ay,az
c
c-----------------------------------------------------------------------
c
c ****** If the point is outside the data limits, return a
c ****** zero value.
c
      if (xv.lt.x(1).or.xv.gt.x(nx).or.
     &    yv.lt.y(1).or.yv.gt.y(ny).or.
     &    zv.lt.z(1).or.zv.gt.z(nz)) then
        fv=0.
        return
      end if
c
      call interp (nx,x,xv,i,ip1,ax,inv%c(1))
      call interp (ny,y,yv,j,jp1,ay,inv%c(2))
      call interp (nz,z,zv,k,kp1,az,inv%c(3))
c
      fv= (one-ax)*( (one-ay)*( (one-az)*f(i  ,j  ,k  )
     &                         +     az *f(i  ,j  ,kp1))
     &              +     ay *( (one-az)*f(i  ,jp1,k  )
     &                         +     az *f(i  ,jp1,kp1)))
     &   +     ax *( (one-ay)*( (one-az)*f(ip1,j  ,k  )
     &                         +     az *f(ip1,j  ,kp1))
     &              +     ay *( (one-az)*f(ip1,jp1,k  )
     &                         +     az *f(ip1,jp1,kp1)))
c
      return
      end
c#######################################################################
      subroutine allocate_trajectory_buffer (xt)
c
c-----------------------------------------------------------------------
c
c ****** Allocate the trajectory buffer XT.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      xt%ndim=3
      xt%size=xt%initial_size
c
c ****** Allocate storage for the trajectory buffer.
c
      allocate (xt%x(xt%ndim))
c
      do i=1,xt%ndim
        allocate (xt%x(i)%f(xt%size))
      enddo
c
c ****** Initialize the current number of points in the
c ****** trajectory buffer.
c
      xt%npts=0
c
      return
      end
c#######################################################################
      subroutine deallocate_trajectory_buffer (xt)
c
c-----------------------------------------------------------------------
c
c ****** Deallocate the trajectory buffer XT.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      if (.not.associated(xt%x)) return
c
      do i=1,xt%ndim
        if (associated(xt%x(i)%f)) then
          deallocate (xt%x(i)%f)
        end if
      enddo
c
      deallocate (xt%x)
c
      xt%npts=0
c
      return
      end
c#######################################################################
      subroutine add_trajectory_point (xt,x)
c
c-----------------------------------------------------------------------
c
c ****** Add the position vector X to the trajectory buffer XT.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
      real(r_typ), dimension(xt%ndim) :: x
c
c-----------------------------------------------------------------------
c
      integer :: i,n
c
c-----------------------------------------------------------------------
c
c ****** Increment the number of points in the trajectory.
c ****** If the buffer is full, expand it.
c
      n=xt%npts
      n=n+1
      if (n.gt.xt%size) call expand_trajectory_buffer (xt)
c
c ****** Add the point to the buffer.
c
      do i=1,xt%ndim
        xt%x(i)%f(n)=x(i)
      enddo
c
      xt%npts=n
c
      return
      end
c#######################################################################
      subroutine expand_trajectory_buffer (xt)
c
c-----------------------------------------------------------------------
c
c ****** Expand the trajectory buffer XT by doubling the number
c ****** of points in the buffer.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use debug
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(:), pointer :: f
      integer :: i,n
c
c-----------------------------------------------------------------------
c
c ****** Double the current buffer size, and copy the current
c ****** contents to the expanded buffer.
c
      n=xt%size
      do i=1,xt%ndim
        allocate (f(2*n))
        f(1:n)=xt%x(i)%f(1:n)
        deallocate (xt%x(i)%f)
        xt%x(i)%f=>f
      enddo
      xt%size=2*n
      if (debug_level.ge.5) then
        write (*,*) 'Expanded a trajectory buffer to ',xt%size
      end if
c
      return
      end
c#######################################################################
      subroutine build_inverse_tables (s,inv)
c
c-----------------------------------------------------------------------
c
c ****** Build the inverse interpolation tables INV for the SDS
c ****** in structure S.
c
c ****** These arrays are used to to increase the efficiency
c ****** of interpolation lookups.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
      type(vtab) :: inv
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
c ****** Use a number of points for the inverse interpolation table
c ****** equal to the number in the original scale.
c
      do i=1,s%ndim
        inv%c(i)%n=s%dims(i)
        allocate (inv%c(i)%f(inv%c(i)%n))
        call getinv (s%scales(i)%f,s%dims(i),inv%c(i))
      enddo
c
      return
      end
c#######################################################################
      subroutine getinv (x,n,tab)
c
c-----------------------------------------------------------------------
c
c ****** Build an inverse interpolation table to increase the
c ****** efficiency of table look-up in a nonuniform mesh.
c
c ****** On input, the table X(N) is specified, together with the
c ****** number of points to use in the inverse interpolation
c ****** table, NU.
c
c ****** The output is a structure TAB with the inverse interpolation
c ****** table.  This structure has the following components:
c
c ******    N:  the number of points in the table;
c ******    D:  the inverse of the uniform table spacing;
c ******    F:  the inverse interpolation table.
c
c-----------------------------------------------------------------------
c
      use number_types
      use invint_def
      use interp_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      type(itab) :: tab
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,k,ip1
      real(r_typ) :: dx,xv,alpha,en
c
c-----------------------------------------------------------------------
c
c ****** Check that the number of points is valid.
c
      if (tab%n.le.1) then
        write (*,*)
        write (*,*) '### ERROR in GETINV:'
        write (*,*) '### Invalid number of points specified'//
     &              ' for the inverse interpolation table.'
        write (*,*)
        write (*,*) 'Number of points = ',tab%n
        call exit (1)
      end if
c
c ****** Set the uniform interval to be used in the inverse
c ****** interpolation.
c
      dx=(x(n)-x(1))/(tab%n-one)
c
      if (dx.le.0.) then
        write (*,*)
        write (*,*) '### ERROR in GETINV:'
        write (*,*) '### Invalid interval for the inverse'//
     &              ' interpolation table.'
        write (*,*)
        write (*,*) 'Interval = ',dx
        call exit (1)
      end if
c
      tab%d=one/dx
c
c ****** Build the inverse interpolation table.
c
      en=n
c
      do k=1,tab%n
        xv=x(1)+(k-one)*dx
        xv=max(xv,x(1))
        xv=min(xv,x(n))
        call interp (n,x,xv,i,ip1,alpha)
        tab%f(k)=i+alpha
        tab%f(k)=max(tab%f(k),one)
        tab%f(k)=min(tab%f(k),en)
      enddo
c
      return
      end
c#######################################################################
      subroutine interp (n,x,xv,i,ip1,alpha,tab)
c
c-----------------------------------------------------------------------
c
c ****** Find the interval I in table X(i), i=1,2,...,N, that encloses
c ****** the value XV, such that X(I).le.XV.le.X(I+1).
c ****** For the special case when N=1, XV must equal X(1) exactly.
c
c ****** This routine uses LOCATE_INTERVAL (from the SPLINE library)
c ****** to do the actual work.  If the interval is not found
c ****** LOCATE_INTERVAL terminates with an error.
c
c ****** This routine does not do the actual interpolation.  However,
c ****** the returned values of I, IP1 (which generally equals I+1),
c ****** and ALPHA can be used to get the interpolant.
c
c ****** The optional inverse interpolation table, TAB, can be
c ****** supplied to improve the efficiency of the search.
c
c-----------------------------------------------------------------------
c
      use number_types
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i
      integer :: ip1
      real(r_typ) :: alpha
      type(itab), optional :: tab
      intent(in) :: n,x,xv,tab
      intent(out) :: i,ip1,alpha
c
c-----------------------------------------------------------------------
c
      if (present(tab)) then
        i=locate_interval(n,x,xv,tab)
      else
        i=locate_interval(n,x,xv)
      end if
c
      if (n.eq.1) then
        ip1=1
        alpha=0.
      else
        ip1=i+1
        if (x(i).eq.x(i+1)) then
          alpha=0.
        else
          alpha=(xv-x(i))/(x(i+1)-x(i))
        end if
      end if
c
      return
      end
c#######################################################################
      subroutine set_parameters
c
c-----------------------------------------------------------------------
c
c ****** Set parameters from the command line arguments.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use syntax
      use paragraph_def
      use get_usage_line_interface
      use print_par_interface
      use delete_par_interface
      use params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Storage the for usage line.
c
      type(paragraph), pointer :: usage
c
c ****** Storage for the error message.
c
      character(72) :: errmsg
c
c-----------------------------------------------------------------------
c
      integer :: ierr
      character(256) :: arg
      logical :: set
c
c-----------------------------------------------------------------------
c
c ****** Define the syntax.
c
      call defarg (GROUP_K ,'-v',' ',' ')
      call defarg (GROUP_A ,'infile',' ',' ')
c
c ****** Parse the command line.
c
      call parse (errmsg,ierr)
c
      if (ierr.ne.0) then
c
        write (*,*)
        write (*,*) '### ',cname,' Version ',cvers,' of ',cdate,'.'
        write (*,*) '### Calculate the field line mapping'//
     &              ' for MAS code runs.'
c
        if (ierr.gt.1) then
          write (*,*)
          write (*,*) errmsg
        end if
c
c ****** Print the usage line.
c
        call get_usage_line (usage)
c
        write (*,*)
        write (*,*) 'Usage:'
        write (*,*)
c
        call print_par (usage)
c
        write (*,*)
        write (*,*) 'Read the parameters from input file <infile>.'
        call delete_par (usage)
c
        call exit (1)
c
      end if
c
c ****** Set the parameters.
c
c ****** Verbose flag.
c
      call fetcharg ('-v',set,arg)
      verbose=set
c
c ****** Input file name.
c
      call fetcharg ('infile',set,arg)
      infile=trim(arg)
c
      return
      end
