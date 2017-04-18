c
c-----------------------------------------------------------------------
c
c ****** Source to build the SDS library.
c ****** These routines are used by Zoran Mikic's tools.
c
c-----------------------------------------------------------------------
c
c        07/29/2003, ZM, Version 1.00:
c
c         - Original version of the SDS library.
c           This library was put together to facilitate the
c           development of ZM's tools.
c           It includes the new read/write routines for
c           scientific data sets (both text and HDF format).
c           The code was cleaned up to use standard FORTRAN90.
c
c        02/20/2004, ZM, Version 1.01:
c
c         - Added the ability to specify the format in writing
c           floating point numbers in routine WRFP.  This is used
c           in writing SDS text files using routine WRTXT.
c           For 32-bit data, WRTXT specifies 7 digits of precision,
c           whereas for 64-bit data, 15 digits of precision are
c           used.
c
c        04/02/2005, ZM, Version 1.02:
c
c         - Added a call to DFSDclear() in routine WRHDF.  This
c           initializes the SDS interface to the default state
c           for each new file.  This is needed to prevent settings
c           from previous files from interfering with each other.
c
c        06/16/2006, ZM, Version 1.03:
c
c         - Fixed some pointer allocation and manipulation issues in
c           the HDF read and write routines to make them obey
c           the FORTRAN 90 standard.  They were misbehaving on
c           the IFORT compiler.
c
c        02/24/2009, ZM, Version 1.04:
c
c         - Made a small change to the way an SDS is deallocated.
c         - Added a routine to initialize an SDS.  This is useful
c           when deallocating SDS structures.
c
c-----------------------------------------------------------------------
c
c#######################################################################
      module sdslib_ident
c
      character(*), parameter :: cname='SDSLIB'
      character(*), parameter :: cvers='1.04'
      character(*), parameter :: cdate='02/24/2009'
c
      end module
c#######################################################################
      module assign_ptr_1d_interface
      interface
        subroutine assign_ptr_1d (from,to)
        use number_types
        implicit none
        real(r_typ), dimension(:), target :: from
        real(r_typ), dimension(:), pointer :: to
        end
      end interface
      end module
c#######################################################################
      module assign_ptr_3d_interface
      interface
        subroutine assign_ptr_3d (from,to)
        use number_types
        implicit none
        real(r_typ), dimension(:,:,:), target :: from
        real(r_typ), dimension(:,:,:), pointer :: to
        end
      end interface
      end module
c#######################################################################
      subroutine assign_ptr_1d (from,to)
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
      real(r_typ), dimension(:), target :: from
      real(r_typ), dimension(:), pointer :: to
c
c-----------------------------------------------------------------------
c
      to=>from
c
      return
      end
c#######################################################################
      subroutine assign_ptr_3d (from,to)
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
      real(r_typ), dimension(:,:,:), target :: from
      real(r_typ), dimension(:,:,:), pointer :: to
c
c-----------------------------------------------------------------------
c
      to=>from
c
      return
      end
c#######################################################################
      subroutine init_sds_pointer_status (s)
c
c-----------------------------------------------------------------------
c
c ****** Disassociate all the pointers in the SDS in structure S.
c
c-----------------------------------------------------------------------
c
c ****** This is useful when subsequently querying the association
c ****** status of these pointers (e.g., when deallocating storage).
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
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
      nullify (s%f)
c
      nullify (s%scales(1)%f)
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
c
      return
      end
c#######################################################################
      subroutine deallocate_sds (s)
c
c-----------------------------------------------------------------------
c
c ****** Deallocate the memory used by the SDS in structure S.
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
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
      if (associated(s%f)) deallocate (s%f)
c
      if (associated(s%scales(1)%f)) deallocate (s%scales(1)%f)
      if (associated(s%scales(2)%f)) deallocate (s%scales(2)%f)
      if (associated(s%scales(3)%f)) deallocate (s%scales(3)%f)
c
      return
      end
c#######################################################################
      subroutine rdhdf_1d (fname,scale,nx,f,x,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D scientific data set from an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDHDF to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(:), pointer :: f
      real(r_typ), dimension(:), pointer :: x
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdhdf (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 1D data set.
c
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_1D:'
        write (*,*) '### The HDF file does not contain a 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      scale=s%scale
      x=>s%scales(1)%f
c
      allocate (f(nx))
      f=s%f(:,1,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdhdf_2d (fname,scale,nx,ny,f,x,y,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 2D scientific data set from an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDHDF to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdhdf (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 2D data set.
c
      if (s%ndim.ne.2) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_2D:'
        write (*,*) '### The HDF file does not contain a 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
c
      allocate (f(nx,ny))
      f=s%f(:,:,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdhdf_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 3D scientific data set from an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDHDF to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y,z
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,nz,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdhdf (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 3D data set.
c
      if (s%ndim.ne.3) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_3D:'
        write (*,*) '### The HDF file does not contain a 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      nz=s%dims(3)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
      z=>s%scales(3)%f
      f=>s%f
c
      return
      end
c#######################################################################
      subroutine wrhdf_1d (fname,scale,nx,f,x,hdf32,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRHDF to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(nx,1,1) :: f
      real(r_typ), dimension(nx) :: x
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,f,x,hdf32
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=1
      s%dims(1)=nx
      s%dims(2)=1
      s%dims(3)=1
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
      else
        nullify (s%scales(1)%f)
      end if
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrhdf (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_1D:'
        write (*,*) '### Could not write the 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrhdf_2d (fname,scale,nx,ny,f,x,y,hdf32,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 2D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRHDF to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(nx,ny,1) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,f,x,y,hdf32
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=2
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=1
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
      end if
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrhdf (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_2D:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrhdf_3d (fname,scale,nx,ny,nz,f,x,y,z,hdf32,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 3D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRHDF to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,nz,f,x,y,z,hdf32
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=3
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=nz
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
        call assign_ptr_1d (z,s%scales(3)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
        nullify (s%scales(3)%f)
      end if
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrhdf (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_3D:'
        write (*,*) '### Could not write the 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine rdsds (fmt,fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a scientific data set from file FNAME into
c ****** SDS structure S.
c
c ****** Use routine RDTXT or RDHDF, depending on the format
c ****** specified by FMT.
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
      character(*) :: fmt
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fmt,fname
      intent(out) :: s,ierr
c
c-----------------------------------------------------------------------
c
      if (fmt.eq.'text') then
        call rdtxt (fname,s,ierr)
      else if (fmt.eq.'hdf') then
        call rdhdf (fname,s,ierr)
      else
        write (*,*)
        write (*,*) '### ERROR in RDSDS:'
        write (*,*) '### Invalid file format specified.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Format = ',fmt
        ierr=5
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrsds (fmt,fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a scientific data set from SDS structure S to
c ****** file FNAME.
c
c ****** Use routine WRTXT or WRHDF, depending on the format
c ****** specified by FMT.
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
      character(*) :: fmt
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fmt,fname,s
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      if (fmt.eq.'text') then
        call wrtxt (fname,s,ierr)
      else if (fmt.eq.'hdf') then
        call wrhdf (fname,s,ierr)
      else
        write (*,*)
        write (*,*) '### ERROR in WRSDS:'
        write (*,*) '### Invalid file format specified.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Format = ',fmt
        ierr=5
        return
      end if
c
      return
      end
c#######################################################################
      subroutine rdhdf (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D, 2D, or 3D scientific data set from an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine allocates the required memory and returns
c ****** pointers to the data and scale arrays.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    HDF data file name to read from.
c
c ****** Output arguments:
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was read
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions found in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    For a 1D data set, DIMS(2)=DIMS(3)=1.
c                    For a 2D data set, DIMS(3)=1.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    were not found; SCALE=.true. means that scales
c                    were found.
c
c          HDF32   : [logical]
c                    Flag to indicate the precision of the data set
c                    read in.  HDF32=.true. means that the data is
c                    32-bit; HDF32=.false. means that the data is
c                    64-bit.
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is undefined otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c ****** The storage for the arrays pointed to by F, and the
c ****** scales (if present) in structure SCALES, is allocated by
c ****** this routine.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname
      intent(out) :: s,ierr
c
c-----------------------------------------------------------------------
c
      integer, parameter :: mxdim=3
      integer :: iret,nscal,i,n
      integer :: hdf_nt
c
c-----------------------------------------------------------------------
c
      integer, external :: ispdp_dsgdisc,ispdp_dsgdata
c
c-----------------------------------------------------------------------
c
c ****** HDF parameters and routines.
c
      integer, parameter :: DFNT_FLOAT32=5
      integer, parameter :: DFNT_FLOAT64=6
c
      integer, external :: DFSDrestart,dsgdims,dsgnt
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Reset the HDF read routines to read from the beginning
c ****** of the file.  This is required, since, if the previous
c ****** read happened to be from the same file, the next
c ****** read will be from the end of the last data set read,
c ****** and not the beginning of the file.
c
      iret=DFSDrestart()
c
c ****** Get the dimensions of the SDS.
c
c ****** Due to a bug in the HDF library (versions 3.3r4 and 4.0r2),
c ****** set S%NDIM to zero to prevent a core dump if the file is
c ****** missing.
c
      s%ndim=0
c
      s%dims(:)=1
c
      iret=dsgdims(fname,s%ndim,s%dims,mxdim)
c
      if (iret.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF:'
        write (*,*) '### Error while retrieving HDF SDS dimensions.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from DSGDIMS) = ',iret,']'
        ierr=1
        return
      end if
c
c ****** Get the precision of the data in the file.
c
      iret=dsgnt(hdf_nt)
c
      if (iret.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF:'
        write (*,*) '### Error while retrieving HDF-file number type.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from DSGNT) = ',iret,']'
        ierr=2
        return
      end if
c
      if (hdf_nt.eq.DFNT_FLOAT32) then
        s%hdf32=.true.
      else if (hdf_nt.eq.DFNT_FLOAT64) then
        s%hdf32=.false.
      else
        write (*,*)
        write (*,*) '### ERROR in RDHDF:'
        write (*,*) '### Unrecognized number type in HDF file.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Number type = ',hdf_nt
        ierr=3
        return
      end if
c
c ****** Get the scales.
c
      nscal=0
c
      do i=1,s%ndim
        allocate (s%scales(i)%f(s%dims(i)))
        iret=ispdp_dsgdisc(s%hdf32,i,s%dims(i),s%scales(i)%f)
        if (iret.eq.0) nscal=nscal+1
      enddo
c
c ****** Check if all the scales are present.
c
      s%scale=nscal.eq.s%ndim
c
c ****** If scales are not present, deallocate any scale arrays
c ****** that may have been allocated.  Then, allocate dummy
c ****** scales (of length 1) so that the pointers to the scales
c ****** are valid.
c
      if (.not.s%scale) then
        do i=1,s%ndim
          deallocate (s%scales(i)%f)
        enddo
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      else
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      end if
c
c ****** Allocate the memory for the array.
c
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
c
c ****** Read the data array.
c
      n=product(s%dims(1:s%ndim))
      iret=ispdp_dsgdata(s%hdf32,fname,s%ndim,s%dims,n,s%f)
c
      if (iret.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF:'
        write (*,*) '### Error while reading the SDS.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from DSGDATA) = ',iret,']'
        ierr=4
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrhdf (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D, 2D, or 3D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    HDF data file name to write to.
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c ****** Output arguments:
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was written
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    Only DIMS(1 .. NDIM) are referenced.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    are not being supplied; SCALE=.true. means that
c                    scales are being supplied.
c
c          HDF32   : [logical]
c                    Flag to specify the precision of the data to
c                    be written to the file.  Set HDF32=.true. to
c                    write 32-bit data, and HDF32=.false. to write
c                    64-bit data.
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is not referenced
c                    otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      integer :: iret,i,n
      integer :: hdf_nt
c
c-----------------------------------------------------------------------
c
      integer, external :: ispdp_dssdisc,ispdp_dspdata
c
c-----------------------------------------------------------------------
c
c ****** HDF parameters and routines.
c
      integer, parameter :: DFNT_FLOAT32=5
      integer, parameter :: DFNT_FLOAT64=6
c
      integer, external :: DFSDclear,dssdims,dssnt
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Clear the SDS interface.  This erases all previously
c ****** set values.
c
      iret=DFSDclear()
c
c ****** Check the number of dimensions.
c
      if (s%ndim.le.0.or.s%ndim.gt.3) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF:'
        write (*,*) '### Could not write the SDS data.'
        write (*,*) 'Invalid number of dimensions.'
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
      iret=dssdims(s%ndim,s%dims)
c
      if (iret.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF:'
        write (*,*) '### Could not set the SDS dimensions.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from DSSDIMS) = ',iret,']'
        ierr=2
        return
      end if
c
c ****** Set the precision of the data according to the flag S%HDF32.
c
      if (s%hdf32) then
        hdf_nt=DFNT_FLOAT32
      else
        hdf_nt=DFNT_FLOAT64
      end if
c
      iret=dssnt(hdf_nt)
c
      if (iret.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF:'
        write (*,*) '### Error while setting the HDF-file number type.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from DSSNT) = ',iret,']'
        ierr=3
        return
      end if
c
c ****** Set the scales.
c
      if (s%scale) then
        do i=1,s%ndim
          iret=ispdp_dssdisc(s%hdf32,i,s%dims(i),s%scales(i)%f)
          if (iret.ne.0) go to 900
        enddo
      end if
c
c ****** Write the data array.
c
      n=product(s%dims(1:s%ndim))
      iret=ispdp_dspdata(s%hdf32,fname,s%ndim,s%dims,n,s%f)
c
      if (iret.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF:'
        write (*,*) '### Could not write the SDS data.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from DSPDATA) = ',iret,']'
        ierr=4
        return
      end if
c
      return
c
  900 continue
c
      write (*,*)
      write (*,*) '### ERROR in WRHDF:'
      write (*,*) '### Could not set the scales.'
      write (*,*) 'File name: ',trim(fname)
      write (*,*) '[Error return (from DSSDISC) = ',iret,']'
      ierr=5
c
      return
      end
c#######################################################################
      function ispdp_dsgdisc (hdf32,idim,n,f)
c
c-----------------------------------------------------------------------
c
c ****** Wrapper routine to call the HDF routine DSGDISC.
c
c ****** This routine converts the data from the file precision
c ****** (specified by HDF32) to the precision that is being
c ****** used for REALs (specified by R_TYP) before loading F.
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
      logical :: hdf32
      integer :: idim
      integer :: n
      real(r_typ), dimension(n) :: f
      integer :: ispdp_dsgdisc
c
c-----------------------------------------------------------------------
c
      real(KIND_REAL_4), dimension(:), allocatable :: f4
      real(KIND_REAL_8), dimension(:), allocatable :: f8
c
c-----------------------------------------------------------------------
c
      integer :: iret
c
c-----------------------------------------------------------------------
c
c ****** HDF routines.
c
      integer, external :: dsgdisc
c
c-----------------------------------------------------------------------
c
      if (hdf32) then
        if (kind(f).ne.KIND_REAL_4) then
          allocate (f4(n))
          iret=dsgdisc(idim,n,f4)
          f(:)=f4(:)
          deallocate (f4)
        else
          iret=dsgdisc(idim,n,f)
        end if
      else
        if (kind(f).ne.KIND_REAL_8) then
          allocate (f8(n))
          iret=dsgdisc(idim,n,f8)
          f(:)=f8(:)
          deallocate (f8)
        else
          iret=dsgdisc(idim,n,f)
        end if
      end if
c
      ispdp_dsgdisc=iret
c
      return
      end
c#######################################################################
      function ispdp_dsgdata (hdf32,fname,ndim,dims,n,f)
c
c-----------------------------------------------------------------------
c
c ****** Wrapper routine to call the HDF routine DSGDATA.
c
c ****** This routine converts the data from the file precision
c ****** (specified by HDF32) to the precision that is being
c ****** used for REALs (specified by R_TYP) before loading F.
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
      logical :: hdf32
      character(*) :: fname
      integer :: ndim
      integer, dimension(ndim) :: dims
      integer :: n
      real(r_typ), dimension(n) :: f
      integer :: ispdp_dsgdata
c
c-----------------------------------------------------------------------
c
      real(KIND_REAL_4), dimension(:), allocatable :: f4
      real(KIND_REAL_8), dimension(:), allocatable :: f8
c
c-----------------------------------------------------------------------
c
      integer :: iret
c
c-----------------------------------------------------------------------
c
c ****** HDF routines.
c
      integer, external :: dsgdata
c
c-----------------------------------------------------------------------
c
      if (hdf32) then
        if (kind(f).ne.KIND_REAL_4) then
          allocate (f4(n))
          iret=dsgdata(fname,ndim,dims,f4)
          f(:)=f4(:)
          deallocate (f4)
        else
          iret=dsgdata(fname,ndim,dims,f)
        end if
      else
        if (kind(f).ne.KIND_REAL_8) then
          allocate (f8(n))
          iret=dsgdata(fname,ndim,dims,f8)
          f(:)=f8(:)
          deallocate (f8)
        else
          iret=dsgdata(fname,ndim,dims,f)
        end if
      end if
c
      ispdp_dsgdata=iret
c
      return
      end
c#######################################################################
      function ispdp_dssdisc (hdf32,idim,n,f)
c
c-----------------------------------------------------------------------
c
c ****** Wrapper routine to call the HDF routine DSSDISC.
c
c ****** This routine converts the data in F from the precision
c ****** that is being used for REALs (specified by R_TYP) to
c ****** the file precision (specified by HDF32) prior to
c ****** writing F.
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
      logical :: hdf32
      integer :: idim
      integer :: n
      real(r_typ), dimension(n) :: f
      integer :: ispdp_dssdisc
c
c-----------------------------------------------------------------------
c
      real(KIND_REAL_4), dimension(:), allocatable :: f4
      real(KIND_REAL_8), dimension(:), allocatable :: f8
c
c-----------------------------------------------------------------------
c
      integer :: iret
c
c-----------------------------------------------------------------------
c
c ****** HDF routines.
c
      integer, external :: dssdisc
c
c-----------------------------------------------------------------------
c
      if (hdf32) then
        if (kind(f).ne.KIND_REAL_4) then
          allocate (f4(n))
          f4(:)=f(:)
          iret=dssdisc(idim,n,f4)
          deallocate (f4)
        else
          iret=dssdisc(idim,n,f)
        end if
      else
        if (kind(f).ne.KIND_REAL_8) then
          allocate (f8(n))
          f8(:)=f(:)
          iret=dssdisc(idim,n,f8)
          deallocate (f8)
        else
          iret=dssdisc(idim,n,f)
        end if
      end if
c
      ispdp_dssdisc=iret
c
      return
      end
c#######################################################################
      function ispdp_dspdata (hdf32,fname,ndim,dims,n,f)
c
c-----------------------------------------------------------------------
c
c ****** Wrapper routine to call the HDF routine DSPDATA.
c
c ****** This routine converts the data in F from the precision
c ****** that is being used for REALs (specified by R_TYP) to
c ****** the file precision (specified by HDF32) prior to
c ****** writing F.
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
      logical :: hdf32
      character(*) :: fname
      integer :: ndim
      integer, dimension(ndim) :: dims
      integer :: n
      real(r_typ), dimension(n) :: f
      integer :: ispdp_dspdata
c
c-----------------------------------------------------------------------
c
      real(KIND_REAL_4), dimension(:), allocatable :: f4
      real(KIND_REAL_8), dimension(:), allocatable :: f8
c
c-----------------------------------------------------------------------
c
      integer :: iret
c
c-----------------------------------------------------------------------
c
c ****** HDF routines.
c
      integer, external :: dspdata
c
c-----------------------------------------------------------------------
c
      if (hdf32) then
        if (kind(f).ne.KIND_REAL_4) then
          allocate (f4(n))
          f4(:)=f(:)
          iret=dspdata(fname,ndim,dims,f4)
          deallocate (f4)
        else
          iret=dspdata(fname,ndim,dims,f)
        end if
      else
        if (kind(f).ne.KIND_REAL_8) then
          allocate (f8(n))
          f8(:)=f(:)
          iret=dspdata(fname,ndim,dims,f8)
          deallocate (f8)
        else
          iret=dspdata(fname,ndim,dims,f)
        end if
      end if
c
      ispdp_dspdata=iret
c
      return
      end
c#######################################################################
      subroutine rdtxt_1d (fname,scale,nx,f,x,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDTXT to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(:), pointer :: f
      real(r_typ), dimension(:), pointer :: x
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdtxt (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 1D data set.
c
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_1D:'
        write (*,*) '### The test file does not contain a 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      scale=s%scale
      x=>s%scales(1)%f
c
      allocate (f(nx))
      f=s%f(:,1,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdtxt_2d (fname,scale,nx,ny,f,x,y,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 2D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDTXT to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdtxt (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 2D data set.
c
      if (s%ndim.ne.2) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_2D:'
        write (*,*) '### The text file does not contain a 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
c
      allocate (f(nx,ny))
      f=s%f(:,:,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdtxt_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 3D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDTXT to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y,z
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,nz,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdtxt (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 3D data set.
c
      if (s%ndim.ne.3) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_3D:'
        write (*,*) '### The text file does not contain a 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      nz=s%dims(3)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
      z=>s%scales(3)%f
      f=>s%f
c
      return
      end
c#######################################################################
      subroutine wrtxt_1d (fname,scale,nx,f,x,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRTXT to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(nx,1,1) :: f
      real(r_typ), dimension(nx) :: x
      integer :: ierr
      intent(in) :: fname,scale,nx,f,x
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=1
      s%dims(1)=nx
      s%dims(2)=1
      s%dims(3)=1
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
      else
        nullify (s%scales(1)%f)
      end if
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrtxt (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_1D:'
        write (*,*) '### Could not write the 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrtxt_2d (fname,scale,nx,ny,f,x,y,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 2D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRTXT to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(nx,ny,1) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,f,x,y
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=2
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=1
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
      end if
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrtxt (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_2D:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrtxt_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 3D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRTXT to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,nz,f,x,y,z
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=3
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=nz
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
        call assign_ptr_1d (z,s%scales(3)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
        nullify (s%scales(3)%f)
      end if
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrtxt (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_3D:'
        write (*,*) '### Could not write the 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine rdtxt (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D, 2D, or 3D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine allocates the required memory and returns
c ****** pointers to the data and scale arrays.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    Text data file name to read from.
c
c ****** Output arguments:
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was read
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions found in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    For a 1D data set, DIMS(2)=DIMS(3)=1.
c                    For a 2D data set, DIMS(3)=1.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    were not found; SCALE=.true. means that scales
c                    were found.
c
c          HDF32   : [logical]
c                    This flag is is not relevant to text files;
c                    it is used for HDF data.
c                    It is arbitrarily set to HDF32=.false..
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is undefined otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c ****** The storage for the arrays pointed to by F, and the
c ****** scales (if present) in structure SCALES, is allocated by
c ****** this routine.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname
      intent(out) :: s,ierr
c
c-----------------------------------------------------------------------
c
      integer, parameter :: mxdim=3
      integer :: ifscale,i,n
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Open the file for reading.
c
      call ffopen (1,fname,'r',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT:'
        write (*,*) '### Could not open the text file.'
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Get the number of dimensions.
c
      call rdint (1,1,s%ndim,ierr)
      if (ierr.ne.0) go to 910
c
      if (s%ndim.lt.1.or.s%ndim.gt.mxdim) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT:'
        write (*,*) '### Invalid number of dimensions in file.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'Maximum number of dimensions = ',mxdim
        ierr=2
        return
      end if
c
c ****** Read the dimensions.
c
      s%dims(:)=1
c
      do i=1,s%ndim
        call rdint (1,1,s%dims(i),ierr)
        if (ierr.ne.0) go to 910
        if (s%dims(i).le.0) go to 920
      enddo
c
c ****** Check if the scales are present.
c
      call rdint (1,1,ifscale,ierr)
      if (ierr.ne.0) go to 910
c
      s%scale=ifscale.ne.0
c
c ****** Allocate memory and read the scales (if present).
c
c ****** If scales are not present,  allocate dummy scales
c ****** (of length 1) so that the pointers to the scales
c ****** are valid.
c
      if (s%scale) then
        do i=1,s%ndim
          allocate (s%scales(i)%f(s%dims(i)))
          call rdfp (1,s%dims(i),s%scales(i)%f,ierr)
          if (ierr.ne.0) go to 910
        enddo
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      else
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      end if
c
c ****** Allocate memory for the array.
c
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
c
c ****** Read the data array.
c
      n=product(s%dims(1:s%ndim))
      call rdfp (1,n,s%f,ierr)
      if (ierr.ne.0) go to 910
c
      s%hdf32=.false.
c
      close (1)
c
      return
c
  910 continue
c
      write (*,*)
      write (*,*) '### ERROR in RDTXT:'
      write (*,*) '### Error while reading text data.'
      write (*,*) 'File name: ',trim(fname)
      ierr=3
      return
c
  920 continue
c
      write (*,*)
      write (*,*) '### ERROR in RDTXT:'
      write (*,*) '### Invalid value for dimension.'
      write (*,*) 'Dimension number = ',i
      write (*,*) 'Dimension value = ',s%dims(i)
      ierr=4
c
      return
      end
c#######################################################################
      subroutine rdint (iun,n,i,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read N words of INTEGER data into array I from unit IUN
c ****** using a free format text read.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: iun
      integer :: n
      integer, dimension(n) :: i
      integer :: ierr
      intent(in) :: iun,n
      intent(out) :: i,ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      read (iun,*,err=100,end=100) i
c
      return
c
  100 continue
c
c ****** Error in reading the data.
c
      ierr=1
c
      return
      end
c#######################################################################
      subroutine rdfp (iun,n,f,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read N words of REAL data into array F from unit IUN
c ****** using a free format text read.
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
      integer :: iun
      integer :: n
      real(r_typ), dimension(n) :: f
      integer :: ierr
      intent(in) :: iun,n
      intent(out) :: f,ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      read (iun,*,err=100,end=100) f
c
      return
c
  100 continue
c
c ****** Error in reading the data.
c
      ierr=1
c
      return
      end
c#######################################################################
      subroutine wrtxt (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D, 2D, or 3D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    Text data file name to write to.
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c ****** Output arguments:
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was written
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    Only DIMS(1 .. NDIM) are referenced.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    are not being supplied; SCALE=.true. means that
c                    scales are being supplied.
c
c          HDF32   : [logical]
c                    Flag that indicates the precision of the data.
c                    This flag is used to determine the format for data
c                    written to the text file.  When HDF32=.TRUE., the
c                    data is assumed to originate from a 32-bit HDF data
c                    file, and is written with 7 digits to the text file.
c                    Otherwise, the data is assumed to originate from a
c                    64-bit HDF data file, and is written with 14 digits
c                    to the text file.
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is not referenced
c                    otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declarations for temporary variables.
c
      integer :: i,n
      character(32) :: fmt
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Open the file for writing.
c
      call ffopen (1,fname,'rw',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT:'
        write (*,*) '### Could not open the text file for writing.'
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Check the number of dimensions.
c
      if (s%ndim.le.0.or.s%ndim.gt.3) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT:'
        write (*,*) '### Could not write the SDS data.'
        write (*,*) 'Invalid number of dimensions.'
        write (*,*) 'NDIM = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Construct the format string for writing floating point
c ****** numbers to the output file.
c
      if (s%hdf32) then
        fmt='(5(1x,1pe13.6))'
      else
        fmt='(3(1x,1pe21.14))'
      end if
c
c ****** Write the number of dimensions.
c
      call wrint (1,1,s%ndim,ierr)
      if (ierr.ne.0) go to 900
c
c ****** Write the dimensions.
c
      do i=1,s%ndim
        call wrint (1,1,s%dims(i),ierr)
        if (ierr.ne.0) go to 900
      enddo
c
c ****** Write the scales.
c
      if (s%scale) then
        call wrint (1,1,1,ierr)
        if (ierr.ne.0) go to 900
        do i=1,s%ndim
          call wrfp (1,s%dims(i),s%scales(i)%f,fmt,ierr)
          if (ierr.ne.0) go to 900
        enddo
      else
        call wrint (1,1,0,ierr)
        if (ierr.ne.0) go to 900
      end if
c
c ****** Write the array.
c
      n=product(s%dims(1:s%ndim))
      call wrfp (1,n,s%f,fmt,ierr)
      if (ierr.ne.0) go to 900
c
      close (1)
c
      return
c
  900 continue
c
      write (*,*)
      write (*,*) '### ERROR in WRTXT:'
      write (*,*) '### Error in writing data to the text file.'
      write (*,*) 'File name: ',trim(fname)
      ierr=2
c
      return
      end
c#######################################################################
      subroutine wrint (iun,n,i,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write N words of INTEGER data from array I to the file
c ****** connected to unit IUN using a free format text write.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: iun
      integer :: n
      integer, dimension(n) :: i
      integer :: ierr
      intent(in) :: iun,n,i
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      write (iun,*,err=100) i
c
      return
c
  100 continue
c
c ****** Error in writing the data.
c
      ierr=1
c
      return
      end
c#######################################################################
      subroutine wrfp (iun,n,f,fmt,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write N words of REAL data from array F to the file
c ****** connected to unit IUN.
c
c ****** FMT specifies the format string to use.
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
      integer :: iun
      integer :: n
      real(r_typ), dimension(n) :: f
      character(*) :: fmt
      integer :: ierr
      intent(in) :: iun,n,f,fmt
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      write (iun,fmt=fmt,err=100) f
c
      return
c
  100 continue
c
c ****** Error in writing the data.
c
      ierr=1
c
      return
      end
