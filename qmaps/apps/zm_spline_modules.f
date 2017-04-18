c#######################################################################
      module spline_def
c
c-----------------------------------------------------------------------
c ****** Definition of cubic spline data structures.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ***** 1D spline structure.
c
      type :: spl1d
        integer :: nx
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: f
        real(r_typ), dimension(:), pointer :: fxx
      end type
c
c ***** 2D spline structure.
c
      type :: spl2d
        integer :: nx
        integer :: ny
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: y
        real(r_typ), dimension(:,:), pointer :: f
        real(r_typ), dimension(:,:), pointer :: fxx
        real(r_typ), dimension(:,:), pointer :: fyy
        real(r_typ), dimension(:,:), pointer :: fxxyy
      end type
c
c ***** 3D spline structure.
c
      type :: spl3d
        integer :: nx
        integer :: ny
        integer :: nz
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: y
        real(r_typ), dimension(:), pointer :: z
        real(r_typ), dimension(:,:,:), pointer :: f
        real(r_typ), dimension(:,:,:), pointer :: fxx
        real(r_typ), dimension(:,:,:), pointer :: fyy
        real(r_typ), dimension(:,:,:), pointer :: fzz
        real(r_typ), dimension(:,:,:), pointer :: fxxyy
        real(r_typ), dimension(:,:,:), pointer :: fxxzz
        real(r_typ), dimension(:,:,:), pointer :: fyyzz
        real(r_typ), dimension(:,:,:), pointer :: fxxyyzz
      end type
c
      end module
c#######################################################################
      module invint_def
c
c-----------------------------------------------------------------------
c ****** Definition of an inverse interpolation table data structure.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      type :: itab
        integer :: n
        real(r_typ), dimension(:), pointer :: f
        real(r_typ) :: d
      end type
c
      end module
c#######################################################################
      module locate_interval_interface
      interface
        function locate_interval (n,x,xv,tab,ierr)
        use number_types
        use invint_def
        implicit none
        integer :: n
        real(r_typ), dimension(n) :: x
        real(r_typ) :: xv
        type(itab), optional :: tab
        integer, optional :: ierr
        integer :: locate_interval
        end
      end interface
      end module
c#######################################################################
      module splint_interface
      interface
        function splint (n,x,f,fpp,xv,tab)
        use number_types
        use invint_def
        use locate_interval_interface
        implicit none
        integer :: n
        real(r_typ), dimension(n) :: x,f,fpp
        real(r_typ) :: xv
        type(itab), optional :: tab
        real(r_typ) :: splint
        end
      end interface
      end module
c#######################################################################
      module evaluate_spline_1d_interface
      interface
        function evaluate_spline_1d (s,x,tab)
        use number_types
        use spline_def
        use invint_def
        use locate_interval_interface
        type(spl1d) :: s
        real(r_typ) :: x
        type(itab), optional :: tab
        real(r_typ) :: evaluate_spline_1d
        end
      end interface
      end module
c#######################################################################
      module evaluate_spline_2d_interface
      interface
        function evaluate_spline_2d (s,x,y,tabx,taby)
        use number_types
        use spline_def
        use invint_def
        use locate_interval_interface
        type(spl2d) :: s
        real(r_typ) :: x,y
        type(itab), optional :: tabx,taby
        real(r_typ) :: evaluate_spline_2d
        end
      end interface
      end module
c#######################################################################
      module evaluate_spline_3d_interface
      interface
        function evaluate_spline_3d (s,x,y,z,tabx,taby,tabz)
        use number_types
        use spline_def
        use invint_def
        use locate_interval_interface
        type(spl3d) :: s
        real(r_typ) :: x,y,z
        type(itab), optional :: tabx,taby,tabz
        real(r_typ) :: evaluate_spline_3d
        end
      end interface
      end module
