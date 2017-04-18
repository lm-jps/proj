c
c-----------------------------------------------------------------------
c
c ****** Source to build the spline library.
c ****** These routines are used by Zoran Mikic's tools.
c
c-----------------------------------------------------------------------
c
c        07/16/2006, ZM, Version 1.00:
c
c         - Original version of the spline library.
c           This library was put together to facilitate the
c           development of ZM's tools.
c           It includes routines to perform cubic spline
c           interpolation in 1D, 2D, and 3D.
c
c        09/29/2006, ZM, Version 1.01:
c
c         - Changed routine LOCATE_INTERVAL so that it handles
c           tables with only one point properly.
c
c        02/16/2009, ZM, Version 1.02:
c
c         - Added the ability to specify periodic boundary conditions
c           for splines.
c
c-----------------------------------------------------------------------
c
c#######################################################################
      module splinelib_ident
c
      character(*), parameter :: cname='SPLINELIB'
      character(*), parameter :: cvers='1.02'
      character(*), parameter :: cdate='02/16/2009'
c
      end module
c#######################################################################
      subroutine compute_spline_1d (nx,x,f,s)
c
c-----------------------------------------------------------------------
c
c ****** Find the cubic spline coefficients for the 1D function
c ****** defined by array F(NX) with scale X(NX).
c
c ****** The spline coefficients are returned in structure S. 
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(nx) :: f
      type(spl1d) :: s
c
c-----------------------------------------------------------------------
c
c ****** Allocate storage for the spline coefficients.
c
      s%nx=nx
c
      allocate (s%x(nx))
      allocate (s%f(nx))
      allocate (s%fxx(nx))
c
c ****** Evaluate the spline coefficients.
c
      s%x=x
      s%f=f
c
      call ezspline (nx,x,s%f,s%fxx)
c
      return
      end
c#######################################################################
      subroutine compute_spline_2d (nx,ny,x,y,f,s)
c
c-----------------------------------------------------------------------
c
c ****** Find the cubic spline coefficients for the 2D function
c ****** defined by array F(NX,NY), with scales X(NX) and Y(NY).
c
c ****** The spline coefficients are returned in structure S. 
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
      type(spl2d) :: s
c
c-----------------------------------------------------------------------
c
      integer :: i,j
      real(r_typ), dimension(nx) :: gx,gppx
      real(r_typ), dimension(ny) :: gy,gppy
c
c-----------------------------------------------------------------------
c
c ****** Allocate storage for the spline coefficients.
c
      s%nx=nx
      s%ny=ny
c
      allocate (s%x(nx))
      allocate (s%y(ny))
      allocate (s%f(nx,ny))
      allocate (s%fxx(nx,ny))
      allocate (s%fyy(nx,ny))
      allocate (s%fxxyy(nx,ny))
c
c ****** Evaluate the spline coefficients.
c
      s%x=x
      s%y=y
      s%f=f
c
      do j=1,ny
        gx(:)=s%f(:,j)
        call ezspline (nx,x,gx,gppx)
        s%fxx(:,j)=gppx(:)
      enddo
c
      do i=1,nx
        gy(:)=s%f(i,:)
        call ezspline (ny,y,gy,gppy)
        s%fyy(i,:)=gppy(:)
      enddo
c
      do i=1,nx
        gy(:)=s%fxx(i,:)
        call ezspline (ny,y,gy,gppy)
        s%fxxyy(i,:)=gppy(:)
      enddo
c
      return
      end
c#######################################################################
      subroutine compute_spline_3d (nx,ny,nz,x,y,z,f,s)
c
c-----------------------------------------------------------------------
c
c ****** Find the cubic spline coefficients for the 3D function
c ****** defined by array F(NX,NY,NZ), with scales X(NX), Y(NY),
c ****** and Z(NZ).
c
c ****** The spline coefficients are returned in structure S. 
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
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
      real(r_typ), dimension(nx,ny,nz) :: f
      type(spl3d) :: s
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
      real(r_typ), dimension(nx) :: gx,gppx
      real(r_typ), dimension(ny) :: gy,gppy
      real(r_typ), dimension(nz) :: gz,gppz
c
c-----------------------------------------------------------------------
c
c ****** Allocate storage for the spline coefficients.
c
      s%nx=nx
      s%ny=ny
      s%nz=nz
c
      allocate (s%x(nx))
      allocate (s%y(ny))
      allocate (s%z(nz))
      allocate (s%f(nx,ny,nz))
      allocate (s%fxx(nx,ny,nz))
      allocate (s%fyy(nx,ny,nz))
      allocate (s%fzz(nx,ny,nz))
      allocate (s%fxxyy(nx,ny,nz))
      allocate (s%fxxzz(nx,ny,nz))
      allocate (s%fyyzz(nx,ny,nz))
      allocate (s%fxxyyzz(nx,ny,nz))
c
c ****** Evaluate the spline coefficients.
c
      s%x=x
      s%y=y
      s%z=z
      s%f=f
c
      do k=1,nz
        do j=1,ny
          gx(:)=s%f(:,j,k)
          call ezspline (nx,x,gx,gppx)
          s%fxx(:,j,k)=gppx(:)
        enddo
      enddo
c
      do k=1,nz
        do i=1,nx
          gy(:)=s%f(i,:,k)
          call ezspline (ny,y,gy,gppy)
          s%fyy(i,:,k)=gppy(:)
        enddo
      enddo
c
      do j=1,ny
        do i=1,nx
          gz(:)=s%f(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fzz(i,j,:)=gppz(:)
        enddo
      enddo
c
      do k=1,nz
        do i=1,nx
          gy(:)=s%fxx(i,:,k)
          call ezspline (ny,y,gy,gppy)
          s%fxxyy(i,:,k)=gppy(:)
        enddo
      enddo
c
      do j=1,ny
        do i=1,nx
          gz(:)=s%fxx(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fxxzz(i,j,:)=gppz(:)
          gz(:)=s%fyy(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fyyzz(i,j,:)=gppz(:)
          gz(:)=s%fxxyy(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fxxyyzz(i,j,:)=gppz(:)
        enddo
      enddo
c
      return
      end
c#######################################################################
      function evaluate_spline_1d (s,x,tab)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the 1D spline in structure S at the
c ****** point X.
c
c ****** The optional argument TAB is an inverse interpolation
c ****** table that can be used to speed up the search for the
c ****** interval that contains X.
c
c ****** The cubic spline coefficients for the spline S can
c ****** be obtained using routine COMPUTE_SPLINE_1D.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(spl1d) :: s
      real(r_typ) :: x
      type(itab), optional :: tab
      real(r_typ) :: evaluate_spline_1d
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
c
c-----------------------------------------------------------------------
c
c ****** Find the index of the cell enclosing point X.
c
      if (present(tab)) then
        i=locate_interval(s%nx,s%x,x,tab)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
c
c ****** Interpolate in x.
c
      evaluate_spline_1d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          s%f(i  ),
     &                          s%f(i+1),
     &                          s%fxx(i  ),
     &                          s%fxx(i+1),
     &                          x)
c
      return
      end
c#######################################################################
      function evaluate_spline_2d (s,x,y,tabx,taby)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the 2D spline in structure S at the
c ****** point (X,Y).
c
c ****** The optional arguments TABX and TABY are inverse
c ****** interpolation tables that can be used to speed up the
c ****** search for the cell that contains (X,Y).
c
c ****** The cubic spline coefficients for the spline S can
c ****** be obtained using routine COMPUTE_SPLINE_2D.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(spl2d) :: s
      real(r_typ) :: x,y
      type(itab), optional :: tabx,taby
      real(r_typ) :: evaluate_spline_2d
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ii
      real(r_typ), dimension(0:1) :: f1,fxx1
      real(r_typ) :: a,b,c,d
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
      real(r_typ), external :: speval_abcd
c
c-----------------------------------------------------------------------
c
c ****** Flag to use the more efficient version of this routine.
c
      logical, parameter :: use_faster=.true.
c
c-----------------------------------------------------------------------
c
c ****** Find the indices of the cell enclosing (X,Y).
c
      if (present(tabx)) then
        i=locate_interval(s%nx,s%x,x,tabx)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
c
      if (present(taby)) then
        j=locate_interval(s%ny,s%y,y,taby)
      else
        j=locate_interval(s%ny,s%y,y)
      end if
c
      if (use_faster) go to 100
c
c ****** Leff efficient version.  This version
c ****** is slower but slightly easier to understand.
c
c ****** Interpolate in y.
c
      do ii=0,1
        f1(ii)=speval(s%y(j  ),
     &                s%y(j+1),
     &                s%f(i+ii,j  ),
     &                s%f(i+ii,j+1),
     &                s%fyy(i+ii,j  ),
     &                s%fyy(i+ii,j+1),
     &                y)
        fxx1(ii)=speval(s%y(j  ),
     &                  s%y(j+1),
     &                  s%fxx(i+ii,j  ),
     &                  s%fxx(i+ii,j+1),
     &                  s%fxxyy(i+ii,j  ),
     &                  s%fxxyy(i+ii,j+1),
     &                  y)
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_2d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
c
  100 continue
c
c ****** More efficient version.  This version
c ****** is slighty faster than the one above.
c
c ****** Interpolate in y.
c
      call speval_get_abcd (s%y(j),s%y(j+1),y,a,b,c,d)
c
      do ii=0,1
        f1(ii)=speval_abcd(a,b,c,d,
     &                     s%f(i+ii,j  ),
     &                     s%f(i+ii,j+1),
     &                     s%fyy(i+ii,j  ),
     &                     s%fyy(i+ii,j+1))
        fxx1(ii)=speval_abcd(a,b,c,d,
     &                       s%fxx(i+ii,j  ),
     &                       s%fxx(i+ii,j+1),
     &                       s%fxxyy(i+ii,j  ),
     &                       s%fxxyy(i+ii,j+1))
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_2d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
      end
c#######################################################################
      function evaluate_spline_3d (s,x,y,z,tabx,taby,tabz)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the 3D spline in structure S at the
c ****** point (X,Y,Z).
c
c ****** The optional arguments TABX, TABY, and TABZ are inverse
c ****** interpolation tables that can be used to speed up the
c ****** search for the cell that contains (X,Y,Z).
c
c ****** The cubic spline coefficients for the spline S can
c ****** be obtained using routine COMPUTE_SPLINE_3D.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(spl3d) :: s
      real(r_typ) :: x,y,z
      type(itab), optional :: tabx,taby,tabz
      real(r_typ) :: evaluate_spline_3d
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,ii,jj
      real(r_typ), dimension(0:1) :: f1,fxx1
      real(r_typ), dimension(0:1,0:1) :: f2,fxx2,fyy2,fxxyy2
      real(r_typ) :: a,b,c,d
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
      real(r_typ), external :: speval_abcd
c
c-----------------------------------------------------------------------
c
c ****** Flag to use the more efficient version of this routine.
c
      logical, parameter :: use_faster=.true.
c
c-----------------------------------------------------------------------
c
c ****** Find the indices of the cell enclosing (X,Y,Z).
c
      if (present(tabx)) then
        i=locate_interval(s%nx,s%x,x,tabx)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
c
      if (present(taby)) then
        j=locate_interval(s%ny,s%y,y,taby)
      else
        j=locate_interval(s%ny,s%y,y)
      end if
c
      if (present(tabz)) then
        k=locate_interval(s%nz,s%z,z,tabz)
      else
        k=locate_interval(s%nz,s%z,z)
      end if
c
      if (use_faster) go to 100
c
c ****** Leff efficient version.  This version
c ****** is slower but slightly easier to understand.
c
c ****** Interpolate in z.
c
      do jj=0,1
        do ii=0,1
          f2(ii,jj)=speval(s%z(k  ),
     &                     s%z(k+1),
     &                     s%f(i+ii,j+jj,k  ),
     &                     s%f(i+ii,j+jj,k+1),
     &                     s%fzz(i+ii,j+jj,k  ),
     &                     s%fzz(i+ii,j+jj,k+1),
     &                     z)
          fxx2(ii,jj)=speval(s%z(k  ),
     &                       s%z(k+1),
     &                       s%fxx(i+ii,j+jj,k  ),
     &                       s%fxx(i+ii,j+jj,k+1),
     &                       s%fxxzz(i+ii,j+jj,k  ),
     &                       s%fxxzz(i+ii,j+jj,k+1),
     &                       z)
          fyy2(ii,jj)=speval(s%z(k  ),
     &                       s%z(k+1),
     &                       s%fyy(i+ii,j+jj,k  ),
     &                       s%fyy(i+ii,j+jj,k+1),
     &                       s%fyyzz(i+ii,j+jj,k  ),
     &                       s%fyyzz(i+ii,j+jj,k+1),
     &                       z)
          fxxyy2(ii,jj)=speval(s%z(k  ),
     &                         s%z(k+1),
     &                         s%fxxyy(i+ii,j+jj,k  ),
     &                         s%fxxyy(i+ii,j+jj,k+1),
     &                         s%fxxyyzz(i+ii,j+jj,k  ),
     &                         s%fxxyyzz(i+ii,j+jj,k+1),
     &                         z)
        enddo
      enddo
c
c ****** Interpolate in y.
c
      do ii=0,1
        f1(ii)=speval(s%y(j  ),
     &                s%y(j+1),
     &                f2(ii,0),
     &                f2(ii,1),
     &                fyy2(ii,0),
     &                fyy2(ii,1),
     &                y)
        fxx1(ii)=speval(s%y(j  ),
     &                  s%y(j+1),
     &                  fxx2(ii,0),
     &                  fxx2(ii,1),
     &                  fxxyy2(ii,0),
     &                  fxxyy2(ii,1),
     &                  y)
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_3d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
c
  100 continue
c
c ****** More efficient version.  This version
c ****** is about 25% faster than the one above.
c
c ****** Interpolate in z.
c
      call speval_get_abcd (s%z(k),s%z(k+1),z,a,b,c,d)
c
      do jj=0,1
        do ii=0,1
          f2(ii,jj)=speval_abcd(a,b,c,d,
     &                          s%f(i+ii,j+jj,k  ),
     &                          s%f(i+ii,j+jj,k+1),
     &                          s%fzz(i+ii,j+jj,k  ),
     &                          s%fzz(i+ii,j+jj,k+1))
          fxx2(ii,jj)=speval_abcd(a,b,c,d,
     &                            s%fxx(i+ii,j+jj,k  ),
     &                            s%fxx(i+ii,j+jj,k+1),
     &                            s%fxxzz(i+ii,j+jj,k  ),
     &                            s%fxxzz(i+ii,j+jj,k+1))
          fyy2(ii,jj)=speval_abcd(a,b,c,d,
     &                            s%fyy(i+ii,j+jj,k  ),
     &                            s%fyy(i+ii,j+jj,k+1),
     &                            s%fyyzz(i+ii,j+jj,k  ),
     &                            s%fyyzz(i+ii,j+jj,k+1))
          fxxyy2(ii,jj)=speval_abcd(a,b,c,d,
     &                              s%fxxyy(i+ii,j+jj,k  ),
     &                              s%fxxyy(i+ii,j+jj,k+1),
     &                              s%fxxyyzz(i+ii,j+jj,k  ),
     &                              s%fxxyyzz(i+ii,j+jj,k+1))
        enddo
      enddo
c
c ****** Interpolate in y.
c
      call speval_get_abcd (s%y(j),s%y(j+1),y,a,b,c,d)
c
      do ii=0,1
        f1(ii)=speval_abcd(a,b,c,d,
     &                     f2(ii,0),
     &                     f2(ii,1),
     &                     fyy2(ii,0),
     &                     fyy2(ii,1))
        fxx1(ii)=speval_abcd(a,b,c,d,
     &                       fxx2(ii,0),
     &                       fxx2(ii,1),
     &                       fxxyy2(ii,0),
     &                       fxxyy2(ii,1))
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_3d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
      end
c#######################################################################
      subroutine spline (n,x,f,ibc0,c0,ibc1,c1,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients.
c
c-----------------------------------------------------------------------
c
c ****** Get the coefficients of a cubic spline interpolant to
c ****** the function values F(i) defined at the points X(i),
c ****** i=1,...,N.
c
c ****** The computed coefficients (which are actually the second
c ****** derivatives of F at the mesh points) are returned in the
c ****** array FPP.
c
c ****** Use routine SPLINT to evaluate the spline at a
c ****** particular position.
c
c-----------------------------------------------------------------------
c
c ****** The boundary conditions at the two ends are specified
c ****** by IBC0 and IBC1, and the coefficient arrays C0 and C1.
c
c ****** These are defined at x=X(1) according to the value
c ****** of IBC0 and C0.  (The conditions at x=X(N) are 
c ****** specified similarly by IBC1 and C1.) 
c
c        IBC0 = 1:  Set the second derivative at one end of the
c                   cell to equal the second derivative at the
c                   other end of the cell. [f''(1)=f''(2)]
c
c        IBC0 = 2:  Set the second derivative to zero.
c                   This corresponds to a "natural spline".
c                   [f''(1)=0.]
c
c        IBC0 = 3:  Set the first derivative to zero.
c                   [f'(1)=0.]
c
c        IBC0 = 4:  Set the second derivative to C0(1).
c                   [f''(1)=C0(1)]
c
c        IBC0 = 5:  Set the first derivative to C0(1).
c                   [f'(1)=C0(1)]
c
c        IBC0 = 6:  Set a linear combination of the first and
c                   second derivatives, according to C0(1),
c                   C0(2), and C0(3).
c                   [C0(2)*f'(1)+C0(3)*f''(1)=C0(1)]
c
c        IBC0 = 7:  Set a linear combination of the second
c                   derivatives at the left and right ends
c                   of the first cell.
c                   [C0(1)*f''(1)+C0(2)*f''(2)=C0(3)]
c
c        IBC0 = 8:  Set a linear combination of the first
c                   derivatives at the left and right ends
c                   of the first cell.
c                   [C0(1)*f'(1)+C0(2)*f'(2)=C0(3)]
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
      real(r_typ), dimension(n) :: x,f,fpp
      integer :: ibc0,ibc1
      real(r_typ), dimension(3) :: c0,c1
c
      intent(in) :: n,x,f,ibc0,c0,ibc1,c1
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ierr
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(n) :: a,b,c
c
c-----------------------------------------------------------------------
c
c ****** Check that there are at least 3 points.
c
      if (n.lt.3) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 3 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
c
c ****** Check that the mesh is monotonic.
c
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the coefficients for the tridiagonal solve at the
c ****** internal points.
c-----------------------------------------------------------------------
c
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the boundary condition at X(1).
c-----------------------------------------------------------------------
c
c ****** Set the unused value to zero.
c
      c(1)=0. 
c
      select case (ibc0)
      case (1)
c
c ****** Second derivatives equal at the left and right ends
c ****** of the first cell.
c
        a(1)=one
        b(1)=-one
        fpp(1)=0.
c
      case (2)
c
c ****** Second derivative is zero ("natural splines").
c
        a(1)=one
        b(1)=0.
        fpp(1)=0.
c
      case (3)
c
c ****** First derivative is zero.
c
        dx=x(2)-x(1)
        a(1)=dx*third
        b(1)=dx*sixth
        fpp(1)=(f(2)-f(1))/dx
c
      case (4)
c
c ****** Second derivative is specified.
c
        a(1)=one
        b(1)=0.
        fpp(1)=c0(1)
c
      case (5)
c
c ****** First derivative is specified.
c
        dx=x(2)-x(1)
        a(1)=dx*third
        b(1)=dx*sixth
        fpp(1)=(f(2)-f(1))/dx-c0(1)
c
      case (6)
c
c ****** A combination of the first and second derivatives
c ****** is specified.
c
        if (c0(2).eq.0..and.c0(3).eq.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### Invalid boundary condition specified'//
     &                ' at X(1).'
          write (*,*) '### Boundary condition type IBC0 = 6.'
          write (*,*) '### It is illegal for both C0(2) and C0(3)'//
     &                ' to be zero.'
          call exit (1)
        end if
c
        if (c0(2).ne.0.) then
          dx=x(2)-x(1)
          a(1)=dx*third-c0(3)/c0(2)
          b(1)=dx*sixth
          fpp(1)=(f(2)-f(1))/dx-c0(1)/c0(2)
        else
          a(1)=one
          b(1)=0.
          fpp(1)=c0(1)/c0(3)
        end if
c
      case (7)
c
c ****** A linear combination of the second derivatives at the left
c ****** and right ends of the first cell is specified.
c
        a(1)=c0(1)
        b(1)=c0(2)
        fpp(1)=c0(3)
c
      case (8)
c
c ****** A linear combination of the first derivatives at the left
c ****** and right ends of the first cell is specified.
c
        dx=x(2)-x(1)
        a(1)=dx*(c0(1)/3-c0(2)/6)
        b(1)=dx*(c0(1)/6-c0(2)/3)
        fpp(1)=(c0(1)+c0(2))*(f(2)-f(1))/dx-c0(3)
c
      case default
c
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid boundary condition specified at X(1).'
        write (*,*) '### IBC0 is invalid.'
        write (*,*) 'IBC0 = ',ibc0
        call exit (1)
c
      end select
c
c-----------------------------------------------------------------------
c ****** Set the boundary condition at X(N).
c-----------------------------------------------------------------------
c
c ****** Set the unused value to zero.
c
      b(n)=0. 
c
      select case (ibc1)
      case (1)
c
c ****** Second derivatives equal at the left and right ends
c ****** of the last cell.
c
        a(n)=one
        c(n)=-one
        fpp(n)=0.
c
      case (2)
c
c ****** Second derivative is zero ("natural splines").
c
        a(n)=one
        c(n)=0.
        fpp(n)=0.
c
      case (3)
c
c ****** First derivative is zero.
c
        dx=x(n)-x(n-1)
        a(n)=dx*third
        c(n)=dx*sixth
        fpp(n)=-(f(n)-f(n-1))/dx
c
      case (4)
c
c ****** Second derivative is specified.
c
        a(n)=one
        c(n)=0.
        fpp(n)=c1(1)
c
      case (5)
c
c ****** First derivative is specified.
c
        dx=x(n)-x(n-1)
        a(n)=dx*third
        c(n)=dx*sixth
        fpp(n)=c1(1)-(f(n)-f(n-1))/dx
c
      case (6)
c
c ****** A combination of the first and second derivatives
c ****** is specified.
c
        if (c1(2).eq.0..and.c1(3).eq.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### Invalid boundary condition specified'//
     &                ' at X(N).'
          write (*,*) '### Boundary condition type IBC1 = 6.'
          write (*,*) '### It is illegal for both C1(2) and C1(3)'//
     &                ' to be zero.'
          call exit (1)
        end if
c
        if (c1(2).ne.0.) then
          dx=x(n)-x(n-1)
          a(n)=dx*third+c1(3)/c1(2)
          c(n)=dx*sixth
          fpp(n)=c1(1)/c1(2)-(f(n)-f(n-1))/dx   
        else
          a(n)=one
          c(n)=0.
          fpp(n)=c1(1)/c1(3)
        end if
c
      case (7)
c
c ****** A linear combination of the second derivatives at the left
c ****** and right ends of the last cell is specified.
c
        a(n)=c1(1)
        c(n)=c1(2)
        fpp(n)=c1(3)
c
      case (8)
c
c ****** A linear combination of the first derivatives at the left
c ****** and right ends of the last cell is specified.
c
        dx=x(n)-x(n-1)
        a(n)=dx*(c1(1)/3-c1(2)/6)
        c(n)=dx*(c1(1)/6-c1(2)/3)
        fpp(n)=c1(3)-(c1(1)+c1(2))*(f(n)-f(n-1))/dx
c
      case default
c
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid boundary condition specified at X(N).'
        write (*,*) '### IBC1 is invalid.'
        write (*,*) 'IBC1 = ',ibc1
        call exit (1)
c
      end select
c
c-----------------------------------------------------------------------
c ****** Solve the tridiagonal system for the second derivative.
c-----------------------------------------------------------------------
c
c ****** The second derivative is returned in FPP.
c
      call trid (n,c,a,b,fpp,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### The tridiagonal matrix relating the'//
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine spline_periodic_type1 (n,x,f,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients for a periodic function.
c
c-----------------------------------------------------------------------
c
c ****** Get the coefficients of a cubic spline interpolant to
c ****** the function values F(i) defined at the points X(i),
c ****** i=1,...,N.
c
c ****** The computed coefficients (which are actually the second
c ****** derivatives of F at the mesh points) are returned in the
c ****** array FPP.
c
c ****** Use routine SPLINT to evaluate the spline at a
c ****** particular position.
c
c-----------------------------------------------------------------------
c
c ****** This routine assumes that the data in F is periodic, such
c ****** that F(N) = F(1), so that the first and last point are
c ****** repeated and represent the same location.  Thus the
c ****** periodicity length is X(N)-X(1).
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
      real(r_typ), dimension(n) :: x,f,fpp
c
      intent(in) :: n,x,f
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ierr,m
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(n-1) :: a,b,c
c
c-----------------------------------------------------------------------
c
c ****** Check that there are at least 3 points.
c
      if (n.lt.3) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 3 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
c
c ****** Check that the mesh is monotonic.
c
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the coefficients for the tridiagonal solve at the
c ****** internal points.
c-----------------------------------------------------------------------
c
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the periodic boundary condition at X(1).
c-----------------------------------------------------------------------
c
      dxp=x(2)-x(1)
      dxm=x(n)-x(n-1)
      dxh=dxp+dxm
      a(1)=dxh*third
      c(1)=dxm*sixth
      b(1)=dxp*sixth
      fpp(1)=(f(2)-f(1))/dxp-(f(n)-f(n-1))/dxm
c
c-----------------------------------------------------------------------
c ****** Solve the (periodic) tridiagonal system for the second
c ****** derivative.
c-----------------------------------------------------------------------
c
c ****** The second derivative is returned in FPP.
c
      m=n-1

      call trid_periodic (m,c,a,b,fpp,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
        write (*,*) '### The tridiagonal matrix relating the'//
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
c
      fpp(n)=fpp(1)
c
      return
      end
c#######################################################################
      subroutine spline_periodic_type2 (n,x,f,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients for a periodic function.
c
c-----------------------------------------------------------------------
c
c ****** Get the coefficients of a cubic spline interpolant to
c ****** the function values F(i) defined at the points X(i),
c ****** i=1,...,N.
c
c ****** The computed coefficients (which are actually the second
c ****** derivatives of F at the mesh points) are returned in the
c ****** array FPP.
c
c ****** Use routine SPLINT to evaluate the spline at a
c ****** particular position.
c
c-----------------------------------------------------------------------
c
c ****** This routine assumes that the data in F is periodic, such
c ****** that F(N-1) = F(1) and F(N)  = F(2), so that two sets of
c ****** points are repeated and represent the same locations.
c ****** Thus the periodicity length is X(N-1)-X(1).
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
      real(r_typ), dimension(n) :: x,f,fpp
c
      intent(in) :: n,x,f
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ierr,m
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(2:n-1) :: a,b,c
c
c-----------------------------------------------------------------------
c
c ****** Check that there are at least 4 points.
c
      if (n.lt.4) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 4 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
c
c ****** Check that the mesh is monotonic.
c
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the coefficients for the tridiagonal solve at the
c ****** internal points.
c-----------------------------------------------------------------------
c
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
c
c-----------------------------------------------------------------------
c ****** Solve the (periodic) tridiagonal system for the second
c ****** derivative.
c-----------------------------------------------------------------------
c
c ****** The second derivative is returned in FPP.
c
      m=n-2

      call trid_periodic (m,c,a,b,fpp(2),ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
        write (*,*) '### The tridiagonal matrix relating the'//
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
c
      fpp(n)=fpp(2)
      fpp(1)=fpp(n-1)
c
      return
      end
c#######################################################################
      subroutine ezspline (n,x,f,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients.
c
c ****** Easy-to-use version of SPLINE.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls SPLINE with boundary conditions
c ****** IBC0=1 and IBC1=1.  See the comments in SPLINE to
c ****** see what this means.
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
      real(r_typ), dimension(n) :: x,f,fpp
c
      intent(in) :: n,x,f
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      integer :: ibc0,ibc1
      real(r_typ), dimension(3) :: c0,c1
c
c-----------------------------------------------------------------------
c
      ibc0=1 
      ibc1=1
c
      c0(:)=0.
      c1(:)=0.
c
      call spline (n,x,f,ibc0,c0,ibc1,c1,fpp)
c
      return
      end
c#######################################################################
      function splint (n,x,f,fpp,xv,tab)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a 1D cubic spline.
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline interpolant at X=XV.
c
c ****** On input, the function values F(i) and the second
c ****** derivatives FPP(i), defined at the points X(i),
c ****** i=1,...,N, need to be specified.
c
c ****** The optional argument TAB is an inverse interpolation
c ****** table that can be used to speed up the search for the
c ****** interval that contains XV.
c
c ****** The routine SPLINE can be used to compute FPP.
c
c ****** This routine uses routine SPEVAL to evaluate the spline.
c
c ****** The value of the spline at XV is returned as the
c ****** function value.
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
      real(r_typ), dimension(n) :: x,f,fpp
      real(r_typ) :: xv
      type(itab), optional :: tab
      real(r_typ) :: splint
c
      intent(in) :: n,x,f,fpp,xv,tab
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
c
c-----------------------------------------------------------------------
c
c ****** Find the mesh interval that encloses XV.
c
      if (present(tab)) then
        i=locate_interval(n,x,xv,tab)
      else
        i=locate_interval(n,x,xv)
      end if
c
c ****** Evaluate the cubic spline.
c
      splint=speval(x(i),x(i+1),f(i),f(i+1),fpp(i),fpp(i+1),xv)
c
      return
      end
c#######################################################################
      function speval (x1,x2,f1,f2,fpp1,fpp2,xv)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline.
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline interpolant at X=XV.
c
c ****** On input, the function values F1 and F2 and the second
c ****** derivatives FPP1 and FPP2, defined at the left and right
c ****** ends of the interval X1 and X2 need to be specified.
c
c ****** The value of the spline at XV is returned as the
c ****** function value.
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
      real(r_typ) :: x1,x2,f1,f2,fpp1,fpp2
      real(r_typ) :: xv
      real(r_typ) :: speval
c
      intent(in) :: x1,x2,f1,f2,fpp1,fpp2,xv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: dx,a,b,c,d
c
c-----------------------------------------------------------------------
c
c ****** Evaluate the cubic spline.
c
      dx=x2-x1
c
      b=(xv-x1)/dx
      a=one-b
c
      c=a*(a**2-one)*dx**2*sixth
      d=b*(b**2-one)*dx**2*sixth
c
      speval=a*f1+b*f2+c*fpp1+d*fpp2
c
      return
      end
c#######################################################################
      subroutine speval_get_abcd (x1,x2,xv,a,b,c,d)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline.
c
c ****** This version splits the calculation into two parts for
c ****** efficiency.  
c
c ****** First call SPEVAL_GET_ABCD, and then call SPEVAL_ABCD.
c
c ****** Typically SPEVAL_GET_ABCD and SPEVAL_GET_ABCD are used
c ****** when multiple spline evaluations are being performed for
c ****** the same position. 
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline interpolant at X=XV.
c
c ****** On input, the left and right limits of the spline interval
c ****** X1 and X2, as well as the position at which the spline is
c ****** being evaluated, XV, need to be specified.
c
c ****** The coefficients of the spline A, B, C, and D are
c ****** returned on output.
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
      real(r_typ) :: x1,x2,xv,a,b,c,d
c
      intent(in) :: x1,x2,xv
      intent(out) :: a,b,c,d
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: dx
c
c-----------------------------------------------------------------------
c
c ****** Evaluate the coefficients of the cubic spline.
c
      dx=x2-x1
c
      b=(xv-x1)/dx
      a=one-b
c
      c=a*(a**2-one)*dx**2*sixth
      d=b*(b**2-one)*dx**2*sixth
c
      return
      end
c#######################################################################
      function speval_abcd (a,b,c,d,f1,f2,fpp1,fpp2)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline.
c
c ****** This version splits the calculation into two parts for
c ****** efficiency.  
c
c ****** First call SPEVAL_GET_ABCD, and then call SPEVAL_ABCD.
c
c ****** Typically SPEVAL_GET_ABCD and SPEVAL_GET_ABCD are used
c ****** when multiple spline evaluations are being performed for
c ****** the same position. 
c
c-----------------------------------------------------------------------
c
c ****** On input, the coefficients A, B, C, and D, and the
c ****** function values F1 and F2 and the second
c ****** derivatives FPP1 and FPP2, need to be specified.
c
c ****** The value of the spline is returned as the
c ****** function value.
c
c ****** The coefficients A, B, C, and D can be obtained using
c ****** routine SPEVAL_GET_ABCD.
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
      real(r_typ) :: a,b,c,d,f1,f2,fpp1,fpp2
      real(r_typ) :: speval_abcd
c
      intent(in) :: a,b,c,d,f1,f2,fpp1,fpp2
c
c-----------------------------------------------------------------------
c
c ****** Evaluate the cubic spline.
c
      speval_abcd=a*f1+b*f2+c*fpp1+d*fpp2
c
      return
      end
c#######################################################################
      function locate_interval (n,x,xv,tab,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Locate a mesh interval.
c
c-----------------------------------------------------------------------
c
c ****** Find the interval I in table X(i), i=1,2,...,N,
c ****** that encloses the value XV, i.e., such that
c ****** X(I).le.XV.le.X(I+1).
c
c ****** For the special case when N=1, XV must equal X(1)
c ****** exactly, otherwise an error occurs.
c
c ****** The optional argument TAB is an inverse interpolation
c ****** table that can be used to speed up the search for the
c ****** interval.
c
c ****** If the optional argument IERR is specified, then this
c ****** routine will return when an error occurs with IERR=1.
c ****** If no error occurs, IERR=0 is returned.  When IERR is not
c ****** specified, this routine will terminate the program
c ****** with a printed error message.
c
c ****** The mesh interval I is returned as the function value.
c
c-----------------------------------------------------------------------
c
      use number_types
      use invint_def
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
      type(itab), optional :: tab
      integer, optional :: ierr
      integer :: locate_interval
c
      intent(in) :: n,x,xv,tab
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,ig
      real(r_typ) :: xi,fiv,alpha
c
c-----------------------------------------------------------------------
c
      if (present(ierr)) then
        ierr=0
      end if
c
c ****** For the special case when the table has only one
c ****** point (N=1), the inverse table is not used.  In this
c ****** case it is necessary for XV to equal X(1) exactly,
c ****** otherwise this routine exits with an error.
c
      if (n.eq.1) then
        if (xv.eq.x(1)) then
          locate_interval=i
          return
        else
          go to 900
        end if
      end if
c
c ****** Search for the interval depending on whether the optional
c ****** inverse interpolation table TAB was specified.
c
      if (.not.present(tab)) then
c
c ****** Search without an inverse interpolation table.
c
        do i=1,n-1
          if (xv.ge.x(i).and.xv.le.x(i+1)) then
            locate_interval=i
            return
          end if
        enddo
c
      else
c
c ****** Search with an inverse interpolation table.
c
c ****** Get an estimate of the nearest grid point location in
c ****** the (uniform) inverse interpolation table.
c
        xi=one+(xv-x(1))*tab%d
        i=xi
        i=max(i,1)
        i=min(i,tab%n-1)
        alpha=xi-i
        fiv=(one-alpha)*tab%f(i)+alpha*tab%f(i+1)
c
c ****** Set IG to be the guess for the nearest grid point.
c
        ig=fiv
        ig=max(ig,1)
        ig=min(ig,n-1)
c
        if (xv.ge.x(ig)) then
c
c ****** Search forwards.
c
          do i=ig,n-1
            if (xv.ge.x(i).and.xv.le.x(i+1)) then
              locate_interval=i
              return
            end if
          enddo
c
        else
c
c ****** Search backwards.
c
          do i=ig-1,1,-1
            if (xv.ge.x(i).and.xv.le.x(i+1)) then
              locate_interval=i
              return
            end if
          enddo
c
        end if
c
      end if
c
  900 continue
c
c ****** Value not found.
c
c ****** If IERR was passed, set IERR=1 and return; otherwise,
c ****** write an error message and terminate the program.
c
      if (present(ierr)) then
        ierr=1
        return
      else
        write (*,*)
        write (*,*) '### ERROR in LOCATE_INTERVAL:'
        write (*,*) '### The value requested was not found in'//
     &              ' the table:'
        write (*,*) 'Value requested = ',xv
        write (*,*) 'Minimum table value = ',x(1)
        write (*,*) 'Maximum table value = ',x(n)
        call exit (1)
      end if
c
      end
c#######################################################################
      subroutine trid (n,c,a,b,d,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Solve the tridiagonal system of equations:
c
c         C(i)*X(i-1) + A(i)*X(i) + B(i)*X(i+1) = D(i)
c
c        for i=2,...,N-1, with
c
c           A(1)*X(1) + B(1)*X(2) = D(1)
c
c        and
c
c           C(N)*X(N-1) + A(N)*X(N) = D(N)
c
c ****** Note that C(1) and B(N) are not referenced.
c
c ****** D is overwritten with the solution.
c
c ****** Return IERR=0 for a successful completion.  If the
c ****** matrix is singular, this routine returns IERR=1 and
c ****** D is invalid.
c
c ****** This routine does not do any pivoting, so the solution
c ****** is not guaranteed to be accurate unless the matrix is
c ****** well-conditioned (e.g., diagonally dominant).
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
      real(r_typ), dimension(n) :: c,a,b,d
      integer :: ierr
c
      intent(in) :: n,c,a,b
      intent(inout) :: d
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
      real(r_typ) :: denom,ace
      real(r_typ), dimension(n) :: aa
c
c-----------------------------------------------------------------------
c
      ierr=1
c
c ****** Copy A to AA, since it will be overwritten during the
c ****** elimination.  This prevents A from being overwritten.
c
      aa=a
c
c ****** Forward elimination.
c
      if (aa(1).eq.0.) return
      d(1)=d(1)/aa(1)
      aa(1)=b(1)/aa(1)
      do i=2,n
        denom=aa(i)-c(i)*aa(i-1)
        if (denom.eq.0.) return
        ace=one/denom
        if (i.ne.n) aa(i)=ace*b(i)
        d(i)=ace*(d(i)-c(i)*d(i-1))
      enddo
c
c ****** Backward substitution.
c
      do i=n-1,1,-1
        d(i)=d(i)-aa(i)*d(i+1)
      enddo
c
c ****** Set the error return flag to indicate successful completion.
c
      ierr=0
c
      return
      end
c#######################################################################
      subroutine trid_periodic (n,c,a,b,d,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Solve the tridiagonal system of equations:
c
c         C(i)*X(i-1) + A(i)*X(i) + B(i)*X(i+1) = D(i)
c
c        for i=2,...,N-1, with
c
c           C(1)*X(N) + A(1)*X(1) + B(1)*X(2) = D(1)
c
c        and
c
c           C(N)*X(N-1) + A(N)*X(N) + B(N)*X(1) = D(N)
c
c ****** D is overwritten with the solution.
c
c ****** Return IERR=0 for a successful completion.  If the
c ****** matrix is singular, this routine returns IERR=1 and
c ****** D is invalid.
c
c ****** This routine does not do any pivoting, so the solution
c ****** is not guaranteed to be accurate unless the matrix is
c ****** well-conditioned (e.g., diagonally dominant).
c
c ****** This system arises for periodic solutions.
c
c ****** This routine uses the Sherman-Morrison formula for
c ****** updating a matrix inverse with a low-rank modification.
c ****** The modification arises from the changes introduced by the
c ****** periodic boundary conditions.
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
      real(r_typ), dimension(n) :: c,a,b,d
      integer :: ierr
c
      intent(in) :: n,c,a,b
      intent(inout) :: d
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
      real(r_typ) :: denom,ace
      real(r_typ), dimension(n) :: aa
      real(r_typ), dimension(n) :: y
      real(r_typ), dimension(n,2) :: z2
      real(r_typ), dimension(2,2) :: t,tinv
      real(r_typ) :: detinv
      real(r_typ), dimension(2) :: tvty
c
c-----------------------------------------------------------------------
c
      ierr=1
c
c ****** First, solve the (non-periodic) system with an RHS
c ****** equal to D.
c
      y=d
c
      call trid (n,c,a,b,y,ierr)
c
      if (ierr.ne.0) return
c
c ****** Next, solve the two systems for the "inhomogenous part".
c
      z2=0.
      z2(1,1)=one
      z2(n,2)=one
c
c ****** Copy A to AA, since it will be overwritten during the
c ****** elimination.  This prevents A from being overwritten.
c
      aa=a
c
c ****** Forward elimination.
c
      if (aa(1).eq.0.) return
      z2(1,:)=z2(1,:)/aa(1)
      aa(1)=b(1)/aa(1)
      do i=2,n
        denom=aa(i)-c(i)*aa(i-1)
        if (denom.eq.0.) return
        ace=one/denom
        if (i.ne.n) aa(i)=ace*b(i)
        z2(i,:)=ace*(z2(i,:)-c(i)*z2(i-1,:))
      enddo
c
c ****** Backward substitution.
c
      do i=n-1,1,-1
        z2(i,:)=z2(i,:)-aa(i)*z2(i+1,:)
      enddo
c
c ****** Invert the 2 x 2 system.
c
      t(1,1)=one+z2(n,1)*b(n)
      t(1,2)=z2(n,2)*b(n)
      t(2,1)=z2(1,1)*c(1)
      t(2,2)=one+z2(1,2)*c(1)
c
      denom=t(1,1)*t(2,2)-t(1,2)*t(2,1)
      if (denom.eq.0.) return
      detinv=one/denom
c
      tinv(1,1)= detinv*t(2,2)
      tinv(2,2)= detinv*t(1,1)
      tinv(1,2)=-detinv*t(1,2)
      tinv(2,1)=-detinv*t(2,1)
c
c ****** Construct the final periodic solution.
c
      tvty(1)=tinv(1,1)*b(n)*y(n)+tinv(1,2)*c(1)*y(1)
      tvty(2)=tinv(2,1)*b(n)*y(n)+tinv(2,2)*c(1)*y(1)
c
      d(:)=y(:)-z2(:,1)*tvty(1)-z2(:,2)*tvty(2)
c
c ****** Set the error return flag to indicate successful completion.
c
      ierr=0
c
      return
      end
