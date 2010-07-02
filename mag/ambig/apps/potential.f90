!***********************************************************************
! subroutine potential                                                 *
!    Use FFTW routines called through the Math Kernel Library to       *
! compute the potential field and some of its gradients.               *
!***********************************************************************

subroutine potential(nxp,nyp)

!-----------------------------------------------------------------------
   Use MKL_DFTI
   use sizes
   use constant
   use pix_size
   use pot_field
   use pad
   use trnsfrm

   implicit none

   integer :: i,j,i1d,im,jm,im1d,nxp,nyp
   type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle
   integer :: Status
   integer :: lengths(2)
   integer :: strides_in(3)
   integer :: strides_out(3)
   real :: bz0,bl0,sigma,xkappa
   real,dimension(:),allocatable :: xi,yi
   complex :: gradlos
   complex,dimension(:),allocatable :: blpad1d,fftbl1d
   complex,dimension(:),allocatable :: db3dl1d,fftdb31d
   complex,dimension(:),allocatable :: bipx1d,bipy1d
   complex,dimension(:),allocatable :: fftbix1d,fftbiy1d
!------------------------------------------------------------------------
! --> Subtract uniform vertical field to achieve flux balance.

   bl0=0.
   do i=1,nxp
      do j=1,nyp
         bl0=bl0+blpad(i,j)
      enddo
   enddo
   bl0=bl0/float(nxp)/float(nyp)
   bz0=bl0/a33

! --> Convert 2D real array of line of sight field to 1D complex array. 

   allocate(blpad1d(nxp*nyp),fftbl1d(nxp*nyp))
   do i=1,nxp
      do j=1,nyp
         i1d=i+(j-1)*nxp
         blpad1d(i1d)=cmplx(blpad(i,j)-bl0,0.)
      enddo
   enddo
   deallocate(blpad)

! --> Fourier transform line of sight component of the field.

   lengths(1) = nxp
   lengths(2) = nyp

   strides_in(1) = 0
   strides_in(2) = 1
   strides_in(3) = nxp

   strides_out(1) = 0
   strides_out(2) = 1
   strides_out(3) = nxp

! --> Create MKL descriptor for 2D real to complex transform
   Status = DftiCreateDescriptor(Desc_Handle, DFTI_SINGLE, DFTI_COMPLEX, 2, lengths)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if

   Status = DftiSetValue(Desc_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if

   Status = DftiSetValue(Desc_Handle,DFTI_INPUT_STRIDES,strides_in)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if

   Status = DftiSetValue(Desc_Handle,DFTI_OUTPUT_STRIDES,strides_out)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if

! --> Set the scale
   sigma=1./(float(nxp)*float(nyp))
   Status = DftiSetValue(Desc_Handle, DFTI_FORWARD_SCALE, sigma)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if

   Status = DftiCommitDescriptor( Desc_Handle )
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if

! --> Compute 2D complex to complex transform
   Status = DftiComputeForward(Desc_Handle,blpad1d,fftbl1d)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if
   deallocate(blpad1d)

! --> Wavenumbers (recall the order in which FFTW returns the frequencies)
   allocate(xi(nxp),yi(nyp))

   do i=1,nxp/2+1
      xi(i)=float(i-1)/float(nxp)/dxi
   enddo
   do i=1,(nxp-1)/2
      xi(nxp-i+1)=-xi(i+1)
   enddo

   do j=1,nyp/2+1
      yi(j)=float(j-1)/float(nyp)/dyi
   enddo
   do j=1,(nyp-1)/2
      yi(nyp-j+1)=-yi(j+1)
   enddo

! --> Multiply by appropriate functions of wavenumber and coordinate
! --> transform. 
   allocate(fftbix1d(nxp*nyp),fftbiy1d(nxp*nyp))
   allocate(fftdb31d(nxp*nyp))

   do i=1,nxp
      do j=1,nyp
         i1d=i+(j-1)*nxp
         xkappa=2.*pi*sqrt(((c11**2+c12**2)*xi(i)**2+2.*(c11*c21+c12*c22)*xi(i)*yi(j)+(c21**2+c22**2)*yi(j)**2))
         gradlos=cmplx(-a33*xkappa,2.*pi*((c11*a13+c12*a23)*xi(i)+(c21*a13+c22*a23)*yi(j)))

! --> Transverse components of potential field.
         fftbix1d(i1d)=fftbl1d(i1d)*cmplx(xi(i)*(c11**2+c12**2)+yi(j)*(c11*c21+c12*c22),xkappa*a31/(2.*pi))/gradlos
         fftbiy1d(i1d)=fftbl1d(i1d)*cmplx(xi(i)*(c11*c21+c12*c22)+yi(j)*(c21**2+c22**2),xkappa*a32/(2.*pi))/gradlos

! --> Vertical derivative of vertical component.
         fftdb31d(i1d)=fftbl1d(i1d)*xkappa**2/gradlos

      enddo
   enddo
   deallocate(xi,yi,fftbl1d)

! --> Force constant term to be 0.
   fftbix1d(1)=cmplx(0.,0.)
   fftbiy1d(1)=cmplx(0.,0.)

! --> Vertical derivative also needs constant term to be 0.
   fftdb31d(1)=cmplx(0.,0.)

! --> Inverse transform to get potential field components

! --> Compute 2D backward transform for x-component
   allocate(bipx1d(nxp*nyp))
   Status = DftiComputeBackward(Desc_Handle,fftbix1d,bipx1d)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if
   deallocate(fftbix1d)

! --> Compute 2D backward transform for y-component
   allocate(bipy1d(nxp*nyp))
   Status = DftiComputeBackward(Desc_Handle,fftbiy1d,bipy1d)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if
   deallocate(fftbiy1d)

! --> Inverse transform to get vertical gradient of the potential field

! --> Compute 2D backward transform for z-component
   allocate(db3dl1d(nxp*nyp))
   Status = DftiComputeBackward(Desc_Handle,fftdb31d,db3dl1d)
   if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
       call dfti_example_status_print(Status)
       goto 999
   end if
   deallocate(fftdb31d)

! --> Extract field and line of sight gradients

   do i=1,nxp
      im=i
      do j=1,nyp
         jm=j
         im1d=im+(jm-1)*nxp

! --> Transverse components of potential field.
         Bpix(i,j)=bz0*a31-2.*pi*aimag(bipx1d(im1d))
         Bpiy(i,j)=bz0*a32-2.*pi*aimag(bipy1d(im1d))

! --> Vertical derivative of vertical component.
         dBpzdz(i,j)=real(db3dl1d(im1d))

      enddo
   enddo
   deallocate(bipx1d,bipy1d)
   deallocate(db3dl1d)

  999 return
end subroutine potential
