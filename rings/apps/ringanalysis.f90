!  Subroutine to combine cart_to_polar, fourier_filter, and hmifit
!  in one program for porting to HMI -- only carries out single mode fit.
!  
!----------------------------------------------------------------
!  
! Call cart_to_polar with following inputs and outputs:
!
! input:  powxy, nkx,nky,nnu,ntht,nk,lnbase,verbose,crpix1,crpix2
!        
! output: powthtk(ntht,nk,nnu) real*4 array (kind=single)
!          ierr    = 0 then subroutine performed correctly.
!                  /= 0 then ierr specifies where it failed 
!-----------------------------------------------------------------              
!
! Call fourier_filter with following inputs and outputs:
!
! input:  powthtk, ntht,nk,nnu,nthts,dnu,mdata,mfilt,verbose
!        
! output: powfilt(nthts,nk,nnu) real*8 array,
!         filter(nthts,nk) real*4 array
!         inumn, inumx
!          ierr    = 0 then subroutine performed correctly.
!                  /= 0 then ierr specifies where it failed 
!-----------------------------------------------------------------
!  Input Constants for Tile to be Fit------------------------------------------
!
!  nkx = nky : number of pixels in x and y dimensions 
!            (assumed to be the same) 
!  nk = number of k values, should be nkx/2
!  nnu = nw : number of frequencies in power spectrum
!  dk = delta_k in Mm^(-1), i.e. the wavenumber per bin
!  dnu = delta_nu in micro Hz, i.e. the frequency per bin
!  ntht = number of thetas for cart_to_polar to unwrap power spectra to
!                  = approx 2*pi*nk
!  nthts = number of thetas after filtering and sub sampling, usually
!                smaller than ntht by a factor of 10, should be divisible by 4
!  lnbase = flag to tell if powxy is ln(power) or just power
!         = real value of e if ln or 0.0 if just power
!  nmin = minimum n value to fit for
!  nmax = maximum n value to fit for
!  ux_guess = starting guess for fit to ux (real*8)
!  uy_guess = starting guess for fit to uy (real*8)
!  kbstrt = array containing the start bins for carrying out fits for each
!		ridge dimensions(nrdk)
!  kbend = array containing which bins to stop at when fitting each ridge
!               dimensions(nrdk)
!  nrdk = size of arrays kbstrt, kbend
!  ma = max # of parameters possible for pow, freq., width, etc. in guessfile
!  np = number of parameters to be fit,(used to be np) currently = 6, 
!            will change for multiridgefit code
!  ikm = amount of smoothing in k to take place for fitting
!         e.g., if ikm=1 then the program will fit the range between k-1 : k+1
!  header = header string for output file
!  strln_* = length of various file names in bytes
!  rgn = record count in drms series 
!  Input File Names ---------------------------------------------------------
!
!  guessfile = file that holds the parameterized guesses
!                  to be used in fitting routine
!
!   Input Flags --------------------------------------------------------------
!
!     doptions(1) = xtol   :  double precision
!       Termination tolerance on X.
!
!     doptions(2) = ftol   :  double precision
!       Termination tolerance on F.
!
!     doptions(3) = feps   :  double precision
!       (Approximately) The accuracy of the function values
!       calculated in FUNC.
!
!     doptions(4) = delta0 :  double precision
!       Initial radius of the trust region.
!
!     ioptions(1) = maxeval  :  integer
!       On input  : Maximum number of function evaluations.
!
!     ioptions(2) = igrad  :  integer
!       If igrad is 0 then F' is calculated numerically.
!       If igrad is 1 then the user supplied subroutine DF1 is called
!       by subroutine dogleg from subroutine hmifits to calculate F'.
!       If igrad is 2 then the user supplied subroutine DF1 is called
!       to calculate F' and initially the results from GRAD are checked
!       by comparing them to the F' obtained using finite differences.
!
!     ioptions(3) = iextrap  :  integer
!       If iextrap is 1 then safeguarded quadratic/cubic interpolation
!       and extrapolation is performed.
!
!     ioptions(4) = idiag  :  integer
!       If idiag is 0 then I is used as the initial approximation to the
!       Hessian.
!       If idiag is 1 then the diagonal of the true Hessian is computed
!       initially using a finite difference approximation.
!
!     ioptions(5) = iprint  :  integer
!       If iprint=1 then information is printed after each iteration.
!
!     verbose = 0 or 1 : 0 = no extra info, 1 = write out more variables
!               for debugging
!
!  Calculated variables--------------------------------------------------
!
!  inumn =  integer frequency bin number for 1.5 mHz 
!  inumx =  integer frequency bin number for 5.5 mHz 

!  Output arrays -------------------------------------------------------
!
!  powfilt(nthts,nk,nnu)  = filtered power spectrum (filtered_pspec in wrapper)
!  powthtk(ntht,nk,nnu) = unwrapped power spectrum (unwrapped_pspec in wrapper)
!  filter(nthts,nk) = filter applied to subsampled powthtk to make powfilt
!   
      SUBROUTINE ringanalysis (powxy, nkx, nky, nnu, ntht, nthts, nk, &
           crpix1, crpix2, ikm, dk, dnu, lnbase, nmin, nmax, np, &
           ntot1, ux_guess, uy_guess, ux, d_ux, uy, d_uy, amp, d_amp,  &
           bg, d_bg, fwhm, d_fwhm, anu, d_nu, delnu, nval, kval, kbin, &
           kbstrt, kbend, kbmin, kbmax, nrdk,  &
           min_func, n_iter, fit_chk, chi, kmin, kmax, modecount, good, &
           powthtk, filter, powfilt, fnu, width, ampg, bkg, verbose, mdata, &
           mfilt, doptions, ioptions, flag, rgn, nrdtot, ierr)
     use cart
! don't need this if you aren't using write_data: use wt_data
     use fourier
     use fits
     implicit none
     integer, parameter :: single = 4
     integer, parameter :: double = 8
     integer :: nkx, nky, nnu, ntht, nthts, nk, nmin, nmax, mdata 
     integer :: verbose, inumn, inumx, np, flag, ntot1
     integer :: ikm, ierr,  mfilt, ik, nnj, rgn, nrdtot
     integer :: i, j, k, inu3, ik3, ik2, ijk, kbmin, kbmax
     integer :: nrdk, modecount
     integer, dimension(5) :: ioptions
     integer, dimension(nrdk) :: kbstrt, kbend
     integer, dimension(nk,ntot1) :: kbin, good, n_iter, fit_chk, nval
     real(kind=single) :: lnbase, crpix1, crpix2, kmin, kmax
     real(kind=single), dimension(nkx,nky,nnu) :: powxy
     real(kind=single), dimension(ntht,nk,nnu) :: powthtk
     real(kind=single), dimension(nthts,nk) :: filter
     real(kind=single), dimension(nthts,nk,nnu) :: powfilt
     real(kind=single), dimension(nk,ntot1) :: ux, uy, d_ux, d_uy
     real(kind=single), dimension(nk,ntot1) :: kval, amp, d_amp, chi, min_func
     real(kind=single), dimension(nk,ntot1) :: fwhm, d_fwhm, bg, d_bg 
     real(kind=single), dimension(nk,ntot1) :: anu, d_nu, delnu 
     real(kind=double), dimension(nk,nrdtot) :: fnu, width, ampg, bkg
     real(kind=double) :: dk, dnu 
     real(kind=double) :: ux_guess, uy_guess
     real(kind=double), dimension(4) :: doptions
     character(len=132) :: filou, filouf

!--------------------------------------------------------------------
     IF (verbose == 1) THEN
        write(*,'(" ring, flag = ",i15)')flag
        write(*,*)nnu, nk, nkx, nky, ntht, nthts
        write(*,*) nrdtot
        inu3=nnu/3
        write(*,'(" rgn =",i5)') rgn
        WRITE(*,10) ((powxy(nk+5,j,k),j=nk+5,nk+7),k=inu3,inu3+2)
  10    format("powxy ",/,3(3pe15.5,2x,/))
        write(*,'(" ring, doptions= ",4(d15.5,2x))')(doptions(i),i=1,4)
        write(*,'(" ring, ioptions= ",5i7)')(ioptions(i),i=1,5)
        write(*,'(" kbstrt = ",<nrdk>(i4))') (kbstrt(i),i=1,nrdk)
        write(*,'(" kbend = ",<nrdk>(i4))') (kbend(i),i=1,nrdk)
        write(*,*) kbmin, kbmax
        write(*,'("in ringanalysis, record count = ",i6)')rgn
!        write(*,'("fnu ",<nmax>(1pd13.3,2x))')((fnu(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
!        write(*,'("width ",<nmax>(1pd13.3,2x))')((width(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
!        write(*,'("amp ",<nmax>(1pd13.3,2x))')((ampg(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
!        write(*,'("bkg",<nmax>(1pd13.3,2x))')((bkg(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
     ENDIF

!-------------------------------------------------------------------------- 
! Carry out unwrapping of power spectrum powxy by interpolating 
!         spectrum powxy to polar coordinates, theta and k, yielding powthtk.
!
!---------------------------------------------------------------------------
      ierr=0
      crpix1=crpix1-0.5
      crpix2=crpix2-0.5
      if (verbose .eq. 1)then
        write(*,'("crpix1,2 = ",2f8.4)') crpix1, crpix2
      endif
      CALL cart_to_polar(powxy,powthtk,nkx,nky,nnu,ntht,nk,crpix1, &
     &                   crpix2,lnbase,verbose,ierr)

      IF (ierr /= 0) RETURN

!-----------------------------------------------------------
!
!  TAKE THESE LINES OUT WHEN INTEGRATED INTO DRMS
!
!      if (verbose == 1) then
!         ik3 = nk/3
!         inu3 = nnu/3
!         write(*,'("powthtk",/,5(1pe13.4,2x))') &
!     &              ((powthtk(1,i,j),i=ik3,ik3+4),j=inu3,inu3+3)
!         filou="powthtk.fits"
!         call write_data(filou,powthtk,ntht,nk,nnu)
!      endif
!-----------------------------------------------------------
!
! Calculate Fourier filter of unwrapped spectrum and subsample in subroutine
! fourier_filter 
!
!----------------------------------------------------------------

     CALL fourier_filter(powthtk,powfilt,filter,ntht,nk,nnu,nthts,dnu, &
                         mdata,mfilt,inumn,inumx,verbose,ierr)
      
      IF (ierr /= 0) RETURN
!-----------------------------------------------------------
!
!  TAKE THESE LINES OUT WHEN INTEGRATED INTO DRMS
!
!      if (verbose == 1) then
!         filouf="powfilt.fits"
!         call write_data(filouf,powfilt,nthts,nk,nnu)
!      endif
!      
!
!  Calculate Lorentzian fits of the spectrum and the associated errors.  
!  The hmifits subroutine writes the fits, the errors, and the goodness
!  of fit parameters specified by hmi ring team to an 
!  ascii file (outfile) to be used for input to the inversion routines
!  

      CALL hmifits (powfilt, nthts, nk, nnu, nmin, nmax, ntot1, dk, dnu, np, &
                 inumn, inumx, ux_guess, uy_guess, kbstrt, kbend, doptions, &
                 ioptions, nrdk, ikm, & 
                 ux, d_ux, uy, d_uy, amp, d_amp, bg, d_bg, fwhm, d_fwhm, & 
                 anu, delnu, d_nu, nval, kval, kbin, modecount, good,  &
                 min_func, n_iter, fit_chk, chi, kmin, kmax, rgn, fnu, &
                 width, ampg, bkg, nrdtot, kbmin, kbmax, verbose, ierr, flag)

 
      IF (ierr /= 0) RETURN
!      if (verbose == 1) then
!        write(*,'(" anu = ", <nk>(1pe12.4))')((anu(ik,nnj),ik=1,nk),nnj=1,nmax+1)
!        write(*,'(" ux = ", <nk>(1pe12.4))')((ux(ik,nnj),ik=1,nk),nnj=1,nmax+1)
!      endif
! 
!  Return to DRMS
!
      RETURN
!
!  This point in program is only reached if there has been an error
!   in any part of the processing
!
!9999  IF (ierr /= 0) THEN
!        IF ((ierr >= 10) .AND. (ierr <= 20)) THEN
!          WRITE(*,'(" Problem unwrapping power spectrum = ",i3)') ierr
!          write(*,'(" Culprit in cart_to_polar subroutine ")')
!        ENDIF
!        IF ((ierr >= 20) .AND. (ierr <= 30)) THEN
!          WRITE(*,'(" Problem subsampling and filtering spectrum = ", &
!     &                i3)') ierr
!        ENDIF
!        IF (ierr >= 30) THEN
!          WRITE(*,'(" Problem fitting power spectrum = ",i3)') ierr
!        ENDIF
!      ENDIF
!      RETURN
      END subroutine
