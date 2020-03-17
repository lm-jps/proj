! module fits contains
!   Subroutine to carry out actual single ridge fits.  Writes fits and errors to outfile.
!   Subroutine hmifits unwrapped rings using  hybrid quasi-Newton technique. 
!     Subroutine dogleg developed by Rasmus Larsen fitting all parameters, 
!     none are held constant.
!     fnua, fwa, and pa are the coefficients of the polynomials for frequencies,
!     widths, and power respectively and are used as first guesses for the
!     fit.  Changed to be able to use first derivative * first derivative hessian 
!           (no second derivs) on Jan. 10, 2007, 
!     If flag=0 use 2nd derivative hessian, if flag=1 just use first derivatives
!     
!------------------------------------------------------------------------------------
!
!  Book Keeping and Dimension Variables-------------------------------------------- 
!
!    nthts, nk    The number of data elements in the theta and k dimensions 
!                   of the subsampled polar coordinate system,  (input)
!
!    nnu          Number of data elements in the frequency dimension (input) 
!
!    ik,itht,inu  Index variables for the appropriate data dimensions. (book keeping)
!
!    ntot         Number of ridges (n values) to be fit
!
!  Variables Determined by the Duration and Size of Tracked Data---------------------------------------------
!        
!    dk           Delta k in inverse Mm.  Keyword read from data (input) REAL*8
!                   determined by size of tracked tile 
!    dnu          Delta nu in microHz. Keyword read from data (input) REAL*8
!                   determined by length of observing run
!       
!  Calculated Variables------------------------------------------------------------------------------------
!
!  emult     Variable by which pow is divided to reduce problems in fitting due to 
!              precision errors.
!  
!  elnmult   Logarithm of emult.
!
!  rdchi     Reduced chisquare test
!
!  ell       Value of the radial order for a given wavenumber       
!
! 
!  File Reading and Writing Variables----------------------------------------------------
!
!       strln_guess      Length of character strings: for guess file, 
!
!       guessfile  File containing polynomial coefficients for determining 
!                    first guesses to ridge frequencies, widths, power, and background.  (input) 
!
!
!  Other inputs------------------------------------------------------------------
!
!    np  number of parameters in fit
!  
!    ikm smoothing parameter = smooth in calling routine = number of k bins on
!         either side of k bin of interest
!
!    ux_guess   initial velocity guess in the x direction
!    uy_guess   initial velocity guess in the y direction
!    kbstrt	starting k values for each ridge
!    kbend	stopping k values for each ridge
!    nrdk	size of kbstrt, kbend
!
!
!  Arrays--------------------------------------------------------------------------------
!
!  pow    Subsampled and filtered power spectrum.  Dimensions(nthts,nk,nnu)  (input)
!            Corresponds to powfilt in main program.         
!  
!  ak          ik*dk, i.e. Wavenumbers in units of Mm^(-1).  Dimensions(nk). Real*8
!
!  anu         inu*dnu i.e. Frequencies in units of micro Hz. Dimensions(nnu). Real*8
!
!  s2      second derivative of function - used for calculation of errors
!             dimensions(np,np)
!  
!  d_pj    corresponding errors for parameters pj(np) 
!
!
!  Routines and Functions -----------------------------------------------------------------------
!
!   Dogleg       Hybrid quasi-Newton fitting technique
!   
!   dpotrf, dpotri   Subroutines to carry out matrix inversion on Hessian
! 
!   FUNK1, FUNK  Routines where the chi-square with 2 degrees of freedom function gets
!                calculated
!   DF1, DF      Where first derivative function of FUNK gets calculated.
!
!   S2           Where second derivative (Hessian) is explicitly calculated. 
!                Inverted to give formal errors.
! 
!  Book Keeping Variables and Constants-----------------------------------------------------
!
!   ni, nj    Used for stepping through radial orders (n values). 
!   npj       total number of data points being fit
!-------------------------------------------------------------------------------------
    module fits
      contains
      SUBROUTINE hmifits (powfilt ,nthts, nk, nnu, nmin, nmax, ntot1, dkn, &
         dnun, np, inumn, inumx, ux_guess, uy_guess, kbstrt, kbend, doptions,  &
         ioptions, nrdk, ikm, ux, d_ux, &
         uy, d_uy, ap, d_ap, bg, d_bg, fwhm, d_fwhm, anu, delnu, d_nu, nval, &
         kval, kbin, modecount, good, min_func, n_iter, fit_chk, chi, kmin, &
         kmax, rgn, fnu, width, amp, bkg, nrdtot, kbmin, kbmax, verbose, ierr, flag)
   
      use ring_pass              
      use factor
      use func1
      use first_deriv_1
      use second_deriv
      IMPLICIT NONE
      integer, parameter :: single=4 
      integer, parameter :: double=8
      integer :: nthts,nk,nnu,ik,itht,inu,nmin,nmax,verbose
      integer :: ikm, rgn
      real(kind=double) :: ux_guess,uy_guess
      integer :: ier,ierr,np,inumn,inumx
      integer :: nkeys
      integer :: kbstn,kbendn,istrtnu,iendnu,flag,fitflag
      integer :: fpix, grp, i , j, maxeval, neval
      integer :: nrdk, nnj, nrdtot, inuave, nk2, ntot1, nk3, nnu3
      integer, dimension(5) :: ioptions
      integer, dimension(nrdk) :: kbstrt, kbend
      integer info,chng,kbmin,kbmax,modecount
      integer, dimension(nk,ntot1) :: good, n_iter, kbin, nval, fit_chk
      real(kind=single), dimension(nk,ntot1) :: d_ux,d_uy,d_ap,d_bg,d_fwhm,d_nu
      real(kind=single), dimension(nk,ntot1) :: chi,fwhm,ap,bg,delnu,ux,uy
      real(kind=single), dimension(nk,ntot1) :: anu, min_func, kval
      real(kind=single) :: kmin, kmax
      real(kind=single), dimension(nthts,nk,nnu) :: powfilt
      real(kind=double), dimension(nthts,nk,nnu) :: pow
      real(kind=double) :: alk,ak,akstm,akend,rchi,rdchi,dkn,dnun
      real(kind=double) :: frng,emult,elnmult,ftol,pi,fmode
      real(kind=double) :: xtol,y1,dtht,theta,twopi,gsun,widmin,widmax
      real(kind=double), dimension(np) :: pj
      real(kind=double), dimension(np,np) :: d_pj, s2
      real(kind=double), dimension(4) :: doptions
      real(kind=double),  dimension(nk,nrdtot) :: width,amp,bkg,fnu
      character(len=1) :: uplo

! 
! initialize modecount array, good counter, and other output arrays
!
        modecount = 0
        good = 0
        ux = 0 
        uy = 0 
        ap = 0 
        bg = 0 
        delnu = 0 
        fwhm = 0 
        anu = 0 
        d_ux = 0 
        d_uy = 0 
        d_ap = 0 
        d_bg = 0 
        d_nu = 0 
        d_fwhm = 0
        nval = 0
        kval = 0
        kbin = 0
        min_func = 0
        n_iter = 0
        fit_chk = 0
        chi = 0

!
!----------------------------------------------------------------
! set up variables to be passed to funk/funk1 via the module ring_pass
!
      v = verbose
      dk = dkn
      dnu = dnun
      ntot = ntot1
      twopi = 2d0*dacos(-1d0)
      dtht = twopi/dble(nthts)
      if (verbose == 1)then
        write(*,'(" hmifits, flag = ",i5)')flag
        write(*,'(" hmifits, doptions= ",4(es15.5,2x))')(doptions(i),i=1,4)
        write(*,'(" hmifits, ioptions= ",5i7)')(ioptions(i),i=1,5)
      endif
  
!
!  Function multfactor calculates and returns a division factor for data 
!  that keeps the maximum of the log(amplitude) at 9.5. This helps
!  keep the division of numerator and denominator from depending too much
!  on machine precision. It returns the factor needed as a real*8 number. 
!  The log of this factor will get added back to the amplitude and background
!  terms before they get written out (or could go in comment section so people
!  know what to add back in themselves).  Since the power spectrum shows
!  up in the maximum likelihood function as 1/power, the factor is actually 
!  multiplied by this term.
!
      ierr=0
      
!      emult = multfactor (powfilt, nthts, nk, nnu, inumn, inumx, &
!                         kbstrt, kbend, kbmin, kbmax, verbose, ierr)
      emult = multfactor (powfilt, nthts, nk, nnu, inumn, inumx, kbmin, &
                          kbmax, verbose, ierr)

      if (ierr .ne. 0) goto 9999
      elnmult = dlog(emult)
      
!
!  Use multiplication factor emult to scale power.
!
!     inuave = (inumx - inumn)/2
     nk3 = nk / 3
     nnu3 = nnu / 3
     if (verbose == 1) then
!      if (inuave < nnu .and. 5 < nthts) then
       write(*,'("powfilt in hmifits",8(1pd13.5,2x))')((powfilt(1,ik,inu), &
                  ik=nk3, nk3+2),inu=nnu3, nnu3+2)
!      endif
     endif
      where (powfilt /= 0.0)
        pow=emult/dble(powfilt)
      end where
     if (verbose == 1) then
!      if (inuave < nnu .and. 5 < nthts) then
       write(*,'("pow in hmifits",8(1pd13.5,2x))')((pow(1,ik,inu), &
                  ik=nk3,nk3+2),inu=nnu3, nnu3+2)
!      endif
     endif
!
!  Calculate minimum and maximum k values to be fit, c program will translate
!   into l values so that both sets of numbers can be included
!   at the start of the fits file, kmin and kmax are in inverse Mm. 
!
      kmin = kbmin * sngl(dk)
      kmax = kbmax * sngl(dk)

      if (verbose .eq. 1) then
        write(*,'("emult = ",1pd13.4)') emult
	write(*,'("kmin = ",1pd13.4,"kmax = ",1pd13.4)') kmin,kmax
      endif
      
!
!  Set up arrays of first guesses as a function of k by reconstructing 
!  the polynomial fits at each k from the parameterizations of frequency, 
!  amplitude,width, and background found in guessfile.
! 

      ierr = 0
   if (verbose == 1) then
      write(*,'("in hmifits, record count = ",i6)')rgn
      write(*,'("fnu ",<nmax>(1pd13.3,2x))')((fnu(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
      write(*,'("width ",<nmax>(1pd13.3,2x))')((width(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
      write(*,'("amp ",<nmax>(1pd13.3,2x))')((amp(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
      write(*,'("bkg",<nmax>(1pd13.3,2x))')((bkg(ik,j),ik=kbstrt(j),kbend(j)),j=1,nmax)
   endif

!      write(34,'("# guess file used = ",a<strln_gus>)') guessfile(1:strln_gus) 
!      write(34,'("# multipLication factor used = ",1pd12.4)') emult
!      write(34,'("# nmin = ",i2,8x," nmax = ",i2)') nmin, nmax
!      write(34,'("# kmin = ",f8.4,3x,"kmax = ",f8.4,3x,"lmin = ",f10.4, &
!               & 3x,"lmax = ",f10.4)') sngl(kmin),sngl(kmax),sngl(lmin),sngl(lmax)
!      write(34,'("# delta_nu = ",1pd13.4,3x,"delta_k = ",1pd13.4)')dnu, dk
!      write(34,'("#n     l        k        nu        d_nu &
!     &        ux         d_ux         uy         d_uy        amp    & 
!     &     d_amp         bg         d_bg        fwhm       d_fwhm    & 
!     &  delnu        d_nu    k_bin  nfe   min_func      rdchi& 
!     &   fit")')
!
!  ioptions(1) starts as the maximum number of iterations allowed for each
!  fit

      maxeval = ioptions(1)
!
!  Main loop for performing fitting of specified regions of power spectrum
!
      do ni = nmin, nmax
         nj = ni+1   !take care of fortran indexing starting at 1
! 
!  Smooth over range of k's when ikm > 0
!
         kbstn = kbstrt(nj) + ikm
         kbendn = kbend(nj) - ikm 
         do ik = kbstn, kbendn
            istk = ik - ikm
            iendk = ik + ikm
            if(istk < 1) istk=1
            if(iendk > nk) iendk=nk
!
!  Calculate ak = wavenumber in Mm^{-1} and its corresponding l value 
!   to put in output file
!
            ak = dble(ik) * dk 
            if (verbose .eq. 1) then
              write(*,'("n =",i3)')ni
              write(*,'("ak =",1pd15.4)')ak
            endif

!
!  Calculate 35 or 40% of the distance to the next ridge in order to determine
!   frequency interval over which the ridge will be fit.  If fnu(ik,nj)=1.0
!   it means that there was no guess for that ridge so the 5d0 is to weed out
!   that possibility.
!    
            if (fnu(ik,nj) < 5d0) then
               write(*,'("no guess for mode n = ",i3,",  ik = ",i5)')ni, ik  
               goto 1000
            endif
            if ((nj+1 <= nrdtot) .and. (fnu(ik,nj+1) > 5d0)) then
              frng = 0.40d0 * (fnu(ik,nj+1) - fnu(ik,nj)) / dnu
            else if (nj+1 > nrdtot) then
              frng = 0.35d0 * (fnu(ik,nj) - fnu(ik,nj-1)) / dnu
            endif
            if ((ni /= 0) .and. (fnu(ik,nj+1) < 5d0)) then
              frng = 0.35d0 * (fnu(ik,nj) - fnu(ik,nj-1)) / dnu
            endif
            istrtnu = nint(fnu(istk,nj) / dnu - frng)
            iendnu = nint(fnu(iendk,nj) / dnu + frng)
       
            if (istrtnu > iendnu) then
              frng = 0.35d0 * (fnu(ik,nj) - fnu(ik,nj-1)) / dnu
              istrtnu = nint(fnu(istk,nj) / dnu - frng)
              iendnu = nint(fnu(iendk,nj) / dnu + frng)
            endif
!
!if starting and stopping end values still inappropriate due to 
!extension of guess table, use asymptotic relationship to give
!it a try
!
            gsun=2.74d+8 !(Mm/microsec^2)
            if (istrtnu > iendnu .and. nj /= 0) then  
              fmode=dble(sqrt(gsun*ak)) / twopi             
              write(*,'("fmode = ", 1pd15.5)') fmode
              frng = 0.4d0*fmode*sqrt(4./3.)*dble(sqrt(real(nj)+1.5)-sqrt(real(nj)+0.5)) / dnu
              write(*,'("asymptotic frng = ", 1pd15.5)') frng
              istrtnu = nint(fnu(istk,nj) / dnu - frng)
              iendnu = nint(fnu(iendk,nj) / dnu + frng)
              write(*,'("bad strt and stop frequency bins")')
              write(*,'("recalculated using asymptotic formula",2(i6,2x))') &
                   &   istrtnu,iendnu
            endif
            if (iendnu > nnu) iendnu = nnu
            if (verbose == 1)then
             write(*,'(" istk, iendk ",2(i6,2x))'),istk,iendk
             write(*,'(" istrtnu,iendnu ",2(i6,2x))'),istrtnu,iendnu
             call flush(6)
            endif
            
!     
!  Set up starting guess for model parameters 
!     
 
            pj(1) = ux_guess
            pj(2) = uy_guess
            pj(3) = amp(ik,nj)
            pj(4) = bkg(ik,nj)
            pj(5) = width(ik,nj)
            pj(6) = 0d0            
            chng = 0
            if (verbose == 1)then
              write(*,'("guess parameters",6(1pe12.4,2x))')(pj(i),i=1,6)
            endif
            
!
!  Reset ioptions(1) for fitting each mode and fitflag
!
            ioptions(1) = maxeval
            fitflag = 0

            if (verbose == 1)then
              write(*,'("nthts nk nnu nrdtot fnu",4i6,1pd14.3)')nthts,nk,nnu, &
                         nrdtot, fnu(ik,nj)
            endif
            call dogleg(np, funk1, df1, pj, y1, doptions, ioptions, &
                       nthts, nk, nnu, pow, fnu, nrdtot ,verbose)

            if (verbose == 1)then
              write(*,'("after dogleg",4i6,1pd14.3)')nthts,nk,nnu, &
                         nrdtot, fnu(ik,nj)
            endif
!
!  ioptions(1) is reset by dogleg to number of iterations actually used
!
            neval = ioptions(1)

            if (verbose .eq. 1) then
              write(*,'(" parameters after call to dogleg = ", &
                 <np>(1pd13.4))') (pj(j),j=1,np)

              write(*,'(" function at minimum = ",1pd13.3)') y1
              write(*,'(" number of function evaluations  = ",i6)') &
                   neval
            endif
!
!  Only use fits that have converged before the maximum number of evaluations
!

            if (neval .lt. maxeval) then

!
!  If parameters are in the wrong ranges try once more using parameters from 
!    previous closest fit.  If still in wrong range go to next mode.
!      (not yet implemented) 

!  set boundaries on widths that varies with temporal resolution of spectrum being fit
              widmin = log(4.0*dnu)
              widmax = log(2.0*frng*0.85*dnu)  
             if (verbose == 1) then
               write(*,'("widmin = ",1pd13.4)') widmin
               write(*,'("widmax = ",1pd13.4)') widmax
             endif
             if((pj(5) .gt. widmax) .or. (pj(5) .lt. widmin)) then 
              if(verbose .eq. 1) then
               write(*,*)'width out of range'
               write(*,'(" log(width) = ",1pd13.3)') pj(5)
!              write(*,*)' trying with update from previous non-zero fit of ridge'
               write(*,*)' continue to next mode'
              endif
               chng = 1
               goto 1000
             endif
             if (pj(3) .lt. pj(4)) then
              if(verbose .eq. 1) then
               write(*,*)' amplitude is less than background'
               write(*,'(" log amplitude = ",1PD13.4)') pj(3)
               write(*,'(" log background = ",1PD13.4)') pj(4)
               write(*,*)' continue to next mode'
              endif
               chng = 1
               goto 1000
             endif
             if ((abs(pj(6)) .gt. 75d0) .and. (dk .lt. 0.1d0)) then 
              if(verbose .eq. 1) then
               write(*,'(" change from starting frequency too large")')
               write(*,'(" change in frequency = ",1PD13.4)') pj(6)
              endif
               chng = 1
               goto 1000
             endif
             if ((dk .ge. 0.1d0) .and. (70d0 < abs(pj(6))) .and. &
                    (abs(pj(6)) < 85d0)) then 
                !write out fit with caveat in fitflag
                chng = 0
                fitflag = 3
             else if ((dk .ge. 0.1d0) .and. (70d0 < abs(pj(6))) .and. &
                    (abs(pj(6)) > 85d0)) then
                chng = 1
                goto 1000
             endif
             rdchi = y1 / dble(npj-np)
!
! Calculate Hessian Error Matrix
!
              call df2(np,pj,flag,nthts,nk,nnu,pow,fnu,nrdtot,s2)
!
!  Invert Matrix.   The 'L' means that only the lower half of the final
!                   matrix holds the correct values
!                   
              uplo='L'
              info=0

              call dpotrf(uplo,np,s2,np,info)

              if (info .ne. 0) then
               if (verbose == 1)then
                 write(*,'(" Error in matrix factorization = ",i5)')info
               endif
                chng=1
                goto 1000
              endif

              call dpotri(uplo,np,s2,np,info)

              if (info .ne. 0)then
               if (verbose == 1)then
                write(*,'(" Error in matrix inversion = ",i5)')info
               endif
                chng=1
                goto 1000
              endif
              
              do i=1,np 
                do j=1,np
                   d_pj(j,i) = sqrt(dabs(s2(j,i)))
                enddo
              enddo
!
!  check to see if background term is reasonable by checking the 
!   size of the background errors but write out in any case
!
              if (d_pj(4,4) .gt. 3d0*abs(pj(4)))then
                if (verbose == 1)then
                 write(*,'(" background error is large")')
                 write(*,'(" background error = ", 1pd13.4)') d_pj(4,4)
                 write(*,'(" absolute background = ", 1pd13.4)') &
     &                 abs(pj(4))
                endif
                 chng = 0
                 fitflag = 2 + fitflag
              endif
!
!  put variables in correct form 
!
              nnj = nj - nmin    ! in case nmin /= 0
              if (chng .eq. 0) then 
                modecount = modecount + 1
                good(ik,nnj) = 1
                anu(ik,nnj) = fnu(ik,nj) + pj(6)
                ux(ik,nnj) = sngl(pj(1))
                uy(ik,nnj) = sngl(pj(2))
                ap(ik,nnj) = sngl(pj(3) + elnmult)
                bg(ik,nnj) = sngl(pj(4) + elnmult)
                fwhm(ik,nnj) = sngl(pj(5)) 
                delnu(ik,nnj) = sngl(pj(6))
                d_ux(ik,nnj) = sngl(d_pj(1,1))
                d_uy(ik,nnj) = sngl(d_pj(2,2))
                d_ap(ik,nnj) = sngl(d_pj(3,3))
                d_bg(ik,nnj) = sngl(d_pj(4,4))
                d_fwhm(ik,nnj) = sngl(d_pj(5,5))
                d_nu(ik,nnj) = sngl(d_pj(6,6))
                min_func(ik,nnj) = sngl(y1)
                chi(ik,nnj) = sngl(rdchi)
                n_iter(ik,nnj) = neval
                fit_chk(ik,nnj) = fitflag
                kbin(ik,nnj) = ik ! k in bin numbers
                nval(ik,nnj) = ni
                kval = kbin * sngl(dk) ! get k in inverse Mm.
                if (verbose .eq. 1) then 
                  write(*,fmt = 400) ni,ik,(pj(j),j = 1,np), &
     &              rdchi
                  write(*,500)ni,ik,((d_pj(j,i),j = 1,np),i = 1,np)
                  call flush(6)
                endif
              endif
            else
              if(verbose == 1)then
                write(*,*) 'Dogleg did not converge.'
              endif
               chng = 2
            endif            
 1000    enddo
!         close(16)
      enddo
      if (modecount == 0) then
        write(*,'(" no fits")')
        ierr = 350
        goto 9999 
      endif
!    if (verbose .eq. 1) then
!      write(*,'(" kval = ", <nk>(1pe12.4))')((kval(ik,nnj),ik=1,nk),nnj=1,nmax+1)
!      write(*,'(" mn_fn = ", <nk>(1pe12.4))')((min_func(ik,nnj),ik=1,nk),nnj=1,nmax+1)
!      write(*,'(" anu = ",1pe12.4)')(anu(ik,nnj)
!      write(*,'(" ux  dx  = ",4(1pe12.4))')pj(1),d_pj(1,1),ux(ik,nnj),d_ux(ik,nnj)
!      write(*,'(" uy  dy  = ",4(1pe12.4))')pj(2),d_pj(2,2),uy(ik,nnj),d_uy(ik,nnj)
!      write(*,'(" ap  d_ap  = ",4(1pe12.4))')pj(3),d_pj(3,3),ap(ik,nnj),d_ap(ik,nnj)
!      write(*,'(" bg  d_bg  = ",4(1pe12.4))')pj(4),d_pj(4,4),bg(ik,nnj),d_bg(ik,nnj)
!      write(*,'(" g  d_g  = ",4(1pe12.4))')pj(5),d_pj(5,5),fwhm(ik,nnj),d_fwhm(ik,nnj)
!      write(*,'(" delnu  d_nu  = ",4(1pe12.4))')pj(6),d_pj(6,6),delnu(ik,nnj),d_nu(ik,nnj)
!      write(*,'(" nev  n_it y1 mf = ",2i6,2(1pe12.4))')neval,n_iter(ik,nnj),y1,min_func(ik,nnj)
!      write(*,'(" rdchi  chi ff fc = ",2(1pe12.4),2i6)')rdchi,chi(ik,nnj),fitflag,fit_chk(ik,nnj)
!    endif 
9999  continue
 300  format(i2,1x,f8.3,1x,f8.5,1x,f9.3,1pe12.4, &
     &       <np>(1pe12.4,1pe12.4),i5,i6,1pe12.4,1pe12.3,2x,i1) 
 400  format(i3,2x,i3,2x,<np>(1pd12.4,1x),d12.4)
 500  format(i2,2x,i3,/,/,<np>(<np>(1PD13.4,1x),/))
      return
      end subroutine
    end module
