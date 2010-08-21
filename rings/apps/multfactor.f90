!
!  Function which calculates a multiplication factor for 1/pow so as to reduce
!    roundoff errors during the fitting and make the initial power estimations 
!    closer to the parameterized power
!
!  Variables passed from module ring_init-------------------------------------- 
!
!    nthts,nk,nnu    sizes of subsampled power array in theta, wavenumber, 
!                     and frequency respectively
!
!    inumn, inumx   min and max of frequency bins to use in determining log(multfactor)
!
!    verbose        = 1, write info for debugging
!
!  Variables passed from module fits--------------------------------------
!
!    pow            filtered, unwrapped power (not in log space) Real*4  
!                                dimensions (nthts,nk,nnu)    
!
!  Input Variables--------------------------------------------------------
!
!    kbstrt, kbend  starting and stopping values of k to be used in fitting each ridge
!                      dimension (nrdg)
!
!  Housekeeping Variables---------------------------------------------------------
!
!    itht, ik, inu, i    step through azimuth, wavenumber, frequency, and ridges
!
   module factor
      integer, parameter :: d2=8
      contains 
!       REAL(kind=d2) FUNCTION MULTFACTOR(powfilt,nthts,nk,nnu,inumn,inumx,kbstrt,kbend,kbmin, &
!                                       &  kbmax,verbose,ierr)
       REAL(kind=d2) FUNCTION MULTFACTOR(powfilt,nthts,nk,nnu,inumn,inumx,kbmin,kbmax,verbose,ierr)
       IMPLICIT NONE
       integer, parameter :: double=8
       integer, parameter :: single=4
       INTEGER :: inu,ik,itht,kbmin,kbmax,i,ierr,ier,verbose
       INTEGER :: nthts,nk,nnu,inumn,inumx,kbminmult,kbmaxmult
       external dlamch
       real(kind=double) :: powmax, safemin, safemax, dlamch
       real(kind=single), dimension(nthts,nk,nnu) :: powfilt
!       integer, dimension(:) :: kbstrt, kbend

!  Now kbmin and kbmax are read in from c program 
!  ignore -Determine kbmin and kbmax from min and max of kbstrt and kbend respectively
!
!       kbmin=minval(kbstrt)
!       kbmax=maxval(kbend)
       kbminmult = kbmin + 6
       kbmaxmult = kbmax - 2
       if (verbose == 1) then  
          WRITE(*,'(" k bin range for determining multiplication &
               &factor",i4," - ",i4)') kbminmult, kbmaxmult
          WRITE(*,'(" nu bin range for determining multiplication &
               &factor",i4," - ",i4)') inumn, inumx
       endif
!
!  Determine multiplication factor in log space - powmax and powmin are for log(power)
!
       powmax = maxval(log(powfilt(:,kbminmult:kbmaxmult,inumn:inumx)))
       if (verbose == 1)then
         WRITE(*,'(" max logarithmic power = ",1pe13.4)') powmax
       endif
!
!  Keep max log power at about 9.0 
!
!  CHECK FOR OVERFLOW OR UNDERFLOW POSSIBILITIES HERE
!
!     Determine underflow and overflow limits 
!
       safemin = log(dlamch('s'))
       safemax = log(1.0d0/dexp(safemin))
!       if (verbose .eq. 1) then
!         write(*,'("safemin = ",1pd13.4)')safemin
!         write(*,'("safemax = ",1pd13.4)')safemax
!       endif

       if ((powmax < safemin) .or. (powmax > safemax)) then
         write(*,'(" Power will cause overflow or underflow ")')
         write(*,'(" Should not be a problem unless data is &
                   &wonky ")')
         write(*,'(" Problem determining multiplication factor")')
         ierr=32
       else
!         write(*,'("multfactor = ",1pe15.4)')exp(powmax)
         multfactor = dble(exp(powmax - 9.0))
        if (verbose == 1)then
         write(*,'("multfactor = ",1pd15.4)')multfactor
        endif
       endif
       RETURN
       END FUNCTION MULTFACTOR
   end module factor
