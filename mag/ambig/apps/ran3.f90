!***********************************************************************************************************************************
function ran3(seed)
!===================================================================================================================================
!     "Minimal standard" pseudo-random number generator of Park and
!     Miller.  Returns a uniform random deviate r s.t. 0 < r < 1.0.
!     Set seed to any non-zero integer value to initialize a sequence,
!     then do not change seed between calls for successive deviates
!     in the sequence.
!
!     References:
!        Park, S. and Miller, K., "Random Number Generators: Good Ones
!           are Hard to Find", Comm. ACM 31, 1192-1201 (Oct. 1988)
!        Park, S. and Miller, K., in "Remarks on Choosing and Imple-
!           menting Random Number Generators", Comm. ACM 36 No. 7,
!           105-110 (July 1993)
!
! This is ran0 from the genetic algorithm code PIKAIA developed by Paul Charbonneau & Barry Knapp
!===================================================================================================================================
! *** Declaration section ***
!
   implicit none
!
!     Input/Output:
   integer :: seed
!
!     Output:
   real :: ran3
!
!     Constants:
   integer,parameter :: A=48271,M=2147483647,Q=44488,R=3399
   real,parameter :: SCALE=1./M,EPS=1.2e-7,RNMX=1.-EPS
!
!     Local:
   integer :: j
!
! *** Executable section ***
!
   j = seed/Q
   seed = A*(seed-j*Q)-R*j
   if (seed .lt. 0) seed = seed+M
   ran3 = min(seed*SCALE,RNMX)

end function ran3
!***********************************************************************************************************************************
