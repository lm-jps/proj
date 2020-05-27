      FUNCTION anticausal_init_ss(f, fsize, z)
c
c+
c - - Purpose: real*8 Function needed by subroutine bspline_ss.f
c - - Usage: a=anticausal_init_ss(f,fsize,z)
c - - Input:  f(fsize) - real*8 array
c - - Input:  fsize - integer value of the dimension of f
c - - Input:  z - real*8 argument
c - - Note:  This function is *only* used inside of subroutine bspline_ss
c-
c   Origins:  This function was written in Fortran by Dave Bercik.
c   It is based on software in C, published by Philippe Thevenaz, 
c   http://bigwww.epfl.ch/thevenaz/interpolation/  
c   originally described in IEEE TRANSACTIONS ON MEDICAL IMAGING, 
c   VOL. 19, NO. 7, JULY 2000.  It was subsequently edited for compatibility 
c   with the PDFI_SS library by George Fisher. 
c  
c   PDFI_SS Electric Field Inversion Software
c   http://cgem.ssl.berkeley.edu/cgi-bin/cgem/PDFI_SS/index
c   Copyright (C) 2015-2018 University of California
c  
c   This software is based on the concepts described in Kazachenko et al. 
c   (2014, ApJ 795, 17).  It also extends those techniques to 
c   spherical coordinates, and uses a staggered, rather than a centered grid.
c   If you use the software in a scientific publication, 
c   the authors would appreciate a citation to this paper and any future papers 
c   describing updates to the methods.
c  
c   This is free software; you can redistribute it and/or
c   modify it under the terms of the GNU Lesser General Public
c   License as published by the Free Software Foundation;
c   either version 2.1 of the License, or (at your option) any later version.
c  
c   This software is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c   See the GNU Lesser General Public License for more details.
c  
c   To view the GNU Lesser General Public License visit
c   http://www.gnu.org/copyleft/lesser.html
c   or write to the Free Software Foundation, Inc.,
c   59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: r8=KIND(1.d0)
      INTEGER, PARAMETER :: i4=KIND(1)
c     Dummy variable declarations
      INTEGER(kind=i4) :: fsize
      REAL(KIND=r8) :: f(fsize)
      REAL(KIND=r8) :: z
c     Local variable declarations
      REAL(KIND=r8) :: anticausal_init_ss
c
      anticausal_init_ss = z/(z*z - 1.0_r8)*(z*f(fsize-1) + f(fsize))
      END FUNCTION anticausal_init_ss
