      FUNCTION causal_init_ss(f, fsize, z)
c
c+
c - - Purpose: real*8 Function needed in bspline_ss.f
c - - Usage: a=causal_init_ss(f,fsize,z)
c - - Input:  f(fsize) - real*8 array
c - - Input:  fsize - integer value of the dimension of f
c - - Input:  z - real*8 argument
c - - Note:  This function is *only* used by subroutine bspline_ss.
c-
c - -  PDFI_SS electric field inversion software
c - -  http://cgem.ssl.berkeley.edu/~fisher/public/software/PDFI_SS
c - -  Copyright (C) 2015-2019 Regents of the University of California
c 
c - -  This software is based on the concepts described in Kazachenko et al.
c - -  (2014, ApJ 795, 17).  A detailed description of the software is in
c - -  Fisher et al. (2019, arXiv:1912.08301 ).
c - -  If you use the software in a scientific 
c - -  publication, the authors would appreciate a citation to these papers 
c - -  and any future papers describing updates to the methods.
c
c - -  This is free software; you can redistribute it and/or
c - -  modify it under the terms of the GNU Lesser General Public
c - -  License as published by the Free Software Foundation,
c - -  version 2.1 of the License.
c
c - -  This software is distributed in the hope that it will be useful,
c - -  but WITHOUT ANY WARRANTY; without even the implied warranty of
c - -  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c - -  See the GNU Lesser General Public License for more details.
c
c - -  To view the GNU Lesser General Public License visit
c - -  http://www.gnu.org/copyleft/lesser.html
c - -  or write to the Free Software Foundation, Inc.,
c - -  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c   Origins:  This function was written in Fortran by Dave Bercik.
c   It is based on software in C, published by Philippe Thevenaz, 
c   http://bigwww.epfl.ch/thevenaz/interpolation/  
c   originally described in IEEE TRANSACTIONS ON MEDICAL IMAGING, 
c   VOL. 19, NO. 7, JULY 2000.  It was subsequently edited for compatibility 
c   with the PDFI_SS library by George Fisher. 
c  
      IMPLICIT NONE
      INTEGER, PARAMETER :: r8=KIND(1.d0)
      INTEGER, PARAMETER :: i4=KIND(1)

c    Dummy variable declarations
      INTEGER(KIND=i4) :: fsize
      REAL(KIND=r8) :: f(fsize)
      REAL(KIND=r8) :: z
c     Local variable declarations
      INTEGER(KIND=i4) :: n, horizon
      REAL(KIND=r8) :: causal_init_ss, zn, iz, z2n
c
      horizon = INT(CEILING(LOG(EPSILON(1.0_r8))/LOG(ABS(z))), KIND=i4)
      IF (horizon < fsize) THEN  
c accelerated loop
        zn = z
        causal_init_ss = f(1)
        DO n=2,horizon
          causal_init_ss = causal_init_ss + zn*f(n)
          zn = zn*z
        END DO
      ELSE  
c full loop
        zn = z
        iz = 1.0_r8/z
        z2n = z**(fsize-1)
        causal_init_ss = f(1) + z2n*f(fsize)
        z2n = z2n*z2n*iz
        DO n=2,fsize-1
           causal_init_ss = causal_init_ss + (zn + z2n)*f(n)
           zn = zn*z
           z2n = z2n*iz
        END DO
        causal_init_ss = causal_init_ss/(1.0_r8 - zn*zn)
      END IF
c
      END FUNCTION causal_init_ss
