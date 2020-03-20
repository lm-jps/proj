      subroutine car2sph_ss(m,n,delth,delx,dely,rsph,a,b,c,d)
c
c+
c - - Purpose:  Compute the radius rsph, and colatitude, longitude
c               domain limits a,b,c,d from m,n,delth,delx,dely.
c               This allows one to perform calculations in Cartesian
c               Coordinates by using the output values of rsph,a,b,c,d
c               in other PDFI_SS subroutines.  The values of a,b will
c               be such that spherical geometry factors will be very near
c               unity. See ../doc/Cartesian-Solutions-with-PDFI_SS.txt for
c               further details.
c
c - - Usage:    call car2sph_ss(m,n,delth,delx,dely,rsph,a,b,c,d)
c
c - - Input:    m,n - integer number of cell interiors in the colat,lon 
c               directions, respectively.  n is also the number of cell
c               interiors in x direction, m is number of cell interiors in y.
c
c - - Input:    delth - real*8 value of colatitude subtended by m*dely.
c               Choosing delth=1d-4 has worked well in the past.  If delth
c               is equal to 0 on input, a value of 1d-4 will be used internal
c               to the subroutine.
c
c - - Input:    delx - real*8 value of thickness in x [km] for one cell.
c
c - - Input:    dely - real*8 value of thickness in y [km] for one cell.
c
c - - Output:   rsph - real*8 value of radius of sphere [km] onto which the
c               Cartesian patch will be placed.
c
c - - Output:   a,b - real*8 values of the colatitude range for the problem 
c               domain, a < b [radians]
c
c - - Output:   c,d - real*8 values of the longitude range for the problem 
c               domain, c < d [radians]
c
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
      implicit none
c
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: delth,delx,dely
c
c - - output variable declarations:
c
      real*8 :: rsph,a,b,c,d
c
c - - local variable declarations:
c
      real*8 :: bma
      real*8 :: halfpi,dum
c
c - - declaration of function from fishpack (fftpack) used to compute PI:
c
      real*8 :: pimach
c
      halfpi=0.5d0*pimach(dum)
c
      if(delth .le. 0.d0) then
        bma=1d-4
      else
        bma=delth
      endif
c
c - - compute a,b:
c
      a=halfpi - 0.5d0*bma
      b=halfpi + 0.5d0*bma
c
c - - compute radius of sphere rsph:
c
      rsph=dely*dble(m)/bma
c
c - - compute c,d (assume c=0)
c
      c=0.d0
      d=delx*dble(n)/rsph
c
c - - we're done
c
      return
      end
