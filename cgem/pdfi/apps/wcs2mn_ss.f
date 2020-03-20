      subroutine wcs2mn_ss(crpix1,crpix2,m,n)
c
c+
c - -    Purpose:  Compute the number of cells in latitude and longitude, m,n
c        from CRPIX1 and CRPIX2 WCS keywords for the COE grid.
c
c        NOTE:  This subroutine assumes the COE grid is in longitude-latitude
c
c - -    Usage:  call wcs2mn_ss(crpix1,crpix2,m,n)
c
c - -    Input: crpix1 - real*8 value of the reference pixel in longitude
c        direction, for the COE grid
c
c - -    Input: crpix2 - real*8 value of the reference pixel in latitude
c        direction, for the COE grid
c
c - -    Output:  m,n - integer number of cell interiors in the lat,lon 
c        directions respectively
c
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
      real*8 :: crpix1,crpix2
c
c - - output variable declarations:
c
      integer :: m,n
c
c - - local variable declarations:
c
      integer :: mtmp,ntmp
      real*8 :: epsilon
c
      epsilon=1.d-6
c
      ntmp=int(2.*crpix1+epsilon) - 2
      mtmp=int(2.*crpix2+epsilon) - 2
c
      n=ntmp
      m=mtmp
c
      return
      end
