      subroutine wcs2abcd_ss(m,n,crval1,crval2,cdelt1,cdelt2,a,b,c,d)
c
c+
c - -    Purpose:  Compute the domain limits a,b (colatitude) and c,d 
c        (longitude), in radians, from the crval1,crval2,cdelt1, and cdelt2
c        WCS keywords, along with m,n.  Units of crval,cdelt assumed in degrees.
c
c        NOTE:  This subroutine assumes the arrays are in longitude-latitude
c        order as e.g. on output from pdfi_wrapper4jsoc_ss.
c
c - -    Usage:  call wcs2abcd_ss(m,n,crval1,crval2,cdelt1,cdelt2,a,b,c,d)
c
c - -    Input:  m,n - integer number of cell interiors in the lat,lon 
c        directions respectively
c
c - -    Input: crval1 - real*8 longitude value of the reference pixel
c        [degrees]
c
c - -    Input: crval2 - real*8 latitude value of the reference pixel
c        [degrees]
c
c - -    Input:  cdelt1 - real*8 - number of degrees per pixel in longitude
c        [degrees/pixel]
c
c - -    Input:  cdelt2 - real*8 - number of degrees per pixel in latitude
c        [degrees/pixel]
c
c - -    Output:  a,b, - real*8 values of colatitude range for the problem 
c        domain, a < b [radians]
c
c - -    Output:  c,d - real*8 values of the longitude range for the problem, 
c        c < d [radians]
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
      real*8 :: crval1,crval2
      real*8 :: cdelt1,cdelt2
c
c - - output variable declarations:
c
      real*8 :: a,b,c,d
c
c - - local variable declarations:
c
      real*8 :: minlon,maxlon,minlat,maxlat,deg2rad,pi,dum
      real*8 :: mincolat,maxcolat
c
c - - declaration of pimach function from FISHPACK/FFTPACK
c
      real*8 :: pimach
c
c - - get value of pi from pimach function:
c
      pi=pimach(dum)
c
      deg2rad=pi/180.d0
c
      minlon=crval1-0.5*dble(n)*cdelt1
      maxlon=crval1+0.5*dble(n)*cdelt1
      minlat=crval2-0.5*dble(m)*cdelt2
      maxlat=crval2+0.5*dble(m)*cdelt2
c
      mincolat=90.d0 - maxlat
      maxcolat=90.d0 - minlat
c
      a=deg2rad*mincolat
      b=deg2rad*maxcolat
      c=deg2rad*minlon
      d=deg2rad*maxlon      
c
      return
      end
