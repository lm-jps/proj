      subroutine abcd2wcs_ss(m,n,a,b,c,d,crpix1,crpix2,crval1,crval2,
     1           cdelt1,cdelt2)
c
c+
c - -    Purpose:  Compute the WCS keywords crpix1, crpix2, crval1,
c        crval2, cdelt1, and cdelt2 from a,b,c,d,m,n.
c        crpix1 and crpix2 are arrays of length 6, and compute the reference
c        pixels for the COE, CO, CE, CEG, TE, and PE grids, in that order.
c
c        NOTE: - This subroutine assumes that the underlying arrays are
c        in longitude-latitude order, as e.g. they are on output from
c        pdfi_wrapper4jsoc_ss.
c
c - -    Usage:  call abcd2wcs_ss(m,n,a,b,c,d,crpix1,crpix2,crval1,crval2,
c             cdelt1,cdelt2)
c
c - -    Input:  m,n - integer number of cell interiors in the lat,lon 
c        directions respectively
c
c - -    Input:  a,b, - real*8 values of colatitude range for the problem 
c        domain, 0 < a < b < pi [radians]
c
c - -    Input:  c,d - real*8 values of the longitude range for the problem, 
c        0 < c < d < 2*pi [radians]
c
c - -    Output: crpix1(6) - real*8 array of reference longitude pixel 
c        locations for all 6 (COE, CO, CE, CEG, TE, PE) grids [pixels]
c
c - -    Output: crpix2(6) - real*8 array of reference latitude pixel 
c        locations for all 6 (COE, CO, CE, CEG, TE, PE) grids [pixels]
c
c - -    Output: crval1 - real*8 longitude value of the reference pixel
c        [degrees]
c
c - -    Output: crval2 - real*8 latitude value of the reference grid
c        [degrees]
c
c - -    Output:  cdelt1 - real*8 - number of degrees per pixel in longitude
c        [degrees/pixel]
c
c - -    Output:  cdelt2 - real*8 - number of degrees per pixel in latitude
c        [degrees/pixel]
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
      real*8 :: a,b,c,d
c
c - - output variable declarations:
c
      real*8 :: crpix1(6),crpix2(6)
      real*8 :: crval1,crval2
      real*8 :: cdelt1,cdelt2
c
c - - local variable declarations:
c
      real*8 :: dtheta,dphi,pi,dum,minlon,maxlon,minlat,maxlat,rad2deg
      integer :: cols(6),rows(6)
c
c - - declaration of pimach function from FISHPACK/FFTPACK
c
      real*8 :: pimach
c
c - - get value of pi from pimach function:
c
      pi=pimach(dum)
c
      rad2deg=180.d0/pi
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
      cdelt1=dphi*rad2deg
      cdelt2=dtheta*rad2deg
c
      minlon=rad2deg*c
      maxlon=rad2deg*d
      minlat=90.d0 - b*rad2deg
      maxlat=90.d0 - a*rad2deg
      crval1=0.5d0*(minlon+maxlon)
      crval2=0.5d0*(minlat+maxlat)
c
c - - 1st dimension of COE, CO, CE, CEG, TE, PE grids (lon-lat order):
c
      cols(:) = [n+1,n-1,n,n+2,n,n+1]
c
c - - 2nd dimension of COE, CO, CE, CEG, TE, PE grids (lon-lat order):
c
      rows(:) = [m+1,m-1,m,m+2,m+1,m]
c
      crpix1(:) = (1.d0 + dble(cols(:)) )*0.5d0
      crpix2(:) = (1.d0 + dble(rows(:)) )*0.5d0
c
      return
      end
