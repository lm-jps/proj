      subroutine cell_ss(m,n,a,b,c,d,lon,lat)
c
c+
c - - Purpose: Compute the longitude and latitude locations of the CE grid
c
c - - Usage:  call cell_ss(m,n,a,b,c,d,lon,lat)
c
c - - Input:  m,n - integer value of number of cell interiors in the 
c             colat,lon directions, respectively
c
c - - Input:  a,b - real*8 values of colatitude range for the problem domain, 
c             a < b [radians]
c
c - - Input:  c,d - real*8 values of the longitude range for the problem, 
c             c < d [radians]
c
c - - Output: lon(n,m) - real*8 array of longitude values of the CE grid
c             [radians]
c
c - - Output: lat(n,m) - real*8 array of latitude values of the CE grid
c             [radians]
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
      real*8 :: lon(n,m)
      real*8 :: lat(n,m)
c
c - - local variable declarations:
c
      real*8 :: dtheta,dphi,hlfdth,hlfdph
      real*8 :: halfpi,dum
      integer :: i,j
c
c - - declaration of function from fishpack (fftpack) used to compute PI:
c
      real*8 :: pimach
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
      hlfdth=0.5d0*dtheta
      hlfdph=0.5d0*dphi
      halfpi=0.5d0*pimach(dum)
c
      do i=1,m
         lat(1:n,m+1-i)=halfpi - (a+hlfdth+dble(i-1)*dtheta)
      end do
c
      do j=1,n
         lon(j,1:m)=c+hlfdph+dble(j-1)*dphi
      end do
c
      return
      end
