      subroutine coell_ss(m,n,a,b,c,d,lon,lat)
c
c+
c - - Purpose:  Compute the longitude and latitude locations of the COE grid
c
c - - Usage:  call coetp_ss(m,n,a,b,c,d,theta,phi)
c
c - - Input:  m,n - number of cell interiors in the colat,lon directions
c - - Input:  a,b,c,d - a,b are colatitude range for the problem domain, a < b
c - - Input:  c,d - c,d are the longitude range for the problem, c < d
c - - Output: lon(n+1,m+1) - the longitude values of the COE grid
c - - Output: lat(n+1,m+1) - the latitude values of the COE grid
c
c-
c   PDFI_SS Electric Field Inversion Software
c   http://cgem.ssl.berkeley.edu/cgi-bin/cgem/PDFI_SS/index
c   Copyright (C) 2015,2016 University of California
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
      implicit none
c
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: a,b,c,d
c
c - - output variable declarations:
c
      real*8 :: lon(n+1,m+1)
      real*8 :: lat(n+1,m+1)
c
c - - local variable declarations:
c
      real*8 :: dtheta,dphi
      real*8 :: halfpi,dum
      integer :: i,j
c
c - - declaration of function from fishpack (fftpack) used to compute PI:
c
      real*8 :: pimach
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
      halfpi=0.5d0*pimach(dum)
c
      do i=1,m+1
         lat(1:n+1,m+2-i)=halfpi - (a+dble(i-1)*dtheta)
      end do
c
      do j=1,n+1
         lon(j,1:m+1)=c+dble(j-1)*dphi
      end do
c
      return
      end
