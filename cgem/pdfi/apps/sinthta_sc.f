      subroutine sinthta_sc(thmin,thmax,m,sinth,sinth_hlf,sinth_gh,
     1 sinth_hlf_gh)
c
c+
c - - Purpose: compute arrays of sin(theta) for the edge locations and at the
c - - half-cell locations.  Compute arrays with and without ghost zones.
c - - Note that this version uses centered grid, not staggered grid.
c
c     Usage: call(sinthta_sc,thmin,thmax,m,sinth,sinth_hlf,sinth_gh,
c    1 sinth_hlf_gh)
c
c - - input: thmin, thmax - minimum and maximum values (in radians) of
c - - co-latitude theta at edges
c
c - - input: m (number of cell interiors in the theta (colatitude) direction
c - - IMPORTANT - m is one less than the array size for the 1st (theta)
c - - dimension of the data arrays!
c
c - - output: sinth, an m+1 length array of cell edge values of sin(theta)
c
c - - output: sinth_gh, an m+3 length array of cell edge values of sin(theta)
c - - including one ghost-zone at each end
c
c - - output:  sinth_hlf, an m length array of cell center values of sin(theta)
c
c - - output:  sinth_hlf_gh, an m+2 length array of cell center values of
c - - sin(theta) including one ghost-zone at each end
c - - Caution!!!! The use of m here is incompatible with the use of m in
c - - PDFI_SS!  Be careful!
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
      integer :: m
      real*8 thmin,thmax
      real*8 sinth(m+1),sinth_gh(m+3),sinth_hlf(m),sinth_hlf_gh(m+2)
c
      integer :: i
      real*8 :: theta(m+1),theta_gh(m+3),theta_hlf(m),
     1 theta_hlf_gh(m+2)
c
      real*8 :: deltheta,deltheta_hlf
      deltheta=(thmax-thmin)/m
      deltheta_hlf=0.5d0*deltheta
      do i=1,m
         theta_hlf(i)=thmin+0.5d0*deltheta_hlf+(i-1)*deltheta
         theta_hlf_gh(i+1)=theta_hlf(i)
         theta(i)=thmin+(i-1)*deltheta
         theta_gh(i+1)=theta(i)
      enddo
      theta_hlf_gh(1)=theta_hlf(1)-deltheta
      theta(m+1)=thmax
      theta_gh(m+2)=thmax
      theta_gh(1)=thmin-deltheta
      theta_gh(m+3)=thmax+deltheta
      theta_hlf_gh(m+2)=thmax+0.5d0*deltheta
      sinth(1:m+1)=sin(theta(1:m+1))
      sinth_gh(1:m+3)=sin(theta_gh(1:m+3))
      sinth_hlf(1:m)=sin(theta_hlf(1:m))
      sinth_hlf_gh(1:m+2)=sin(theta_hlf_gh(1:m+2))
c
      return
      end
