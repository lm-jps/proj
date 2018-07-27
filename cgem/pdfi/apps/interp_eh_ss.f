      subroutine interp_eh_ss(m,n,et,ep,etc,epc)
c+
c - - Purpose: Interpolate et and ep from edges to corners.  Designed so that
c              centered grid gradients and divergences computed at corners
c              are consistent with averages of adjacent staggered grid 
c              gradients and divergences
c
c - - Usage:   call interp_eh_ss(m,n,et,ep,etc,epc)
c
c - - Input:   m,n are number of cell centers in theta, phi directions, resp.
c
c - - Input:   et(m,n+1) - real*8 array of theta component of electric field
c              at PE grid locations. [G km/sec or V/cm]
c
c - - Input:   ep(m+1,n) - real*8 array of phi component of electric field
c              at TE grid locations. [G km/sec or V/cm]
c
c - - Output:  etc(m-1,n-1) - real*8 array of theta component of electric field
c              located on interior corners (CO grid) [G km/sec or V/cm]
c
c - - Output:  epc(m-1,n-1) - real*8 array of phi component of electric field
c              located on interior corners (CO grid) [G km/sec or V/cm]
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
c - - input variables:
c
      integer :: m,n
      real*8 :: et(m,n+1),ep(m+1,n)
c
c - - output variables:
c
      real*8 :: etc(m-1,n-1),epc(m-1,n-1)
c
c - - local variables:
c
      integer :: i,j
c
      do i=1,m-1
         do j=1,n-1
            etc(i,j)=0.5d0*(et(i+1,j+1)+et(i,j+1))
            epc(i,j)=0.5d0*(ep(i+1,j+1)+ep(i+1,j))
         enddo
      enddo
c
      return
      end
