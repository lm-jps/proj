      subroutine bryeetp2ll_ss(m,n,brtpce,brllce)
c
c+
c - - Purpose: To transpose B_r data array from theta,phi to lon,lat order.
c - - Usage:  call bryeetp2ll_ss(m,n,brtpce,brllce)
c - - Input:  m,n - number of cell centers in the theta (lat), and phi (lon)
c             directions, respectively.
c - - Input:  brtpce(m,n) - real*8 array of the radial component of the
c             magnetic field evaluated at CE locations in theta,phi order
c             (cell-centers). [G]
c - - Output: brllce(n,m) - real*8 array of radial magnetic field component
c             stored in lon,lat index order. [G]
c-
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
      implicit none
c
      integer :: m,n
      real*8 :: brtpce(m,n)
c
      real*8 :: brllce(n,m)
c
      integer :: i,j
      do i=1,m
         do j=1,n
            brllce(j,i)=brtpce(m+1-i,j)
         enddo
      enddo
c
      return
      end
