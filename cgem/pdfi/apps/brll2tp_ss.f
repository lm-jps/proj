      subroutine brll2tp_ss(m,n,brllcoe,brtpcoe)
c
c+
c - - Purpose: To transpose B_r data array from lon,lat to theta,phi order.
c             Can also be used to transpose other scalar COE arrays.
c
c - - Usage:  call brll2tp_ss(m,n,brllcoe,brtpcoe)
c
c - - Input:  m,n - number of cell centers in the theta (lat), and phi (lon)
c             directions, respectively.
c
c - - Input:  brllcoe(n+1,m+1) - real*8 array of the radial component of the
c             magnetic field evaluated at COE locations
c             (corners plus exterior corners on boundary). [G]
c
c - - Output: brtpcoe(m+1,n+1) - real*8 array of radial magnetic field component
c             at COE locations, stored in theta,phi index order. [G]
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
      integer :: m,n
      real*8 :: brllcoe(n+1,m+1)
c
      real*8 :: brtpcoe(m+1,n+1)
c
      integer :: i,j
      do i=1,m+1
         do j=1,n+1
            brtpcoe(i,j)=brllcoe(j,m+2-i)
         enddo
      enddo
c
      return
      end
