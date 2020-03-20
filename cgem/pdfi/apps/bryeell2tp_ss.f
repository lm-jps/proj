      subroutine bryeell2tp_ss(m,n,brllce,brtpce)
c
c+
c - - Purpose: To transpose B_r data array from lon,lat to theta,phi order.
c
c - - Usage:  call bryeell2tp_ss(m,n,brllce,brtpce)
c
c - - Input:  m,n - number of cell centers in the theta (lat), and phi (lon)
c             directions, respectively.
c
c - - Input:  brllce(n,m) - real*8 array of radial magnetic field component
c             at cell centers, stored in lon,lat index order. [G]
c
c - - Output: brtpce(m,n) - real*8 array of the radial component of the
c             magnetic field evaluated at CE locations in theta,phi order
c             (cell-centers). [G]
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
      real*8 :: brtpce(m,n)
c
      real*8 :: brllce(n,m)
c
      integer :: i,j
      do i=1,m
         do j=1,n
            brtpce(m+1-i,j)=brllce(j,i)
         enddo
      enddo
c
      return
      end
