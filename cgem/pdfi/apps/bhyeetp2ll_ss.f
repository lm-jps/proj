      subroutine bhyeetp2ll_ss(m,n,btte,bppe,blon,blat)
c
c+
c - - Purpose: To transpose B_h data arrays from theta,phi to lon,lat order
c              and flip sign to get B_lat.  Array values are assumed at
c              staggered Yee grid locations.
c
c - - Usage:   call bhyeetp2ll_ss(m,n,btte,bppe,blon,blat)
c
c - - Input:   m,n - integer number of cell centers in the theta (lat), 
c              and phi (lon) directions, respectively.
c
c - - Input:   btte (m+1,n),bppe(m,n+1) - real*8 arrays of the co-latitudinal 
c              and azimuthal components of the magnetic field evaluated at
c              TE and PE locations (theta and phi edges, resp.) [G]
c
c - - Output:  blon(n+1,m),blat(n,m+1) - real*8 arrays of longitudinal and
c              latitudinal components of magnetic field, stored in lon,lat
c              index order. [G]
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
c - - input variables:
      integer :: m,n
      real*8 :: btte(m+1,n),bppe(m,n+1)
c - - output variables:
      real*8 :: blat(n,m+1),blon(n+1,m)
c - - local variables:
      integer :: i,j
c
c - - note that colat, lat unit vectors have opposite sign, so must change 
c     sign of latitude component from colat component.
c
      do i=1,m+1
         do j=1,n
            blat(j,i)=-btte(m+2-i,j)
         enddo
      enddo
c
      do i=1,m
         do j=1,n+1
            blon(j,i)=bppe(m+1-i,j)
         enddo
      enddo
c
      return
      end
