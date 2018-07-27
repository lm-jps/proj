      subroutine bhll2tp_ss(m,n,bloncoe,blatcoe,btcoe,bpcoe)
c
c+
c - - Purpose: To transpose B_h data arrays from lon,lat to theta,phi order
c              and flip sign to get B_theta.  This same subroutine can be
c              used to convert the components of v_lon and v_lat returned from
c              FLCT into components of v_theta and v_phi.  COE grid locations.
c
c - - Usage:   call bhll2tp_ss(m,n,bloncoe,blatcoe,btcoe,bpcoe)
c
c - - Input:   m,n - number of cell centers in the theta (lat), and phi (lon)
c              directions, respectively.
c
c - - Input:   bloncoe(n+1,m+1),blatcoe(n+1,m+1) - real*8 arrays of the 
c              longitudinal and latitudinal components of the magnetic field 
c              evaluated at COE locations (corners plus exterior corners on 
c              boundary) [G].
c
c - - Output:  btcoe(m+1,n+1),bpcoe(m+1,n+1) - real*8 arrays of colatitudinal 
c              and azimuthal components of magnetic field, stored in theta,phi
c              index order, COE grid locations. [G]
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
c - - input variables:
      integer :: m,n
      real*8 :: bloncoe(n+1,m+1),blatcoe(n+1,m+1)
c - - output variables:
      real*8 :: btcoe(m+1,n+1),bpcoe(m+1,n+1)
c - - local variables:
      integer :: i,j
c
c - - Note lat, colat unit vectors have opposite direction, so must change
c - - sign of theta component from lat component:
c
      do i=1,m+1
         do j=1,n+1
            bpcoe(i,j)=bloncoe(j,m+2-i)
            btcoe(i,j)=-blatcoe(j,m+2-i)
         enddo
      enddo
c
      return
      end
