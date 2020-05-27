      subroutine interp_ehcoe_ss(m,n,etcoe,epcoe,et,ep)
c+
c - -   Purpose: To interpolate a 2-d electric field vector from COE grid 
c                to its proper locations on the TE, PE.
c                Interpolation is done using straightforward linear
c                averages, motivated by simplicity and noise reduction.
c
c - -   Usage:   call interp_ehcoe_ss(m,n,etcoe,epcoe,et,ep)
c
c - -    Input:  m,n - integer values of the number of cell centers in colat
c                and lon directions, respectively.
c
c - -    Input:  etcoe(m+1,n+1) - real*8 array of theta component of E located
c                on the COE grid (corners, including boundary corners) 
c                [G km/sec or V/cm]
c
c - -    Input:  epcoe(m+1,n+1) - real*8 array of phi component of E located
c                on the COE grid [G km/sec or V/cm]
c
c - -   Output:  et(m,n+1) - real*8 array of theta component of E located
c                at theta edge locations (PE grid) [G km/sec or V/cm]
c
c - -   Output:  ep(m+1,n) - real*8 array of phi component of E located at 
c                phi edge locations (TE grid) [G km/sec or V/cm]
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
      real*8 :: etcoe(m+1,n+1),epcoe(m+1,n+1)
c
c - - output variable declarations:
c
      real*8 :: ep(m+1,n),et(m,n+1)
c
c - - local variable declarations:
c
      integer :: i,j,jp1,jph,ip1,iph
c
c - - interpolate ep at theta edges:
c
      do i=1,m+1
         do j=1,n
            jph=j
            jp1=j+1
            ep(i,jph)=0.5d0*(epcoe(i,j)+epcoe(i,jp1))
         enddo
      enddo
c
c - - interpolate et at phi edges:
c
      do i=1,m
         iph=i
         ip1=i+1
         do j=1,n+1
            et(iph,j)=0.5d0*(etcoe(i,j)+etcoe(ip1,j))
         enddo
      enddo
c
      return
      end
