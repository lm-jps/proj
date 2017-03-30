      subroutine delh2_sc(m,n,psi,rsun,sinth_gh,sinth_hlf_gh,dtheta,
     1 dphi,delh2)
c
c+
c - - Purpose:  Compute horizontal Laplacian in the centered (rather than
c               staggered) spherical coordinate system.
c               Use the same formulation as assumed
c               in Fishpack.
c - - Usage:  call delh2_sc(m,n,psi,rsun,sinth_gh,sinth_hlf_gh,dtheta,dphi,
c    1 delh2)
c
c - - Use standard 2nd order finite differences for interior
c - - cells. No one-sided derivatives used.  It is assumed that psi does
c - - include ghost-zones.  It is also assumed that sinth_gh (cell-edge values
c - - sin(theta) has ghost-zones, and that
c - - sinth_hlf_gh (cell-interior values of sin(theta))
c - - includes ghost-zones.  The output array, delh2, is computed only for
c - - active zones, so this array will have dimensions that are 2 smaller
c - - than the dimensions for the input array psi.
c - - Input:  rsun - units for the radius of the Sun.
c - - Input:  sinth_gh, sin(theta) computed at cell edges
c - - over the co-latitude range.  Ghost zones assumed included.
c - - Input:  sinth_hlf_gh, cell interior values of
c - - sin(theta) computed over the co-latitude range.  Ghostzones assumed
c - - included.
c - - Input:  dtheta,dphi - the spacing of cells in theta,phi directions
c - - Output: delh2 - The 2d horizontal Laplacian array, whose size is two less
c - - in each direction than the psi array.
c - - NOTE:  Plate Carree grid spacing (dtheta, dphi are constants) assumed!
c - - NOTE: The variables m,n are not the same as m,n in PDFI_SS!
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
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: psi(m+3,n+3),delh2(m+1,n+1),sinth_gh(m+3),
     1          sinth_hlf_gh(m+2)
c
      integer ntmax,npmax,j
      real*8 delth2inv,delph2inv,rsun2inv
      real*8, allocatable :: delh2tmp(:,:)

      ntmax=m+3
      npmax=n+3
c
      allocate(delh2tmp(ntmax,npmax))
c
      delth2inv=1.d0/(dtheta**2)
      delph2inv=1.d0/(dphi**2)
      rsun2inv=1.d0/(rsun**2)
c
      do j=2,npmax-1
         delh2tmp(2:ntmax-1,j)=delth2inv*(sinth_hlf_gh(2:ntmax-1)*
     1   (psi(3:ntmax,j)-psi(2:ntmax-1,j))-sinth_hlf_gh(1:ntmax-2)*
     2   (psi(2:ntmax-1,j)-psi(1:ntmax-2,j)))/sinth_gh(2:ntmax-1)+
     3   delph2inv*(psi(2:ntmax-1,j+1)+psi(2:ntmax-1,j-1)-
     4   2.d0*psi(2:ntmax-1,j))/(sinth_gh(2:ntmax-1))**2
      enddo
c 
c - - Trim Laplacian array from temporary array
c
      delh2(1:m+1,1:n+1)=delh2tmp(2:ntmax-1,2:npmax-1)
      delh2(1:m+1,1:n+1)=delh2(1:m+1,1:n+1)*rsun2inv
      deallocate(delh2tmp)
c
      return
      end
