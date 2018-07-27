      subroutine delh2_sc(m,n,psi,rsun,sinth_gh,sinth_hlf_gh,dtheta,
     1 dphi,delh2)
c
c+
c - - Purpose:  Compute horizontal Laplacian in the centered (rather than
c               staggered) spherical coordinate system.
c               Use the same formulation as assumed in Fishpack.
c
c               NOTE:  This subroutine *is not valid* for the staggered grid,
c               it is valid *only* within the context of the relax_psi_3d_ss
c               subroutine, which uses a centered grid formalism.
c
c - - Usage:    call delh2_sc(m,n,psi,rsun,sinth_gh,sinth_hlf_gh,dtheta,dphi,
c               delh2)
c
c - - Input:    m,n - integer values of the number of grid points in theta,
c               phi directions.  
c
c               NOTE:  These values of m,n *are not* the
c               same values of m,n that are in the staggered grid!  They are
c               each smaller by two than the staggered grid values of m,n.
c
c - - Input:    psi(m+3,n+3) - real*8 array of scalar potential [G km^2/sec]
c
c - - Input:    rsun - real*8 array of radius of the Sun [km].  Normally 6.96d5
c
c - - Input:    sinth_gh(m+3) - real*8 array of values of sin(colatitude),
c               including two ghost zone values.  Computed in subroutine
c               sinthta_sc, the centered grid version of sinthta.
c
c - - Input:    sinth_hlf_gh(m+2) - real*8 array of values of sin(colatitude),
c               at half-way points between the values in sinth_gh.
c               Computed in subroutine sinthta_sc, the centered grid version
c               of sinthta.
c
c - - Input:    dtheta,dphi - real*8 values of the angular separation in
c               colatitude and longitude, respectively. [radians]
c
c - - Output:   delh2(m+1,n+1) - real*8 array of the horizontal Laplacian
c               evaluated on interior grid points. [G/sec]
c
c     Note:     Uses standard 2nd order finite differences for interior
c               cells. No one-sided derivatives used.  It is assumed that psi 
c               does include ghost-zones.  
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
