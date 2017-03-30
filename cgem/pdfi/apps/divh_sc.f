      subroutine divh_sc(m,n,bt,bp,rsun,sinth,dtheta,dphi,div)
c
c+
c - - Purpose: Compute divergence of the vector bt,bp using centered grid
c - - (rather than staggered grid) formalism
c
c - - Usage:  call(divh_sc(m,n,bt,bp,rsun,sinth,dtheta,dphi,div)
c
c - - compute divergence of horizontal field B_h, given rectangular arrays
c - - of bt (B_theta), bp( B_phi).  Arrays assumed to have 1st index in
c - - theta direction, 2nd index in phi direction.
c - - Use standard 2nd order finite differences for interior
c - - cells, use one-sided derivatives at edges computed from 2nd order
c - - Lagrange interpolating polynomials.  rsun is assumed to be radius of
c - - the Sun, dtheta and dphi are the cell sizes in the theta and phi
c - - directions (in radians).  No ghost-zone values are used here, so
c - - make sure sinth includes no ghost-zones.
c - - Input: bt(m+1,n+1), bp(m+1,n+1):  (Input) Arrays of B_theta, B_phi.
c - - Input: rsun:  user input value for radius of Sun.
c - - Input: dtheta, dphi - distance between grid points in theta,phi
c - - Input sinth: array of sin(theta).  Must be computed from the 
c - - colatitudes corresponding to the 1st index of bt and bp.
c - - Output div(m+1,n+1): The two-dimensional divergence 
c - - NOTE:  Plate Carree grid spacing (dtheta, dphi are constants) assumed!
c
c - - NOTE:  Centered grid (not staggered!) assumed here.
c - - NOTE:  The use of m,n here is *not* consistent with our use of m,n in
c - - the PDFI_SS problem!  What we refer to as m-1 (PDFI_SS) will be the same
c - - as m+1 (centered grid).  Similarly n-1 (PDFI_SS) will be n+1 (centered
c - - grid).  Use Caution if using directly in staggered PDFI software!
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
      real*8 :: bt(m+1,n+1),bp(m+1,n+1),div(m+1,n+1),sinth(m+1)
c
      integer :: ntmax,npmax,i,j
      real*8 :: halfodt,halfodp,rsuninv
c
      ntmax=m+1
      npmax=n+1
c
      halfodt=0.5d0/dtheta
      halfodp=0.5d0/dphi
      rsuninv=1.d0/rsun
c
c - - divergence at interior points:  explicit loop over j, implicit loop
c - - over i
c
      do j=2,npmax-1
         div(2:ntmax-1,j)=
     1   halfodt*(bt(3:ntmax,j)*sinth(3:ntmax)-
     2   bt(1:ntmax-2,j)*sinth(1:ntmax-2))/sinth(2:ntmax-1)
     3   +halfodp*(bp(2:ntmax-1,j+1)-bp(2:ntmax-1,j-1))/
     4   sinth(2:ntmax-1)
      enddo

c
c - - one sided derivatives theta edges top and bottom:
c
c      write(6,*) 'divh_sc:59 sum(div) ',sum(div)
      do j=2,npmax-1
         div(1,j)=halfodp*(bp(1,j+1)-bp(1,j-1))/sinth(1)+
     1   (4.*bt(2,j)*sinth(2)-3.*bt(1,j)*sinth(1)-bt(3,j)
     2   *sinth(3))*halfodt/sinth(1)   
         div(ntmax,j)=halfodp*(bp(ntmax,j+1)-bp(ntmax,j-1))
     1   /sinth(ntmax)
     2   +(3.*bt(ntmax,j)*sinth(ntmax)+bt(ntmax-2,j)*sinth(ntmax-2)
     3   -4.*bt(ntmax-1,j)*sinth(ntmax-1))*halfodt/sinth(ntmax)
      enddo
c      write(6,*) 'divh_sc: 70 sum(div) ',sum(div)
c      write(6,*) 'divh_sc: 60 tot bt,bp',sum(bt),sum(bp)

c
c - - one sided derivatives of phi edges left and right:
      do i=2,ntmax-1
         div(i,1)=halfodt*(bt(i+1,1)*sinth(i+1)-bt(i-1,1)*sinth(i-1))
     1   /sinth(i)+(4.*bp(i,2)-3.*bp(i,1)-bp(i,3))*halfodp/sinth(i)
c
         div(i,npmax)=halfodt*(bt(i+1,npmax)*sinth(i+1)-bt(i-1,npmax)
     1   *sinth(i-1))/sinth(i)
     2   +(3.*bp(i,npmax)+bp(i,npmax-2)-4.*bp(i,npmax-1))*halfodp
     3   /sinth(i)
      enddo
c
c - - corners
c
      div(1,1)=halfodt*
     1   (4.*bt(2,1)*sinth(2)-3.*bt(1,1)*sinth(1)-bt(3,1)*sinth(3))
     2   /sinth(1)
     3   +halfodp*(4.*bp(1,2)-3.*bp(1,1)-bp(1,3))/sinth(1)
c
      div(1,npmax)=halfodt*
     1   (4.*bt(2,npmax)*sinth(2)-3.*bt(1,npmax)*sinth(1)-
     2   bt(3,npmax)*sinth(3))/sinth(1)+
     3   halfodp*(3.*bp(1,npmax)+bp(1,npmax-2)-4.*bp(1,npmax-1))
     4   /sinth(1)
c 
      div(ntmax,1)=halfodt*(3.*bt(ntmax,1)*sinth(ntmax)
     1   +bt(ntmax-2,1)*sinth(ntmax-2)-4.*bt(ntmax-1,1)*
     2   sinth(ntmax-1))/sinth(ntmax)+
     3   halfodp*(4.*bp(ntmax,2)-3.*bp(ntmax,1)-bp(ntmax,3))
     4   /sinth(ntmax)
c
      div(ntmax,npmax)=halfodt*(3.*bt(ntmax,npmax)*sinth(ntmax)+
     1   bt(ntmax-2,npmax)*sinth(ntmax-2)-4.*bt(ntmax-1,npmax)*
     2   sinth(ntmax-1))/sinth(ntmax)+
     3   halfodp*(3.*bp(ntmax,npmax)+bp(ntmax,npmax-2)-4.*
     4   bp(ntmax,npmax-1))/sinth(ntmax)
c   
      div(:,:)=div(:,:)*rsuninv
      return
      end
