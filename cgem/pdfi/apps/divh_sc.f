      subroutine divh_sc(m,n,et,ep,rsun,sinth,dtheta,dphi,div)
c
c+
c - - Purpose: Compute divergence of the vector et,ep using centered grid
c              (rather than staggered grid) formalism
c
c              NOTE:  This subroutine is valid *ONLY* within the centered
c              grid context used in subroutine relax_psi_3d_ss.f  *DO NOT USE*
c              with staggered grid variables.
c
c - - Usage:   call(divh_sc(m,n,et,ep,rsun,sinth,dtheta,dphi,div)
c 
c - - Method:  Use centered difference expressions for interior points,
c              one-sided difference approximations to derivatives for
c              edge and corner points derived from 2nd order Lagrange
c              interpolating polynomials.
c
c - - Input:   m,n - integer number of gridpoints in colatitude and longitude
c              directions, respectively.
c
c              NOTE:  m,n here are smaller by two than the value of m,n
c              used in the standard staggered mesh formalism of PDFI_SS.
c
c - - Input:   et(m+1,n+1) - real*8 array of theta component of c times
c              the electric field [G km/sec]
c
c - - Input:   ep(m+1,n+1) - real*8 array of phi component of c times the
c              electric field [G km/sec]
c
c - - Input:   rsun - real*8 value of radius of Sun. [km] Normally 6.96d5.
c
c - - Input:   sinth(m+1) - real*8 array of sin(colatitude), for the 
c              colatitude grid range.  Computed by subroutine sinthta_sc.
c
c - - Input:   dtheta,dphi - real*8 values of the angular gridpoint separation
c              in the colatitude, and longitude directions, respectively.
c              [radians]
c 
c - - Output:  div(m+1,n+1) - real*8 array of values of the horizontal
c              divergence of cE_h, the horizontal components of cE.
c              [G/s]
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
c
      implicit none
c
c - - input variables:
c
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: et(m+1,n+1),ep(m+1,n+1),sinth(m+1)
c
c - - output variables:
c
      real*8 :: div(m+1,n+1)
c
c - - local variables:
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
     1   halfodt*(et(3:ntmax,j)*sinth(3:ntmax)-
     2   et(1:ntmax-2,j)*sinth(1:ntmax-2))/sinth(2:ntmax-1)
     3   +halfodp*(ep(2:ntmax-1,j+1)-ep(2:ntmax-1,j-1))/
     4   sinth(2:ntmax-1)
      enddo

c
c - - one sided derivatives theta edges top and bottom:
c
      do j=2,npmax-1
         div(1,j)=halfodp*(ep(1,j+1)-ep(1,j-1))/sinth(1)+
     1   (4.*et(2,j)*sinth(2)-3.*et(1,j)*sinth(1)-et(3,j)
     2   *sinth(3))*halfodt/sinth(1)   
         div(ntmax,j)=halfodp*(ep(ntmax,j+1)-ep(ntmax,j-1))
     1   /sinth(ntmax)
     2   +(3.*et(ntmax,j)*sinth(ntmax)+et(ntmax-2,j)*sinth(ntmax-2)
     3   -4.*et(ntmax-1,j)*sinth(ntmax-1))*halfodt/sinth(ntmax)
      enddo
c
c - - one sided derivatives of phi edges left and right:
c
      do i=2,ntmax-1
         div(i,1)=halfodt*(et(i+1,1)*sinth(i+1)-et(i-1,1)*sinth(i-1))
     1   /sinth(i)+(4.*ep(i,2)-3.*ep(i,1)-ep(i,3))*halfodp/sinth(i)
c
         div(i,npmax)=halfodt*(et(i+1,npmax)*sinth(i+1)-et(i-1,npmax)
     1   *sinth(i-1))/sinth(i)
     2   +(3.*ep(i,npmax)+ep(i,npmax-2)-4.*ep(i,npmax-1))*halfodp
     3   /sinth(i)
      enddo
c
c - - corners (one-sided derivs in both directions)
c
      div(1,1)=halfodt*
     1   (4.*et(2,1)*sinth(2)-3.*et(1,1)*sinth(1)-et(3,1)*sinth(3))
     2   /sinth(1)
     3   +halfodp*(4.*ep(1,2)-3.*ep(1,1)-ep(1,3))/sinth(1)
c
      div(1,npmax)=halfodt*
     1   (4.*et(2,npmax)*sinth(2)-3.*et(1,npmax)*sinth(1)-
     2   et(3,npmax)*sinth(3))/sinth(1)+
     3   halfodp*(3.*ep(1,npmax)+ep(1,npmax-2)-4.*ep(1,npmax-1))
     4   /sinth(1)
c 
      div(ntmax,1)=halfodt*(3.*et(ntmax,1)*sinth(ntmax)
     1   +et(ntmax-2,1)*sinth(ntmax-2)-4.*et(ntmax-1,1)*
     2   sinth(ntmax-1))/sinth(ntmax)+
     3   halfodp*(4.*ep(ntmax,2)-3.*ep(ntmax,1)-ep(ntmax,3))
     4   /sinth(ntmax)
c
      div(ntmax,npmax)=halfodt*(3.*et(ntmax,npmax)*sinth(ntmax)+
     1   et(ntmax-2,npmax)*sinth(ntmax-2)-4.*et(ntmax-1,npmax)*
     2   sinth(ntmax-1))/sinth(ntmax)+
     3   halfodp*(3.*ep(ntmax,npmax)+ep(ntmax,npmax-2)-4.*
     4   ep(ntmax,npmax-1))/sinth(ntmax)
c   
      div(:,:)=div(:,:)*rsuninv
      return
      end
