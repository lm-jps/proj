      subroutine curlahpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,atpot,appot,
     1 btpot,bppot,brpot)
c
c+ - - Purpose: Compute potential field values of bt,bp,br given the 3-d arrays
c              atpot,appot of the potential-field vector potential, by taking
c              the curl of A_h.
c
c - -  Usage:  call curlahpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,atpot,appot,btpot,
c              bppot,brpot)
c
c - -  Input:  m,n,p -  integer values of numbers of cell centers in 
c              theta (colat), phi, and r directions
c - -  Input:  a,b - real*8 values of min, max of colatitude. a < b [radians]
c - -  Input:  c,d - real*8 values of min, max of longitude. c < d [radians]
c - -  Input:  rsun,rssmrs -  real*8 values of radius of sun, and distance 
c              from phot to source surface. [km].  Normally rsun=6.96d5.
c - -  Input:  atpot(m,n+1,p+1) - real*8 array of theta-comp vector potential
c              [G km]
c - -  Input:  appot(m+1,n,p+1) - real*8 array of phi-comp vector potential
c              [G km]
c - -  Output: btpot(m+1,n,p) - real*8 array of theta-comp magnetic field [G]
c - -  Output: bppot(m,n+1,p) - real*8 array of phi-comp magnetic field [G]
c - -  Output: brpot(m,n,p+1) - real*8 array of r-comp magnetic field [G]
c - -  Note:   btpot,bppot are computed on theta and phi edges, and mid-way
c              between radial shells. brpot is computed on radial shells, at
c              cell-center in theta and phi. atpot is located on radial shells,
c              on the PE (phi-edge) grid, and appot is located on radial shells,
c              on the TE (theta-edge) grid.
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
c - - input variables:
c
      integer :: m,n,p
      real*8 :: a,b,c,d,rsun,rssmrs
      real*8 :: atpot(m,n+1,p+1),appot(m+1,n,p+1)
c 
c - - output variables:
c
      real*8 :: btpot(m+1,n,p)
      real*8 :: bppot(m,n+1,p)
      real*8 :: brpot(m,n,p+1)
c
c - - local variables:
c
      integer :: q,qmh,qph
      real*8 :: sinth(m+1),sinth_hlf(m),r(p+1),rh(p)
      real*8 :: datdr(m,n+1),dapdr(m+1,n),curl(m,n)
      real*8 :: dphi,dtheta,delr
c 
      dtheta=(b-a)/m
      dphi=(d-c)/n
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
      delr=rssmrs/p
c
c - - define edge and cell-center radius arrays:
c
      r(1)=rsun
      do q=2,p+1
         qmh=q-1
         r(q)=r(1)+(q-1)*delr
         rh(qmh)=rsun+0.5*delr+(q-2)*delr
      end do
c
c - - loop to get B_t, B_p:
c
      do qph=1,p
c - -    compute necessary radial derivatives of A_theta and A_phi:
         q=qph
c        1/r d/dr(r A_theta):
         datdr(1:m,1:n+1)=(atpot(1:m,1:n+1,q+1)*r(q+1)
     1   -atpot(1:m,1:n+1,q)*r(q))/(rh(qph)*delr)
c        1/r d/dr(r A_phi):
         dapdr(1:m+1,1:n)=(appot(1:m+1,1:n,q+1)*r(q+1)
     1   -appot(1:m+1,1:n,q)*r(q))/(rh(qph)*delr)
c        B_t = - 1/r d/dr (r A_phi)
         btpot(1:m+1,1:n,qph) = -dapdr(1:m+1,1:n)
c        B_p  = 1/r d/dr (r A_theta)
         bppot(1:m,1:n+1,qph)=datdr(1:m,1:n+1)
      end do
c
c - - loop to get B_r:
c
      do q=1,p+1
         call curlh_ce_ss(m,n,atpot(1:m,1:n+1,q),appot(1:m+1,1:n,q),
     1        r(q),sinth,sinth_hlf,dtheta,dphi,curl(1:m,1:n))
         brpot(1:m,1:n,q)=curl(1:m,1:n)
      end do
c
c
      return
      end
