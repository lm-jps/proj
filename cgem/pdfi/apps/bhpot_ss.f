       subroutine bhpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,scrb3d,btpot,bppot)
c
c+ - - Purpose: compute potential field value of bt,bp given the 3-d array
c              scrb3d by taking horizontal gradient of d scriptB / dr.
c              scrb3d is computed by subroutine scrbpot_ss.
c
c - -  Usage:  call bhpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,scrb3d,btpot,bppot)
c
c - -  Input:  m,n,p: integer values of numbers of cell centers in theta,
c              phi, and r directions
c - -  Input:  a,b,c,d:  real*8 values of min, max colatitude, min, max
c              values of longitude. [radians]
c - -  Input:  rsun,rssmrs: real*8 values of radius of sun, and distance 
c              from phot to source surface. [km] Normally, rsun=6.96d5
c - -  Input:  scrb3d(m,n,p+1): real*8 array of poloidal potential scribtb
c              [G km^2]
c - -  Output: btpot(m+1,n,p): real*8 array of theta-comp magnetic field [G]
c - -  Output: bppot(m,n+1,p): real*8 array of phi-comp magnetic field [G]
c - -  Note:   btpot,bppot are computed on theta and phi edges, and mid-way
c              between radial shells.
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
c
c - - input variables:
c
      integer :: m,n,p
      real*8 :: a,b,c,d,rsun,rssmrs
      real*8 :: scrb3d(m,n,p+1)
c 
c - - output variables:
c
      real*8 :: btpot(m+1,n,p)
      real*8 :: bppot(m,n+1,p)
c
c - - local variables:
c
      integer :: q,qmh,qph
      real*8 :: scrb(m+2,n+2),bdas(n),bdbs(n),bdcs(m),bdds(m)
      real*8 :: sinth(m+1),sinth_hlf(m),r(p+1),rh(p)
      real*8 :: gradp(m,n+1),gradt(m+1,n)
      real*8 :: dphi,dtheta,delr
c
      bdas(1:n)=0.d0
      bdbs(1:n)=0.d0
      bdcs(1:m)=0.d0
      bdds(1:m)=0.d0
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
      do qph=1,p
c - - get d script B / dr mid-way between shells, set equal to scrb array:
         q=qph
         scrb(2:m+1,2:n+1)=(scrb3d(1:m,1:n,q+1)-scrb3d(1:m,1:n,q))/delr
c - - fill in ghost zones in the mid-shell surface:
c - - (NOTE that for global solutions extending to poles, this may need to
c - - be altered).
         scrb(1,2:n+1)=scrb(2,2:n+1)-1.d0*dtheta*bdas(1:n)
         scrb(m+2,2:n+1)=scrb(m+1,2:n+1)+1.d0*dtheta*bdbs(1:n)
c        scrb(2:m+1,1)=scrb(2:m+1,2)-1.d0*dphi*bdcs(1:m)
c        scrb(2:m+1,n+2)=scrb(2:m+1,n+1)+1.d0*dphi*bdds(1:m)
c - -    Periodic BC in phi:
         scrb(2:m+1,1)=scrb(2:m+1,n+1)
         scrb(2:m+1,n+2)=scrb(2:m+1,2)
c
         call gradh_ce_ss(m,n,scrb,rh(qph),sinth_hlf,dtheta,dphi,gradt,
     1        gradp)
         btpot(1:m+1,1:n,qph)=gradt(1:m+1,1:n)
         bppot(1:m,1:n+1,qph)=gradp(1:m,1:n+1)
      end do
c
      return
      end
