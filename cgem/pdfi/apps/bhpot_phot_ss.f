      subroutine bhpot_phot_ss(m,n,p,a,b,c,d,rsun,rssmrs,brce,scrb3d,
     1            btphot,bpphot)
c
c+ - - Purpose: compute potential field values of bt,bp 
c              at the photosphere given the 3-d array
c              scrb3d by taking horizontal gradient of d scriptB / dr 
c              (interpolated to the photosphere).  scrb3d is computed by 
c              subroutine scrbpot_ss.
c      
c
c  - - Usage:  call bhpot_phot_ss(m,n,p,a,b,c,d,rsun,rssmrs,brce,scrb3d,
c              btphot,bpphot)
c
c  - - Input:  m,n,p: integer values of numbers of cell centers in theta,
c              phi, and r directions
c
c  - - Input:  a,b,c,d:  real*8 values of min, max colatitude, min, max
c              values of longitude. [radians]
c
c  - - Input:  rsun,rssmrs: real*8 values of radius of sun, and distance 
c              from phot to source surface. [km] Normally, rsun=6.96d5.
c
c  - - Input:  brce(m,n): real*8 array of photospheric magnetic field values [G]
c
c  - - Input:  scrb3d(m,n,p+1): real*8 array of poloidal potential scribtb
c              [G km^2]
c
c  - - Output: btpot(m+1,n): real*8 array of theta-comp magnetic field at phot
c              [G]
c
c  - - Output: bppot(m,n+1): real*8 array of phi-comp magnetic field at phot
c              [G]
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
      real*8 :: brce(m,n),scrb3d(m,n,p+1)
c 
c - - output variables:
c
      real*8 :: btphot(m+1,n)
      real*8 :: bpphot(m,n+1)
c
c - - local variables:
c
      integer :: q
      real*8 :: scrb(m+2,n+2),bdas(n),bdbs(n),bdcs(m),bdds(m)
      real*8 :: scrbmh(m+2,n+2),scrbph(m+2,n+2)
      real*8 :: sinth(m+1),sinth_hlf(m)
      real*8 :: gradp(m,n+1),gradt(m+1,n),brcegh(m+2,n+2)
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
      q=1
c - - Get d script B /dr half a zone above the photosphere:
      scrb(2:m+1,2:n+1)=(scrb3d(1:m,1:n,q+1)-scrb3d(1:m,1:n,q))/delr
c
c - - brcegh is just brce with ghost zones filled in, using same boundary
c - - conditions as applied to scrb:
c
      brcegh(2:m+1,2:n+1)=brce(1:m,1:n)
c - - fill in ghost zones for scrb and brce:
      scrb(1,2:n+1)=scrb(2,2:n+1)-1.d0*dtheta*bdas(1:n)
      brcegh(1,2:n+1)=brcegh(2,2:n+1)-1.d0*dtheta*bdas(1:n)
      scrb(m+2,2:n+1)=scrb(m+1,2:n+1)+1.d0*dtheta*bdbs(1:n)
      brcegh(m+2,2:n+1)=brcegh(m+1,2:n+1)+1.d0*dtheta*bdbs(1:n)
c     scrb(2:m+1,1)=scrb(2:m+1,2)-1.d0*dphi*bdcs(1:m)
c     scrb(2:m+1,n+2)=scrb(2:m+1,n+1)+1.d0*dphi*bdds(1:m)
c - - Periodic BC in phi:
      scrb(2:m+1,1)=scrb(2:m+1,n+1)
      brcegh(2:m+1,1)=brcegh(2:m+1,n+1)
      scrb(2:m+1,n+2)=scrb(2:m+1,2)
      brcegh(2:m+1,n+2)=brcegh(2:m+1,2)
c
c - - Now use the fact that d/dr (d script B / dr) = B_r to find 
c - - d script B / dr half a zone below the photosphere:
c
      scrbmh(1:m+2,1:n+2)=scrb(1:m+2,1:n+2)-delr*brcegh(1:m+2,1:n+2)
c
c - - Now take average of scrb half a zone below phot and half a zone
c - - below it to interpolate scrb to the photosphere:
c
      scrbph(1:m+2,1:n+2)=0.5d0*(scrb(1:m+2,1:n+2)
     1 +scrbmh(1:m+2,1:n+2))
c
c - - Take gradient of scrb at photosphere to find btpot,bppot at phot.     
c
      call gradh_ce_ss(m,n,scrbph,rsun,sinth_hlf,dtheta,dphi,gradt,
     1        gradp)
         btphot(1:m+1,1:n)=gradt(1:m+1,1:n)
         bpphot(1:m,1:n+1)=gradp(1:m,1:n+1)
c
      return
      end
