      subroutine divb3d_ss(m,n,p,a,b,c,d,rsun,rssmrs,
     1 btpot,bppot,brpot,divb)
c
c+ - - Purpose: Diagnostic, to compute divergence of the 3D magnetic field
c               within spherical wedge coronal domain, using potential 
c               field arrays of btpot,bppot,brpot.  Result is divergence of 
c               B as a 3D array.
c
c - -  Usage:  call divb3d_ss(m,n,p,a,b,c,d,rsun,rssmrs,btpot,
c              bppot,brpot,divb)
c
c - -  Input:  m,n,p -  integer values of numbers of cell centers in 
c              theta (colat), phi, and r directions
c
c - -  Input:  a,b - real*8 values of min, max of colatitude. a < b [radians]
c
c - -  Input:  c,d - real*8 values of min, max of longitude. c < d [radians]
c
c - -  Input:  rsun,rssmrs -  real*8 values of radius of sun, and distance 
c              from phot to source surface. [km].  Normally rsun=6.96d5.
c
c - -  Input:  btpot(m+1,n,p) - real*8 array of theta-comp magnetic field [G]
c
c - -  Input:  bppot(m,n+1,p) - real*8 array of phi-comp magnetic field [G]
c
c - -  Input:  brpot(m,n,p+1) - real*8 array of r-comp magnetic field [G]
c
c - - Output:  divb(m,n,p) - real*8 array of divergence of B in volume
c              [G / km]
c
c - -  Note:   btpot,bppot are computed on theta and phi edges, and mid-way
c              between radial shells. brpot is computed on radial shells, at
c              cell-center in theta and phi.  Divb evaluated at voxel centers. 
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
c
c - - input variables:
c
      integer :: m,n,p
      real*8 :: a,b,c,d,rsun,rssmrs
      real*8 :: btpot(m+1,n,p)
      real*8 :: bppot(m,n+1,p)
      real*8 :: brpot(m,n,p+1)
c 
c - - output variables:
c
      real*8 :: divb(m,n,p)
c
c - - local variables:
c
      integer :: q,qmh,qph
      real*8 :: sinth(m+1),sinth_hlf(m),r(p+1),rh(p)
      real*8 :: divbh(m,n),divbr(m,n)
      real*8 :: dphi,dtheta,delr
      real*8 :: rqph,rq,rqp1,rqph2inv,delrinv
c     integer :: iph,jph
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
      delr=rssmrs/dble(p)
      delrinv=1.d0/delr
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
c - - Loop to calculate divB.  
c
      do qph=1,p
         q=qph
c
c - - get radii at qph, q, and q+1:
c
         rqph=rh(qph)
         rqp1=r(q+1)
         rq=r(q)
c
c - - 1/rqph^2:
c
         rqph2inv=1.d0/(rqph**2)
c
c - - For each radial layer, first compute
c - - use divh_ce_ss to compute horizontal divergence of B_h and get result
c - - as a 2D array
c
         call divh_ce_ss(m,n,btpot(:,:,qph),bppot(:,:,qph),rqph,sinth,
     1        sinth_hlf,dtheta,dphi,divbh)
c
c - - Evaluate div B at each voxel center:
c
         divbr(:,:) = delrinv*rqph2inv*
     1   (brpot(:,:,q+1)*(rqp1**2) - brpot(:,:,q)*(rq**2))
         divb(:,:,qph)=divbr(:,:)+divbh(:,:)
c
      end do
c
c - - we're done
c
      return
      end
