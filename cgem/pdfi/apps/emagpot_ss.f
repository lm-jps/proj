      subroutine emagpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,
     1 btpot,bppot,brpot,emag)
c
c+ - - Purpose: Compute potential magnetic energy within spherical wedge
c               coronal domain, using potential field arrays of bt,bp,br.  
c               Result is volume integral of B^2/(8 pi). 
c
c - -  Usage:  call emagpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,btpot,
c              bppot,brpot,emag)
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
c - - Output:  emag - real*8 value of total magnetic energy in coronal volume
c              [erg]
c
c - -  Note:   btpot,bppot are computed on theta and phi edges, and mid-way
c              between radial shells. brpot is computed on radial shells, at
c              cell-center in theta and phi. 
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
      real*8 :: emag
c
c - - local variables:
c
      integer :: q,qmh,qph
      real*8 :: sinth(m+1),sinth_hlf(m),r(p+1),rh(p)
      real*8 :: dphi,dtheta,delr,pi,dum,pi8inv,emagtmp
      real*8 :: bp,bt,br,dv
      integer :: i,iph,j,jph
c 
c - - function declaration for pimach (from FISHPACK)
c
      real*8 :: pimach
c
      pi=pimach(dum)
      pi8inv=1.d0/(8.d0*pi)
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
      emagtmp=0.d0
c
c - - 3D loop to do volume integral:
c
      do qph=1,p
         q=qph
c
         do jph=1,n
           j=jph
c
           do iph=1,m
             i=iph
c
c - - get magnetic field components at center of voxel from voxel faces:
c
             bt=0.5*(btpot(i+1,jph,qph)+btpot(i,jph,qph))
             bp=0.5*(bppot(iph,j+1,qph)+bppot(iph,j,qph))
             br=0.5*(brpot(iph,jph,q+1)+0.5*brpot(iph,jph,q))
c
c - - dv is volume of voxel (factor of 1d15 converts km^3 to cm^3):
c
             dv=(rh(qph)**2)*dtheta*dphi*sinth_hlf(iph)*delr
     1          *1d15
c
c - - emagtmp is running total of magnetic energy
c
             emagtmp=emagtmp+pi8inv*dv*(bt**2+bp**2+br**2)
           end do
         end do
      end do
c
c - - assign output value to emagtmp:
c
      emag=emagtmp
c
c - - we're done
c
      return
      end
