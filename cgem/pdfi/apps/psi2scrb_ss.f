      subroutine psi2scrb_ss(m,n,p,rsun,rssmrs,psi3d,scrb3d)
c
c+ - Purpose: Integrate 3-d scalar potential field within a spherical wedge
c             domain, ie fixed co-latitude and azimuth boundaries, 
c             plus inner and outer radial shells at Rsun and
c             Rsun plus distance to source surface, from the source surface
c             toward the photosphere.
c             Solution (poloidal potential scrb3d) is returned
c             on radial voxel faces.  From scrb3d, vector potential and magnetic
c             field components can be recovered.
c
c - - Method: Interpolate scalar potential
c             in radial direction to voxel centers using 2nd and
c             3rd order Lagrange
c             interpolation polynomials, and then
c             use trapezoidal rule to do the integration.
c             Results in conservative relationship 
c             psi(q + 1/2)=-(scrb(q+1)-scrb(q))/ (Delta r)
c
c   - Usage:  call psi2scrb_ss(m,n,p,rsun,rssmrs,psi3d,scrb3d)
c
c - - Input:  m,n,p - integer number of cell interiors in the theta,phi,r 
c             directions, respectively.
c
c - - Input:  rsun - real*8 value of radius of Sun [km].  Normally 6.96d5.
c
c - - Input:  rssmrs - real*8 distance from surface to source surface [km]
c
c - - Input:  psi3d(m,n,p+1) - real*8 3D array of psi (scalar potential)
c             for a potential magnetic field, evaluated at at cell centers in 
c             theta,phi, and at radial shells. [G km]
c
c - - Output: scrb3d(m,n,p+1) - real*8 3D array of scrb (poloidal potential)
c             evaluated at radial voxel faces [G km^2]
c
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
      real*8 :: rsun,rssmrs
      real*8 :: psi3d(m,n,p+1)
c
c - - output variables:
c
      real*8 :: scrb3d(m,n,p+1)
c
c - - Internal variables and arrays:
c
      integer :: q,qph
      real*8 :: delr
      real*8 :: r(p+1),rce(p)
      real*8 :: wt0i,wt1i,wt2i,wt3i,wt1p,wt2p,wt3p,wt0s,wt1s,wt2s
      real*8 :: wthlf
      real*8 :: psii3d(m,n)
c
c - - radial spacing:
c
      delr=rssmrs/p
c
c - - define edge and cell-center radius arrays:
c
      r(1)=rsun
      do q=2,p+1
         r(q)=r(1)+(q-1)*delr
         rce(q-1)=rsun+0.5*delr+delr*(q-2)
      end do
c
c - - define outermost radial point of poloidal potential:
c - - technically the source surface is half a voxel beyond r(p+1),
c - - so the following expression accounts for the half-voxel contribution,
c - - approximately:
c
      scrb3d(:,:,p+1)=0.5d0*delr*psi3d(:,:,p+1)
c
c - - interpolation weights for interior radial points:
c
      wt0i=-1.d0/16.d0
      wt1i=9.d0/16.d0
      wt2i=9.d0/16.d0
      wt3i=-1.d0/16.d0
c
c - - interpolation weights for radial point near photosphere
c
      wt1p=3.d0/8.d0
      wt2p=3.d0/4.d0
      wt3p=-1.d0/8.d0
c
c - - interpolation weights for radial point near source-surface:
c
      wt0s=-1.d0/8.d0
      wt1s=3.d0/4.d0
      wt2s=3.d0/8.d0
c
c     wthlf=0.5d0
c
c - - Define psi interpolated to midpoints of voxels in radius:
c
      do qph=p,1,-1
c Note loop is going in negative direction, from source-surf to photosphere:
        q=qph
        if(qph .eq. 1) then
c - photospheric boundary: (2nd order Lagrange interp. polynomial)
          psii3d(:,:)=wt1p*psi3d(:,:,q)+wt2p*psi3d(:,:,q+1)
     1                 +wt3p*psi3d(:,:,q+2)
        else if (qph .eq. p) then
c - source surface "boundary" (2nd order Lagrange interp. polynomial)
          psii3d(:,:)=wt0s*psi3d(:,:,q-1)+wt1s*psi3d(:,:,q)
     1                 +wt2s*psi3d(:,:,q+1)
        else
c - interior locations: (avg of two 2nd order Lagrange interp. polynomials)
c   (which is also the 3rd order Lagrange interp. polynomial):
          psii3d(:,:)=wt0i*psi3d(:,:,q-1)+wt1i*psi3d(:,:,q)
     1                 +wt2i*psi3d(:,:,q+1)+wt3i*psi3d(:,:,q+2)
        endif
c
c - - Use recursion relationship
c - - to define rest of points, using integral relationship with psi evaluated
c - - at radial voxel centers:
c
         scrb3d(:,:,q)=scrb3d(:,:,q+1)+delr*psii3d(:,:)
c
      end do
c - - we're done:
c
      return
      end
