      subroutine e_voxels3d_ss(m,n,rsun,sinth_hlf,dtheta,dphi,et,ep,
     1 scrb,dscrbdr,dr,ettop,etbot,eptop,epbot)
c
c+
c - - Purpose:  To compute the 3d electric field on the upper and lower rails 
c               of a layer of spherical voxels, given the horizontal electric
c               fields et,ep located mid-way between the layers, and ptd
c               solutions for scrb and dscrbdr.
c
c - - Method:   Use electric field computed in a spherical surface, plus 
c               the radial derivative of the time derivative of the poloidal 
c               potential in that same surface, to compute the electric field 
c               on all the rails of the voxels.  
c               It is assumed that the surface is placed half-way (in radius) 
c               between the top and the bottom layers of the voxels.  
c               We shall refer to this surface layer as the "photosphere".  
c               The radial derivative is derived from dscrbdr by taking a 
c               curl of dscrbdr times the rhat unit vector.
c
c - - Usage:    call e_voxels3d_ss(m,n,rsun,sinth_hlf,dtheta,dphi,et,ep,
c               scrb,dscrbdr,dr,ettop,etbot,eptop,epbot)
c
c - - Input:    m,n - integers describing the number of radial voxel face 
c               centers in the colatitude, and longitudinal directions, 
c               respectively.
c
c - - Input:    rsun:  real*8 value of radius of the Sun [km] Normally 6.96d5.
c
c - - Input:    sinth_hlf(m) - real*8 array of  sin(colatitude) computed 
c               at cell centers (computed from subroutine sinthta_ss)
c
c - - Input:    dtheta,dphi: real*8 values of cell thickness in colatitude, 
c               longitude directions [radians]

c - - Input:    et(m,n+1),ep(m+1,n) - real*8 arrays of theta and phi components
c               of cE computed from PTD or PDFI techniques at the photosphere.  
c               On PE, TE grid locations, respectively. [G km/s].
c
c - - Input:    scrb(m+2,n+2) - real*8 array of time derivative of the
c               poloidal potential returned from ptdsolve_ss subroutine.
c               CE grid locations plus ghost zones. [G km^2/sec]
c
c - - Input:    dscrbdr(m+2,n+2) - real*8 array of time derivative of radial 
c               derivative of poloidal potential returned from the 
c               ptdsolve_ss subroutine.  CE grid locations plus ghost zones.
c               [G km/sec]
c
c - - Input:    dr - real*8 value [km] of the depth of the radial legs
c               of the voxels.  The upper rails will be 0.5*dr above the 
c               photosphere, while the lower rails will be 0.5*dr below the 
c               photosphere.
c
c - - Output:   ettop(m,n+1),etbot(m,n+1) - real*8 arrays of E_theta, multiplied
c               by the speed of light, on the PE grid edges, at the top layer, 
c               and bottom layer, respectively [G km/s].
c
c - - Output:   eptop(m+1,n),epbot(m+1,n) - real*8 arrays of E_phi, multiplied 
c               by the speed of light, on the TE grid edges, at the top and 
c               bottom layers, respectively  [G km/s].
c
c - - NOTE:     the vertical rails will have the value of er along them, 
c               but since er is unchanged, we do not include er in the calling 
c               arguments.
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
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: et(m,n+1),ep(m+1,n)
      real*8 :: scrb(m+2,n+2),dscrbdr(m+2,n+2)
      real*8 :: rsun,dtheta,dphi,dr
      real*8 :: sinth_hlf(m)
c
c - - output variable declarations:
c
      real*8 :: ettop(m,n+1),etbot(m,n+1)
      real*8 :: eptop(m+1,n),epbot(m+1,n)
c
c - - local variable declarations:
c
      real*8 :: curlt(m,n+1),curlp(m+1,n),et_ptd(m,n+1),ep_ptd(m+1,n)
      real*8 :: detdr(m,n+1),depdr(m+1,n)
c
      call curl_psi_rhat_ce_ss(m,n,scrb,rsun,sinth_hlf,dtheta,dphi,
     1     curlt,curlp)
      et_ptd=-curlt
      ep_ptd=-curlp
      call curl_psi_rhat_ce_ss(m,n,dscrbdr,rsun,sinth_hlf,dtheta,dphi,
     1     curlt,curlp)
      detdr=-curlt
      depdr=-curlp
c
c - - dE_h/dr = -curl (dscrbdr rhat) - E_h/r.  [To see this, let E_h
c     = -curl (scrb rhat), and then evaluate 1/r (d/dr) r E_h.]  Then use
c     dE_h/dr and E_h to find E_h on the top and bottom rails of the voxels.
c
      ettop(1:m,1:n+1)=et(1:m,1:n+1)+0.5d0*dr*
     1     (detdr(1:m,1:n+1)-et_ptd(1:m,1:n+1)/rsun)
      etbot(1:m,1:n+1)=et(1:m,1:n+1)-0.5d0*dr*
     1     (detdr(1:m,1:n+1)-et_ptd(1:m,1:n+1)/rsun)
      eptop(1:m+1,1:n)=ep(1:m+1,1:n)+0.5d0*dr*
     1      (depdr(1:m+1,1:n)-ep_ptd(1:m+1,1:n)/rsun)
      epbot(1:m+1,1:n)=ep(1:m+1,1:n)-0.5d0*dr*
     1      (depdr(1:m+1,1:n)-ep_ptd(1:m+1,1:n)/rsun)
c
c - - we're done
c
      return
      end
