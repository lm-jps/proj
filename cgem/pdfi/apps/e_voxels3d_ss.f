      subroutine e_voxels3d_ss(m,n,rsun,sinth_hlf,dtheta,dphi,et,ep,
     1 dscrbdr,dr,ettop,etbot,eptop,epbot)
c
c+
c - - Purpose:  To compute the 3d electric field on the upper and lower rails 
c     of a layer of spherical voxels.
c
c - - Method:  Use electric field computed in a spherical surface, plus 
c     the radial derivative of the time derivative of the poloidal potential 
c     in that same
c     surface, to compute the electric field on all the rails of the voxels.  
c     It is assumed that the surface is placed half-way (in radius) between
c     the top and the bottom layers of the voxels.  We shall refer to this
c     surface layer as the "photosphere".  The radial derivative is derived
c     from dscrbdr by taking a curl of dscrbdr times the rhat unit vector.
c
c - - Usage:  call e_voxels3d_ss(m,n,rsun,sinth_hlf,dtheta,dphi,et,ep,
c     dscrbdr,dr,ettop,etbot,eptop,epbot)
c
c - - Input:  m,n - integers describing the number of radial voxel face centers
c     in the colatitude, and longitudinal directions, respectively.
c - - Input:  rsun:  Assumed radius of the Sun [in km]
c - - Input:  sinth_hlf(m) : sin(colatitude) computed at cell centers
c             (computed from subroutine sinthta_ss)
c - - Input:  dtheta,dphi: cell thickness in colatitude, longitude directions
c - - Input:  et(m,n+1),ep(m+1,n) - real*8 electric field variables 
c     computed from PTD, or PDFI techniques at the photosphere.
c - - NOTE:  the vertical rails will have the value of er along them, but
c     since er is unchanged, we do not include er in the calling arguments.
c - - Input:  dscrbdr(m+2,n+2) - real*8 array returned from the ptdsolve_ss
c     subroutine.
c - - Input:  dr - a real*8 scalar [in km] that provides the depth of 
c     the radial legs
c     of the voxels.  The upper rails will be 0.5*dr above the photosphere,
c     while the lower rails will be 0.5*dr below the photosphere.
c - - Output:  ettop(m,n+1),etbot(m,n+1) - real*8 values of E_theta on the 
c     PE grid edges, at the top layer, and bottom layer, respectively.
c - - Output:  eptop(m+1,n),epbot(m+1,n) - real*8 values of E_phi on the TE grid
c     edges, at the top and bottom layers, respectively
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
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: et(m,n+1),ep(m+1,n)
      real*8 :: dscrbdr(m+2,n+2)
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
      real*8 :: curlt(m,n+1),curlp(m+1,n)
c
      call curl_psi_rhat_ce_ss(m,n,dscrbdr,rsun,sinth_hlf,dtheta,dphi,
     1     curlt,curlp)
c
c - - Note that dE_h/dr has opposite sign from curlt,curlp:
c
      ettop(1:m,1:n+1)=et(1:m,1:n+1)-0.5d0*dr*curlt(1:m,1:n+1)
      etbot(1:m,1:n+1)=et(1:m,1:n+1)+0.5d0*dr*curlt(1:m,1:n+1)
      eptop(1:m+1,1:n)=ep(1:m+1,1:n)-0.5d0*dr*curlp(1:m+1,1:n)
      epbot(1:m+1,1:n)=ep(1:m+1,1:n)+0.5d0*dr*curlp(1:m+1,1:n)
c
c - - we're done
c
      return
      end
