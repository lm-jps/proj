      subroutine br_voxels3d_ss(m,n,rsun,sinth,sinth_hlf,dtheta,dphi,
     1 bt,bp,br,dr,brtop,brbot)
c
c+
c - - Purpose: To compute the radial magnetic field (or its time derivative)
c             on the upper and lower faces of a layer of spherical voxels, 
c             given magnetic field components halfway in radius between 
c             the faces.
c
c - - Method: Use div B=0 (or div B dot = 0) to find radial derivative of 
c             Br (or Br dot) from measured 
c             values of divh cdot Bh:  (d/dr)(r**2 B_r) = -r**2 div_h B_h,
c             where B_h is the horizontal magnetic field (or its time 
c             derivative).  Once radial derivative of Br (or Br dot) is known 
c             at photosphere, use it to estimate radial magnetic field half 
c             a voxel above and below the photosphere.
c             It is assumed that the photosphere is placed half-way 
c             (in radius) between the top and the bottom layers of the voxels.
c
c - - Usage:  call br_voxels3d_ss(m,n,rsun,sinth,sinth_hlf,dtheta,dphi,bt,bp,br,
c             dr,brtop,brbot)
c
c - - Input:  m,n - integers describing the number of radial voxel face centers
c             in the colatitude, and longitudinal directions, respectively.
c
c - - Input:  rsun:  real*8 radius of the Sun [km].  Normally 6.96e5.
c
c - - Input:  sinth(m+1) - real*8 array of sin(colatitude) computed at 
c             theta cell edges (computed from subroutine sinthta_ss)
c
c - - Input:  sinth_hlf(m) - real*8 array of sin(colatitude) computed at 
c             theta cell centers (computed from subroutine sinthta_ss)
c
c - - Input:  dtheta,dphi: - real*8 cell thickness in colatitude, longitude 
c             directions [radians]
c
c - - Input:  bt(m+1,n),bp(m,n+1),br(m,n) - real*8 magnetic field variables 
c             (or their time derivatives) on staggered mesh 
c             (br on CE grid, bt on TE grid, bp on PE grid). [G or G/sec]
c
c - - Input:  dr - a real*8 scalar [km] that provides the depth of 
c             the radial legs of the voxels.  The upper face will be 0.5*dr 
c             above the photosphere, while the lower face will be 0.5*dr 
c             below the photosphere.
c
c - - Output: brtop(m,n),brbot(m,n) - real*8 values of B_r (or B_r dot)
c             on the CE grid, at the top layer, and bottom layer, respectively.
c             [G or G/sec]
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
      real*8 :: bt(m+1,n),bp(m,n+1),br(m,n)
      real*8 :: rsun,dtheta,dphi,dr
      real*8 :: sinth_hlf(m),sinth(m+1)
c
c - - output variable declarations:
c
      real*8 :: brtop(m,n),brbot(m,n)
c
c - - local variable declarations:
c
      real*8 :: divbh(m,n)
c
      call divh_ce_ss(m,n,bt,bp,rsun,sinth,sinth_hlf,dtheta,dphi,
     1     divbh)
c
c - - Note that dB_r/dr has opposite sign from divbh:
c
      brtop(1:m,1:n)=((rsun**2)*br(1:m,1:n)-0.5d0*dr*(rsun**2)*
     1 divbh(1:m,1:n))/((rsun+0.5d0*dr)**2)
      brbot(1:m,1:n)=((rsun**2)*br(1:m,1:n)+0.5d0*dr*(rsun**2)*
     1 divbh(1:m,1:n))/((rsun-0.5d0*dr)**2)
c
c - - we're done
c
      return
      end
