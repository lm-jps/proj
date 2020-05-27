      subroutine enudge3d_ss(m,n,a,b,c,d,rsun,dr,btt,bpt,brt,
     1           ettop,etbot,eptop,epbot,er)
c+
c - - Purpose:  Given time derivatives of the 3 magnetic field components
c               (btt,bpt,brt) through the corresponding voxel faces in a 
c               layer of voxels, computes the electric field on all edges of 
c               the spherical voxels. 
c
c - - Usage:    call enudge3d_ss(m,n,a,b,c,d,rsun,dr,btt,bpt,brt,ettop,etbot,
c               eptop,epbot,er)
c
c - - Input:    m,n - integers denoting the numbers of cell interiors in the
c               colatitude and longitude directions, respectively.
c
c - - Input:    a,b - real*8 values of colatitude at the northern and southern
c               edges of the boundary [radians]
c
c - - Input:    c,d - real*8 variables of longitude at the boundary edges
c               [radians]
c
c - - Input:    rsun - real*8 value of the assumed solar radius [km]
c
c - - Input:    dr - real*8 value [km] of the voxel depth (in radial direction)
c
c - - Input:    btt(m+1,n) - real*8 array of time derivatives of B_theta
c               at TE grid locations [G/sec]
c
c - - Input:    bpt(m,n+1) - real*8 array of time derivatives of B_phi
c               at PE grid locations [G/sec]
c
c - - Input:    brt(m,n) - real*8 array of time derivatives of B_r
c               at CE grid locations [G/sec]
c
c - - Output:   ettop(m,n+1),etbot(m,n+1) - arrays of E_theta, multiplied by 
c               c (the speed of light), at top and bottom
c               edges of the voxels (PE grid) [G km/s].
c
c - - Output:   eptop(m+1,n),epbot(m+1,n) - arrays of E_phi, multiplied by 
c               c (the speed of light), at top and bottom edges of the voxels 
c               (TE grid) [G km/s].
c
c - - Output:   er(m+1,n+1) - array of E_r, multiplied by 
c               c (the speed of light), on the vertical rails of the voxels
c               (COE grid) [G km/s]
c
c - - NOTE:     This assumes the photospheric data is located halfway 
c               (in radius) through the layer of voxels.
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
c - - Input variable declarations:
c
      integer :: m,n
      real*8 :: a,b,c,d,rsun,dr
      real*8 :: btt(m+1,n),bpt(m,n+1),brt(m,n)
c
c - - Output variable declarations:
c
      real*8 :: ettop(m,n+1),etbot(m,n+1),eptop(m+1,n),epbot(m+1,n)
      real*8 :: er(m+1,n+1)
c
c - - Local variable declarations:
c
      real*8 dtheta,dphi
      real*8 sinth(m+1),sinth_hlf(m),scrb(m+2,n+2),dscrbdr(m+2,n+2)
      real*8 scrj(m+1,n+1)
      real*8 et(m,n+1),ep(m+1,n)
c
c - - compute sin(theta) arrays:
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - compute dtheta,dphi:
c
      dtheta=(b-a)/m
      dphi=(d-c)/n
c
c - - Solve ptd equations for time derivatives of scrb,dscrbdr,scrj:
c
      call ptdsolve_ss(m,n,btt,bpt,brt,rsun,sinth,sinth_hlf,a,b,c,d,
     1     scrb,dscrbdr,scrj)
c
c - - Compute PTD solution for E from scrb,scrj:
c
      call e_ptd_ss(m,n,scrb,scrj,rsun,sinth_hlf,dtheta,dphi,et,ep,er)
c
c - - Get E_h on top and bottom layers by calling e_voxels3d_ss:
c
      call e_voxels3d_ss(m,n,rsun,sinth_hlf,dtheta,dphi,et,ep,
     1     scrb,dscrbdr,dr,ettop,etbot,eptop,epbot)
c
c - - we're done
c
      return
      end
