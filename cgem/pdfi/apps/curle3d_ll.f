      subroutine curle3d_ll(m,n,a,b,c,d,rsun,dr,elontop,elonbot,elattop,
     1 elatbot,erll,blont,blatt,brtll)
c
c+
c - - Purpose:  Given cE on all the rails of a layer of spherical voxels,
c               compute the curl of cE, and express the results as blont,
c               blatt, and brllt, the negative of the curl of E in the 
c               lon,lat, and radial directions.  The results can be 
c               viewed as the time derivative of the lon ,lat, and radial 
c               components of B through Faraday's law.
c - - Usage:    call curle3d_ll(m,n,a,b,c,d,rsun,dr,elontop,elonbot,elattop,
c               elatbot,erll,blont,blatt,brllt)
c - - Input:    m,n - integer values of the number of cell centers in the
c               colatitude (latitude) direction, and longitude direction, resp.
c - - Input:    a,b,c,d - real*8 values of the limits in colatitude 
c               (a and b; 0 <= a < b <= pi) and longitude (c <=0 < d
c               <= 2*pi) of the domain in spherical coordinates. [radians]
c - - Input:    rsun - real*8 value of the radius of the Sun [km].  Usually
c               6.96d5.
c - - Input:    dr - real*8 value [km] of the radial distance between the top 
c               and bottom surface of the voxels.  The top is assumed to be 
c               0.5*dr above the photosphere, and the bottom is 0.5*dr below 
c               the photosphere.
c - - Input:    elontop(n,m+1),elonbot(n,m+1) - real*8 arrays of the longitude
c               component of cE [G km/sec], at the top and bottom
c               faces, and on the TE grid, of the layer of voxels.
c - - Input:    elattop(n+1,m),elatbot(n+1,m) - real*8 arrays of lat component
c               of cE [G km/sec], at the top and bottom faces, and on the PE
c               grid, of the layer of voxels.
c - - Input:    erll(n+1,m+1) - real*8 array of the radial component of cE
c               [G km/sec], at the photosphere (mid-way in radius through
c               the layer of voxels).  erll is on COE grid in lon,lat order.
c - - Output:   blont(n+1,m),blatt(n,m+1),brtll(n,m) - real*8 arrays of the time
c               derivative of the magnetic field components in the theta,phi,
c               and radial directions, respectively, by setting B dot = 
c               - curl cE. [G/sec] Output in lon,lat order.
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
      integer :: m,n
      real*8 :: a,b,c,d,rsun,dr
      real*8 :: elontop(n,m+1),elonbot(n,m+1)
      real*8 :: elattop(n+1,m),elatbot(n+1,m)
      real*8 :: erll(n+1,m+1)
c
c - - Output variables:
c
      real*8 :: blont(n+1,m),blatt(n,m+1),brtll(n,m)
c
c - - Local variables:
c
      real*8 :: btt(m+1,n),bpt(m,n+1),brt(m,n)
      real*8 :: ettop(m,n+1),etbot(m,n+1),eptop(m+1,n),epbot(m+1,n)
      real*8 :: er(m+1,n+1)
c
c - - Transpose and flip input electric field arrays:
c
      call ehyeell2tp_ss(m,n,elontop,elattop,ettop,eptop)
      call ehyeell2tp_ss(m,n,elonbot,elatbot,etbot,epbot)
      call eryeell2tp_ss(m,n,erll,er)
c
c - - call spherical polar version of curle3d_ss:
c
      call curle3d_ss(m,n,a,b,c,d,rsun,dr,ettop,etbot,eptop,epbot,
     1     er,btt,bpt,brt)
c
c - - transpose and flip magnetic field time derivative arrays to lon,lat order:
c
      call bhyeetp2ll_ss(m,n,btt,bpt,blont,blatt)
      call bryeetp2ll_ss(m,n,brt,brtll)
c
c - - we're done
c
      return
      end
