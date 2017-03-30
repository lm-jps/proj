      subroutine srtot_ss(m,n,rsun,sinth_hlf,dtheta,dphi,sr,srtot)
c
c+
c - - Purpose:  To integrate the Poynting flux density over photosphere
c
c - - Usage: call srtot_ss(m,n,rsun,sinth_hlf,dtheta,dphi,sr,srtot)
c
c - - Input:  m,n - number of cell centers in theta, phi directions, resp.
c - - Input:  rsun:  Assumed radius of the Sun [in km]
c - - Input:  sinth_hlf(m) : sin(colatitude) computed at cell centers
c - - Input:  dtheta,dphi: cell thickness in colatitude, longitude
c - - Input:  sr(m,n):  radial Poynting flux [erg cm^-2 s^-1] on CE grid
c - - Output: srtot: real*8 value of area-integrated Poynting flux [erg s^-1]
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
c - - declare pimach function (FISHPACK/FFTPACK) to compute pi
c     real*8 :: pimach
c
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: sinth_hlf(m)
      real*8 :: sr(m,n)
c
      real*8 :: srtot
c
c - - local variable declarations:
c
      real*8 :: convfact
      integer :: iph,jph
c
c
c - - srtot = sum deltaA_i,j sr_i,j
c - - Conversion factor to get srtot units into erg s^-1:  
c - - sr in units of erg cm^-2 s^-1, but units of delta A in km^2.
c - - need to multiply by (1e5 cm/km)^2 = 1e10
c - - therefore,convfact=1.d10
c
      convfact=1.d10
c 
c - - iph,jph indices at cell center, i,j,ip1,jp1 indices at cell edges
c
      srtot=0.d0
c
      do iph=1,m
         do jph=1,n
            srtot=srtot+
     1      sr(iph,jph)*dtheta*sinth_hlf(iph)*dphi
         enddo
      enddo
c
c - - multiply srtot by rsun^2 * convfact
c
      srtot=srtot*rsun*rsun*convfact
c
      return
      end
