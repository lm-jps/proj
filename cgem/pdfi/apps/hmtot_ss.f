      subroutine hmtot_ss(m,n,rsun,sinth_hlf,dtheta,dphi,hm,hmtot)
c
c+
c - - Purpose:  To integrate the Helicity flux density over spherical wedge
c               domain to find helicity injection rate
c
c - - Usage:    call hmtot_ss(m,n,rsun,sinth_hlf,dtheta,dphi,hm,hmtot)
c
c - - Input:    m,n - integer number of cell centers in theta, phi directions, 
c               respectively.
c
c - - Input:    rsun - real*8 value of radius of Sun [km].  Normally 6.96d5.
c
c - - Input:    sinth_hlf(m) - real*8 array of sin(colatitude) computed 
c               at cell centers
c
c - - Input:    dtheta,dphi - real*8 values of cell thickness in colatitude, 
c               longitude [radians]
c
c - - Input:    hm(m,n) - real*8 array of Helicity flux density computed by
c               subroutine hm_ss on CE grid. [Mx^2 cm^-2 s^-1]
c
c - - Output:   hmtot: real*8 value of area-integrated Helicity injection rate 
c               [Mx^2 s^-1]
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
c - - input variables:
c
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: sinth_hlf(m)
      real*8 :: hm(m,n)
c
c - - output variable:
c
      real*8 :: hmtot
c
c - - local variable declarations:
c
      real*8 :: convfact
      integer :: iph,jph
c
c
c - - hmtot = sum deltaA_i,j hm_i,j
c - - Conversion factor to get hmtot units into Mx^2 s^-1:  
c - - hm in units of Mx^2 cm^-2 s^-1, but units of delta A in km^2.
c - - need to multiply by (1e5 cm/km)^2 = 1e10 to get area in cm^2.
c - - therefore, convfact=1.d10
c
      convfact=1.d10
c 
c - - iph,jph indices at cell center, i,j,ip1,jp1 indices at cell edges
c
      hmtot=0.d0
c
      do iph=1,m
         do jph=1,n
            hmtot=hmtot+
     1      hm(iph,jph)*dtheta*sinth_hlf(iph)*dphi
         enddo
      enddo
c
c - - multiply hmtot by rsun^2 * convfact
c
      hmtot=hmtot*rsun*rsun*convfact
c
      return
      end
