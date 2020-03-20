      subroutine hm_ss(m,n,rsun,sinth_hlf,dtheta,dphi,et,ep,scrb,hm)
c
c+
c - - Purpose:  To compute the magnetic helicity flux density.
c
c - - Usage:    call hm_ss(m,n,rsun,sinth_hlf,dtheta,dphi,et,ep,scrb,hm)
c
c - - Input:    m,n - number of cell centers in theta, phi directions, resp.
c
c - - Input:    rsun -  real*8 value of radius of Sun [km].  Normally 6.96d5.
c
c - - Input:    sinth_hlf(m) - real*8 array of sin(colatitude) computed at 
c               cell centers.
c
c - - Input:    dtheta,dphi - real*8 values of cell thickness in colatitude, 
c               longitude. [radians]
c
c - - Input:    et(m,n+1): real*8 array of theta component electric field 
c               at PE grid locations [V/cm]
c
c - - Input:    ep(m+1,n): real*8 array of phi component electric field at 
c               TE grid locations [V/cm]
c
c - - Input:    scrb(m+2,n+2) - real*8 array of the poloidal potential 
c               at cell center locations (CE grid plus ghost zones). 
c               [G km^2]
c
c - - Output:   hm(m,n) - real*8 array of Magnetic helicity flux density 
c               computed on CE grid [Mx^2 cm^-2 s^-1] 
c
c - - NOTE:     If you want to call this subroutine using components of cE
c               in units of [G km/sec], divide by 1d3 to convert to [V/cm]
c               before calling this subroutine.
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
c - - declare pimach function (FISHPACK/FFTPACK) to compute pi
c     real*8 :: pimach
c
c - - input variables:
c
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: sinth_hlf(m)
      real*8 :: et(m,n+1),ep(m+1,n),scrb(m+2,n+2)
c
c - - output variables:
c
      real*8 :: hm(m,n)
c
c - - local variable declarations:
c
      real*8 :: convfact
      integer :: i,j,iph,jph
      real*8 :: at(m,n+1),ap(m+1,n)
      real*8 :: atbar,apbar,etbar,epbar
c
      call curl_psi_rhat_ce_ss(m,n,scrb,rsun,sinth_hlf,dtheta,
     1 dphi,at,ap)
c
c - - hm = 2 c E_h X A_h
c - - Conversion factor to get hm units into Mx^2 cm^-2 s^-1:  
c - - Vector potential will be output from curl_psi_rhat_ce_ss in units
c - - of G-km.  To convert A_h to G-cm, multiply by 1e5.
c - - To convert E_h in units of V/cm into c E_h in units of G cm s^-1
c - - multiply by a factor of 1e8.  Product of the two is a factor of 1d13:
c
      convfact=1.d13
c 
c - - iph,jph indices at cell center, i,j,ip1,jp1 indices at cell edges
c
      do iph=1,m
         i=iph
         do jph=1,n
            j=jph
c - - Idea here is to average each component of (cE_h X A_h)
c - - to cell  center.  Note that because each element of the cross product is
c - - not co-located, must first interpolate all quantities to cell center
c - - before computing cross product.
c - - Note that doing it this way there are no sin(theta) 
c - - geometrical factors that come into the calculation, 
c - - just a straight arithmetic average between each pair of edges.
            atbar=0.5d0*(at(iph,j)+at(iph,j+1))
            etbar=0.5d0*(et(iph,j)+et(iph,j+1))
            apbar=0.5d0*(ap(i,jph)+ap(i+1,jph))
            epbar=0.5d0*(ep(i,jph)+ep(i+1,jph))
            hm(iph,jph)=2.d0*(etbar*apbar-epbar*atbar)
     1             *convfact
         enddo
      enddo
c
      return
      end
