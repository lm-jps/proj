      subroutine emagpot_srf_ss(m,n,a,b,c,d,rsun,at,ap,bt,bp,emag)
c
c+
c - - Purpose:  To compute the magnetic energy for a potential magnetic field
c               with just a photospheric surface integral involving horizontal 
c               components of the vector potential and the magnetic field
c
c - - Usage:    call emagpot_srf_ss(m,n,a,b,c,d,rsun,at,ap,bt,bp,emag)
c
c - - Input:    m,n - integer number of cell centers in theta, phi directions, 
c               respectively.
c
c - - Input:    a,b - real*8 values of min, max of colatitude. a < b [radians]
c
c - - Input:    c,d - real*8 values of min, max of longitude. c < d [radians]
c
c - - Input:    rsun - real*8 value of radius of Sun [km].  Normally 6.96d5.
c
c - - Input:    at(m,n+1): real*8 array theta component vector potential on
c               PE grid at the photosphere [G km]
c
c - - Input:    ap(m+1,n): real*8 array of phi component vector potential on 
c               TE grid  at the photosphere [G km]
c
c - - Input:    bt(m+1,n): real*8 array of theta component pot. magnetic field
c               on TE grid at the photosphere [G]
c
c - - Input:    bp(m,n+1): real*8 array of phi component pot. magnetic field on 
c               PE grid at the photosphere [G]
c
c - - Output:   emag:  real*8 value of the pot. field magnetic energy [erg]
c
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
c
c - - input variables:
c
      integer :: m,n
      real*8 :: at(m,n+1),ap(m+1,n),bt(m+1,n),bp(m,n+1)
      real*8 :: a,b,c,d,rsun
c
c - - output variables:
c
      real*8 :: emag
c
c - - local variable declarations:
c
      real*8 :: dum,eightpim1,convfact,dphi,dtheta
      real*8 :: sinth(m+1),sinth_hlf(m)
      real*8 :: emagtmp,emagden(m,n)
      integer :: i,j,iph,jph
c
c - - declare pimach function (FISHPACK/FFTPACK) to compute pi
c
      real*8 :: pimach
c
      eightpim1=1./(8.d0*pimach(dum))
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - emag = integral of (- A_h X B_h) /(8 pi) dot rhat.
c - - Conversion factor to get emag  units into erg:  convert km^3 to cm^3:
c - - conversion factor = (1d5)^3 = 1d15
c
      convfact=1.d15
c 
c - - iph,jph indices at cell center, i,j,ip1,jp1 indices at cell edges
c
      do iph=1,m
         i=iph
         do jph=1,n
            j=jph
c - - Idea here is to average each component of -(A X B)/(8 pi) from edges 
c - - to cell  center.  Note doing it this way there are no sin(theta) 
c - - geometrical factors that come into the calculation, 
c - - just a straight average.
            emagden(iph,jph)= - 
     1             (0.5d0*(at(iph,j)*bp(iph,j)+at(iph,j+1)*bp(iph,j+1))
     2             -0.5d0*(ap(i+1,jph)*bt(i+1,jph)+ap(i,jph)*bt(i,jph)))
     3             *eightpim1
         enddo
      enddo
            emagtmp=0.d0
c
      do iph=1,m
         do jph=1,n
            emagtmp=emagtmp+
     1      emagden(iph,jph)*dtheta*sinth_hlf(iph)*dphi
         enddo
      enddo
c
c - - multiply emagtmp by rsun^2 * convfact
c
      emag=emagtmp*rsun*rsun*convfact
c
c
      return
      end
