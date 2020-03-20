      subroutine sr_ss(m,n,et,ep,bt,bp,sr)
c
c+
c - - Purpose:  To compute the Poynting flux of magnetic energy
c
c - - Usage:    call sr_ss(m,n,et,ep,bt,bp,sr)
c
c - - Input:    m,n - integer number of cell centers in theta, phi directions, 
c               respectively.
c
c - - Input:    et(m,n+1): real*8 array theta component electric field on 
c               PE grid [V/cm]
c
c - - Input:    ep(m+1,n): real*8 array of phi component electric field on 
c               TE grid [V/cm]
c
c - - Input:    bt(m+1,n): real*8 array of theta component magnetic field
c               on TE grid [G]
c
c - - Input:    bp(m,n+1): real*8 array of phi component magnetic field on 
c               PE grid [G]
c
c - - Output:   sr(m,n):  real*8 array of Poynting flux computed on CE grid
c               [erg/(cm^2 sec)]
c
c - - NOTE:     If you want to use this subroutine with components of cE in
c               units of [G km/sec], you will first need to divide et,ep by 1d3
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
c
c - - input variables:
c
      integer :: m,n
      real*8 :: et(m,n+1),ep(m+1,n),bt(m+1,n),bp(m,n+1)
c
c - - output variables:
c
      real*8 :: sr(m,n)
c
c - - local variable declarations:
c
      real*8 :: dum,fourpim1,convfact
      integer :: i,j,iph,jph
c
c - - declare pimach function (FISHPACK/FFTPACK) to compute pi
c
      real*8 :: pimach
c
      fourpim1=1./(4.d0*pimach(dum))
c
c - - sr = c E_h X B_h /(4 pi).
c - - Conversion factor to get sr units into erg cm^2 s^-1:  One factor of 1d3
c - - to convert E in V/cm to cE in units of G km/sec, and
c - - another factor of 1d5 to convert cE in G km/sec to cE in G cm/sec,
c - - total conversion factor = 1d3*1d5=1d8.
c
      convfact=1.d8
c 
c - - iph,jph indices at cell center, i,j,ip1,jp1 indices at cell edges
c
      do iph=1,m
         i=iph
         do jph=1,n
            j=jph
c - - Idea here is to average each component of (cE X B)/(4 pi) from edges 
c - - to cell  center.  Note doing it this way there are no sin(theta) 
c - - geometrical factors that come into the calculation, 
c - - just a straight average.
            sr(iph,jph)=
     1             (0.5d0*(et(iph,j)*bp(iph,j)+et(iph,j+1)*bp(iph,j+1))
     2             -0.5d0*(ep(i+1,jph)*bt(i+1,jph)+ep(i,jph)*bt(i,jph)))
     3             *fourpim1*convfact
         enddo
      enddo
c
      return
      end
