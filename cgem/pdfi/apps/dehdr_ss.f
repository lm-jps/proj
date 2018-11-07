      subroutine dehdr_ss(m,n,rsun,sinth_hlf,dtheta,dphi,
     1 scrb,dscrbdr,detdr,depdr)
c
c+
c - - Purpose:  To compute the radial derivatives of the horizontal electric
c               field components, given the solutions for scrb and dscrbdr.
c
c - - Method:   Use ptd electric field computed in a spherical surface, plus 
c               the radial derivative of the time derivative of the poloidal 
c               potential in that same surface, to compute the electric field 
c               radial derivative at the photosphere.
c               The radial derivative is derived from dscrbdr by taking a 
c               curl of dscrbdr times the rhat unit vector plus a small 
c               spherical geometry contribution.  Contributions from the
c               gradient of a scalar are ignored [but these contribute nothing
c               to the curl of E].
c
c - - Usage:    call dehdr_ss(m,n,rsun,sinth_hlf,dtheta,dphi,
c               scrb,dscrbdr,detdr,depdr)
c
c               [Also see Note 2 below]
c
c - - Input:    m,n - integers describing the number of cell
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

c - - Input:    scrb(m+2,n+2) - real*8 array of time derivative of the
c               poloidal potential returned from ptdsolve_ss subroutine.
c               CE grid locations plus ghost zones. [G km^2/sec]
c
c - - Input:    dscrbdr(m+2,n+2) - real*8 array of time derivative of radial 
c               derivative of poloidal potential returned from the 
c               ptdsolve_ss subroutine.  CE grid locations plus ghost zones.
c               [G km/sec]
c
c - - Output:   detdr(m,n+1) - real*8 array of c*d E_theta / dr, 
c               on the PE grid edges, at the photosphere [G/s].
c
c - - Output:   depdr(m+1,n) - real*8 arrays of c*d E_phi / dr,
c               on the TE grid edges, at the photosphere [G/s].
c
c - - Note 1:   To convert detdr and depdr to units of [V/cm^2], 
c               multiply both arrays by 1d-8.
c
c - - Note 2:   You can use output arrays detdr and depdr to compute upper and
c               lower rail values for voxels that extend -0.5*dr below the
c               photosphere and 0.5*dr above the photosphere, along with 
c               photospheric values of the pdfi solution arrays et and ep
c               as follows:
c
c               ettop(:,:) = et(:,:) + 0.5*dr*detdr(:,:)
c               etbot(:,:) = et(:,:) - 0.5*dr*detdr(:,:)
c               eptop(:,:) = ep(:,:) + 0.5*dr*depdr(:,:)
c               epbot(:,:) = ep(:,:) - 0.5*dr*depdr(:,:)
c
c               This is an alternative to the use of e_voxels3d_ss.
c
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
      real*8 :: scrb(m+2,n+2),dscrbdr(m+2,n+2)
      real*8 :: rsun,dtheta,dphi
      real*8 :: sinth_hlf(m)
c
c - - output variable declarations:
c
      real*8 :: detdr(m,n+1)
      real*8 :: depdr(m+1,n)
c
c - - local variable declarations:
c
      real*8 :: curlt(m,n+1),curlp(m+1,n),et_ptd(m,n+1),ep_ptd(m+1,n)
      real*8 :: etdriv(m,n+1),epdriv(m+1,n)
c
      call curl_psi_rhat_ce_ss(m,n,scrb,rsun,sinth_hlf,dtheta,dphi,
     1     curlt,curlp)
      et_ptd=-curlt
      ep_ptd=-curlp
      call curl_psi_rhat_ce_ss(m,n,dscrbdr,rsun,sinth_hlf,dtheta,dphi,
     1     curlt,curlp)
      etdriv=-curlt
      epdriv=-curlp
c
c - - dE_h/dr = -curl (dscrbdr rhat) - E_h/r.  [To see this, let E_h
c     = -curl (scrb rhat), and then evaluate 1/r (d/dr) r E_h.]
c
      detdr(1:m,1:n+1)= etdriv(1:m,1:n+1)-et_ptd(1:m,1:n+1)/rsun
      depdr(1:m+1,1:n)= epdriv(1:m+1,1:n)-ep_ptd(1:m+1,1:n)/rsun
c
c - - we're done
c
      return
      end
