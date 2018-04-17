      subroutine curle_ss(m,n,rsun,a,b,c,d,et,ep,scrj,dscrbdr,btt,bpt,
     1 brt)
c
c+
c - - Purpose: Compute 3 components of the curl of cE.
c              This can be used to check inductivity of solutions for E.
c              This subroutine uses photospheric values of et,ep, plus
c              arrays of toroidal potential plus radial derivative of
c              poloidal potential, to compute all 3 components of the curl
c              of cE.  This is in contrast to curle3d_ss, which uses E
c              values on all the rails of the spherical voxels to compute
c              the 3 components of the curl.
c
c - - Usage:   call curle_ss(m,n,rsun,a,b,c,d,et,ep,scrj,dscrbdr,btt,bpt,brt)
c
c - - Input:   m,n - number of cell centers in theta, phi directions, resp.
c - - Input:   rsun - real*8 value of assumed radius of Sun [km].
c              Normally 6.96d5.
c - - Input:   a,b - real*8 colatitude limits of wedge domain in radians (a < b)
c              [radians]
c - - Input:   c,d - real*8 longitude limits of wedge domain in radians (c < d)
c              [radians]
c - - Input:   et(m,n+1) - real*8 array of the theta component of electric 
c              field (cE) multiplied by the speed of light [G km/s]
c - - Input:   ep(m+1,n) - real*8 array of the phi component of electric field 
c              (cE) multiplied by the speed of light [G km/s]
c - - Input:   scrj(m+1,n+1) - real*8 array of time derivative of the toroidal
c              potential returned from ptdsolve_ss [G km/s]
c - - Input:   dscrbdr(m+2,n+2) - real*8 array of the time deriv. of the radial
c              derivative of the poloidal potential returned from ptdsolve_ss.
c              [G km/s]
c - - Output:  btt(m+1,n) - real*8 array of theta component of curl c*E [G/s]
c - - Output:  bpt(m,n+1) - real*8 array of phi component of curl c*E [G/s]
c - - Output:  brt(m,n) - real*8 array of r component of curl c*E [G/s]
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
c - - Input variables:
c
      integer :: m,n
      real*8 :: rsun,a,b,c,d
      real*8 :: et(m,n+1),ep(m+1,n),scrj(m+1,n+1),dscrbdr(m+2,n+2)
c
c - - Output variables:
c
      real*8 :: btt(m+1,n),bpt(m,n+1),brt(m,n)
c
c - - Local variables:
c
      real*8 :: curlr(m,n),curlt(m+1,n),curlp(m,n+1)
      real*8 :: gradt(m+1,n),gradp(m,n+1)
      real*8 :: sinth(m+1),sinth_hlf(m),dtheta,dphi
c
      dtheta=(b-a)/m
      dphi=(d-c)/n
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
      call curlh_ce_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,
     1 dphi, curlr)
c - - note that E_r (inductive) is minus scrj
      call curl_psi_rhat_co_ss(m,n,scrj,rsun,sinth,dtheta,dphi,curlt,
     1 curlp)
c
      call gradh_ce_ss(m,n,dscrbdr,rsun,sinth_hlf,dtheta,dphi,gradt,
     1 gradp)
c - - these quantities should be minus the time derivatives of the magnetic 
c - - field components
      brt(1:m,1:n)=curlr(1:m,1:n)
      btt(1:m+1,1:n)=-curlt(1:m+1,1:n)-gradt(1:m+1,1:n)
      bpt(1:m,1:n+1)=-curlp(1:m,1:n+1)-gradp(1:m,1:n+1)
c
      return
      end
