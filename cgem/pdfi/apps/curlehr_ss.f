      subroutine curlehr_ss(m,n,rsun,a,b,c,d,et,ep,brt)
c
c+
c - - Purpose: Compute radial component of the curl of c*E_h.
c     This can be used to check inductivity of solutions for E_h.
c
c - - Usage:  call curlehr_ss(m,n,rsun,a,b,c,d,et,ep,brt)
c
c - - Input:  m,n - number of cell centers in theta, phi directions, resp.
c - - Input:  rsun - assumed radius of Sun in km.
c - - Input:  a,b - colatitude limits of wedge domain in radians (a < b)
c - - Input:  c,d - longitude limits of wedge domain in radians (c < d)
c - - Input: et(m,n+1) - theta component of ptd electric field
c - - Input: ep(m+1,n) - phi component of ptd electric field
c - - Output: brt(m,n) - the radial component of curl c*E_h
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
c - - Input variables:
c
      integer :: m,n
      real*8 :: rsun,a,b,c,d
      real*8 :: et(m,n+1),ep(m+1,n)
c
c - - Output variables:
c
      real*8 :: brt(m,n)
c
c - - Local variables:
c
      real*8 :: curlr(m,n)
      real*8 :: sinth(m+1),sinth_hlf(m),dtheta,dphi
c
      dtheta=(b-a)/m
      dphi=(d-c)/n
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
      call curlh_ce_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,
     1 dphi, curlr)
c
      brt(1:m,1:n) = curlr(1:m,1:n)
c
      return
      end
