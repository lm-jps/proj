      subroutine divh_ce_ss(m,n,bt,bp,rsun,sinth,sinth_hlf,dtheta,
     1 dphi,div)
c
c+ 
c - -  Purpose: Compute horizontal divergence of the vector with components 
c               bt (TE grid), bp (PE grid) and evaluated at cell centers 
c               (CE grid).
c
c - -   Usage:  call divh_ce_ss (m,n,bt,bp,rsun,sinth,sinth_hlf,dtheta,dphi,div)
c
c - -   Input:  m,n - integer number of cell centers in theta and phi 
c               directions, resp.
c
c - -   Input:  bt(m+1,n) - real*8 array of theta component of B, or its time
c               derivative, located on TE grid (theta edges). [G or G/sec]
c
c - -   Input:  bp(m,n+1) - real*8 array of phi component of B or its time
c               derivative, located on PE grid (phi edges). [G or G/sec]
c
c - -   Input:  rsun - real*8 value of radius of Sun [km]. Normally 6.96d5
c
c - -   Input:  sinth(m+1) - real*8 array of theta cell-edge values of 
c               sin(theta), spanning the domain of the problem.
c
c - -   Input:  sinth_hlf(m) - real*8 array of theta-cell center values of
c               sin(theta).
c
c - -   Input:  dtheta,dphi - real*8 values of distance between theta edges
c               and phi edges, respectively. [radians]
c
c - -  Output:  div(m,n) - real*8 array array of horizontal divergence of 
c               bt,bp evaluated at cell-centers (CE grid). [G/km or G/(km-sec)]
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
      real*8 :: rsun,dtheta,dphi
      real*8 :: sinth_hlf(m),sinth(m+1)
      real*8 :: bt(m+1,n),bp(m,n+1)
c
c - - output variable declarations:
c
      real*8 :: div(m,n)
c
c - - local variable declarations:
c
      integer :: mp1,np1,i,j,iph,jph,ip1,jp1
      real*8 :: oneodt,oneodp,rsuninv
c
      mp1=m+1
      np1=n+1
      oneodt=1.d0/dtheta
      oneodp=1.d0/dphi
      rsuninv=1.d0/rsun
c
      do jph=1,n
         j=jph
         jp1=j+1
         do iph=1,m
            i=iph
            ip1=i+1
            div(iph,jph)=oneodt*(bt(ip1,jph)*sinth(ip1)-bt(i,jph)
     1                  *sinth(i))/sinth_hlf(iph) +
     2                  oneodp*(bp(iph,jp1)-bp(iph,j))/sinth_hlf(iph)
         enddo
      enddo
c
c - - divide by rsun:
c
      div(:,:)=div(:,:)*rsuninv
c
      return
      end
