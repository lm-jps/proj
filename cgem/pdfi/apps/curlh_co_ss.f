      subroutine curlh_co_ss(m,n,bt,bp,rsun,sinth,sinth_hlf,dtheta,
     1 dphi,curl)
c
c+     Purpose: Compute rhat cdot curl of vector with components bt,bp.  
c               The bt array is on the TE grid, bp array is on the PE grid, and
c               the result is evaluated on the CO grid (interior corners).
c
c - -  Usage:   call curlh_co_ss(m,n,bt,bp,rsun,sinth,sinth_hlf,
c               dtheta,dphi,curl)
c - -  Input:   m,n - integer number of cell centers in theta, phi directions, 
c               respectively.
c - -  Input:   bt(m+1,n) - real*8 array of theta component of vector 
c               [G or G/sec].
c - -  Input:   bp(m,n+1) - real*8 array of phi component of vector 
c               [G or G/sec].
c - -  Input:   rsun - real*8 value of radius of Sun [km].  Normally 6.96d5.
c - -  Input:   sinth(m+1) - real*8 array of sin(theta) evaluated at 
c               theta edges.
c - -  Input:   sinth_hlf(m) - real*8 array of sin(theta) evaluated 
c               cell centers.
c - -  Input:   dtheta - real*8 value of angular distance between 
c               theta edges [radians]
c - -  Input:   dphi - real*8 value of angular distance between 
c               phi edges [radians]
c - - Output:   curl(m-1,n-1) - real*8 array containing the
c               curl evaluated at interior cell corners (CO grid). 
c               [G/km or G/(km-sec)]
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
      integer :: m,n,mm1,nm1,j,jmh,jph,i,imh,iph
      real*8 :: rsun,dtheta,dphi,rsuninv,oneodt,oneodp
      real*8 :: sinth(m+1),sinth_hlf(m),bt(m+1,n),bp(m,n+1),
     1          curl(m-1,n-1)
c
      nm1=n-1
      mm1=m-1
c
      oneodt=1.d0/dtheta
      oneodp=1.d0/dphi
      rsuninv=1.d0/rsun
c
      do j=1,nm1
         jmh=j
         jph=j+1
         do i=1,mm1
            imh=i
            iph=i+1
            curl(i,j)=oneodt*(bp(iph,j+1)*sinth_hlf(iph) -
     1      bp(imh,j+1)*sinth_hlf(imh))/sinth(i+1) -
     2      oneodp*(bt(i+1,jph)-bt(i+1,jmh))/sinth(i+1)
         enddo
      enddo
c
c - - divide by rsun:
c
      curl(:,:)=curl(:,:)*rsuninv
c
      return
      end
