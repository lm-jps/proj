      subroutine divh_co_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,dphi,
     1 div)
c
c+ 
c - -  Purpose: Compute horizontal divergence of the vector cE with theta,phi 
c               components et, (PE grid) and ep, (TE grid) evaluated at 
c               interior cell corners (CO grid).
c
c - -  Usage:   call divh_co_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,dphi,div)
c
c - -  Input:   m,n - integer number of cell centers in theta and phi 
c               directions, resp.
c
c - -  Input:   et(m,n+1) - real*8 array of theta component of cE, defined 
c               on phi edges (PE grid) [G km/sec]
c
c - -  Input:   ep(m+1,n) - real*8 array of phi component of cE, defined
c               on theta edges (TE grid) [G km/sec]
c
c - -  Input:   rsun - real*8 value of radius of Sun [km]. Normally 6.96d5
c
c - -  Input:   sinth(m+1) - real*8 array of sin(theta) evaluated at 
c               theta edges
c
c - -  Input:   sinth_hlf(m) - real*8 array of sin(theta) evaluated at
c               cell centers.
c
c - -  Input:   dtheta,dphi - real*8 values of angular distance between theta 
c               and phi edges, respectively. [radians]
c
c - - Output:   div(m-1,n-1) - real*8 array of horizontal div. of et,ep 
c               evaluated at interior cell corners (CO grid). [G/sec]
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
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: sinth(m+1),sinth_hlf(m),et(m,n+1),ep(m+1,n)
c
c - - output variable declarations:
c
      real*8 :: div(m-1,n-1)
c
c - - local variable declarations:
c
      integer :: nm1,mm1,imh,iph,jmh,jph,i,j
      real*8 :: rsuninv,oneodt,oneodp
c
      mm1=m-1
      nm1=n-1
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
            div(i,j)=oneodt*(et(iph,j+1)*sinth_hlf(iph)-
     1               et(imh,j+1)*sinth_hlf(imh))/sinth(i+1) +
     2               oneodp*(ep(i+1,jph)-ep(i+1,jmh))/sinth(i+1)
         enddo
      enddo
c
c - - divide by rsun:
c
      div(:,:)=div(:,:)*rsuninv
c
      return
      end
