      subroutine divh_co_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,dphi,
     1 div)
c
c+ spherical staggered grid version of divh evaluated on interior cell corners
c  Purpose: Compute horizontal divergence of the vector with components et,ep
c           evaluated at interior cell corners (CO grid).
c
c  Usage:  call divh_co_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,dphi,div)
c  Input:  m,n - number of cell centers in theta and phi directions, resp.
c  Input:  et - double prec. array of theta component of vector, dims m,n+1
c  Input:  ep - double prec. array of phi component of vector, dims m+1,n
c  Input:  rsun - double prec., assumed value of radius of Sun
c  Input:  sinth - double prec. array of sin(theta) evaluated at theta edges,
c          dimension must be m+1
c  Input:  sinth_hlf - double prec. array of sin(theta) evaluated cell centers.
c          dimension must be m.
c  Input:  dtheta - double prec. value of angular distance between theta edges
c  Input:  dphi - double prec. value of angular distance between phi edges
c Output:  div - double prec. array of horizontal div. of et,ep evaluated
c          at active cell corners.  Dimensions m-1,n-1.
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
      integer :: m,n,nm1,mm1,imh,iph,jmh,jph,i,j
      real*8 :: rsun,dtheta,dphi,rsuninv,oneodt,oneodp
      real*8 :: sinth(m+1),sinth_hlf(m),et(m,n+1),ep(m+1,n),
     1          div(m-1,n-1)
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
