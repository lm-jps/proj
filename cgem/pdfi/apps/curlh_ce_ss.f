      subroutine curlh_ce_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,
     1 dphi,curl)
c
c+ staggered spherical grid version of curlh: evaluated at cell-centers
c  Purpose: Compute rhat cdot curl of the vector with components et and ep.
c           Result is evaluated at cell centers (CE grid).
c
c  Usage:  call curlh_ce_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,dphi,curl)
c  Input:  m,n - numbers of cell-centers in theta, phi directions, resp.
c  Input:  et - double prec. theta component of vector, dimensioned m,n+1
c  Input:  ep - double prec. phi component of vector, dimensioned m+1,n
c  Input:  rsun - double prec., assumed units of radius of Sun
c  Input:  sinth - double prec. array of theta cell-edge values of sin(theta),
c          must be dimensioned m+1
c  Input:  sinth_hlf - double prec. array of theta-cell center values of
c          sin(theta), must be dimensioned m
c  Input:  dtheta - double prec. value of distance between theta edges
c  Input:  dphi - double prec. value of distance between phi edges
c Output:  curl - double prec. array dimensioned m,n - cell-center values
c          of rhat dot curl(et,ep)
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
      integer :: m,n,i,iph,ip1,j,jph,jp1
      real*8 :: rsun,dtheta,dphi,rsuninv,oneodp,oneodt
      real*8 :: sinth(m+1),sinth_hlf(m),et(m,n+1),ep(m+1,n),curl(m,n)
c
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
            curl(iph,jph)=oneodt*(ep(ip1,jph)*sinth(ip1)-ep(i,jph)
     1      *sinth(i))/sinth_hlf(iph) -
     2      oneodp*(et(iph,jp1)-et(iph,j))/sinth_hlf(iph)
         enddo
      enddo
c
c - - divide by rsun:
c
      curl(:,:)=curl(:,:)*rsuninv
c
      return
      end
