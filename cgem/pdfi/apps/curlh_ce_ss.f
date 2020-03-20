      subroutine curlh_ce_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,
     1 dphi,curl)
c
c+     Purpose: Compute rhat dot curl of the vector with components et and ep.
c               et is on the PE grid, ep is on the TE grid.
c               Result is evaluated at cell centers (CE grid).
c
c - -  Usage:   call curlh_ce_ss(m,n,et,ep,rsun,sinth,sinth_hlf,
c               dtheta,dphi,curl)
c
c - -  Input:   m,n - integer numbers of cell-centers in theta, phi directions, 
c               resp.
c
c - -  Input:   et(m,n+1) - real*8 array of theta component of cE vector.
c               [G km/sec]
c
c - -  Input:   ep(m+1,n) - real*8 array of phi component of cE vector. 
c               [G km/sec]
c
c - -  Input:   rsun - real*8 value of radius of Sun. [km] Normally 6.96d5.
c
c - -  Input:   sinth(m+1) - real*8 array of theta cell-edge values of
c               sin(theta).
c
c - -  Input:   sinth_hlf(m) - real*8 array of theta-cell center values of
c               sin(theta).
c
c - -  Input:   dtheta - real*8 value of distance between theta edges. [radians]
c
c - -  Input:   dphi - real*8 value of distance between phi edges. [radians]
c
c - - Output:   curl(m,n) - real*8 array of cell-center values
c               of rhat dot curl(et,ep) [G/sec]
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
