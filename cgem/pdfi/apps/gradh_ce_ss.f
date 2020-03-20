      subroutine gradh_ce_ss(m,n,psi,rsun,sinth_hlf,dtheta,dphi,gradt,
     1 gradp)
c
c+
c - - Purpose: To compute the horizontal gradient of the function psi, where
c              psi is assumed to located at cell centers (CE grid 
c              plus ghost-zones). Theta component (gradt) located on theta
c              edges (TE grid), and phi component (gradp) located on phi edges
c              (PE grid).
c
c - - Usage:   call gradh_ce_ss(m,n,psi,rsun,sinth_hlf,dtheta,dphi,gradt,gradp)
c
c - - Input:   m,n - integer number of cell centers in theta, phi directions, 
c              respectively.
c
c - - Input:   psi(m+2,n+2) - real*8 array of the scalar potential at 
c              cell-centers (CE grid), plus ghost-zones. [G km^2/sec]
c
c - - Input:   rsun - real*8 value of radius of the Sun. [km].  Normally 6.96d5
c
c - - Input:   sinth_hlf(m) - real*8 array of sin(colatitude) computed over 
c              the co-latitude range, at cell centers.  
c
c - - Input:   dtheta,dphi - real*8 values of the spacing of cells in 
c              theta,phi directions. [radians]
c
c - - Output:  gradt(m+1,n) - real*8 array equal to gradient in
c              theta direction.  The output is at TE grid locations. [G km/sec]
c
c - - Output:  gradp(m,n+1) - real*8 array equal to gradient in phi direction.  
c              The output is at PE grid locations. [G km/sec]
c
c - - NOTE:    Plate Carree grid spacing (dtheta, dphi are constants) assumed
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
c - - input variables:
c
      integer :: m,n
      real*8 :: psi(m+2,n+2)
      real*8 :: sinth_hlf(m)
      real*8 :: rsun,dtheta,dphi
c
c - - output variables:
c
      real*8 :: gradt(m+1,n),gradp(m,n+1)
c
c - - local variables:
c
      real*8 :: oneodt,oneodp,rsuninv
      integer :: mp1,np1,iph,imh,i,j,jph,jmh
c
      oneodt=1.0d0/dtheta
      oneodp=1.0d0/dphi
      rsuninv=1.d0/rsun
      np1=n+1
      mp1=m+1
c
c - - gradp loop:
c
      do iph=1,m
         do j=1,np1
            jmh=j
            jph=j+1
            gradp(iph,j)=oneodp*(psi(iph+1,jph)-psi(iph+1,jmh)) /
     1      sinth_hlf(iph)
         enddo
      enddo
c
c - - gradt loop:
c
      do jph=1,n
         do i=1,mp1
            imh=i
            iph=i+1
            gradt(i,jph)=oneodt*(psi(iph,jph+1)-psi(imh,jph+1))
         enddo
      enddo
c 
c - - Finally, divide all results by rsun
c
      gradt(:,:)=gradt(:,:)*rsuninv
      gradp(:,:)=gradp(:,:)*rsuninv
c
      return
      end
