      subroutine gradh_ce_ss(m,n,psi,rsun,sinth_hlf,dtheta,dphi,gradt,
     1 gradp)
c
c+
c - - Purpose: To compute the horizontal gradient of the function psi, where
c     psi is assumed to located at cell centers (plus ghost-zones).
c
c - - Usage:  call gradh_ce_ss(m,n,psi,rsun,sinth_hlf,dtheta,dphi,gradt,gradp)
c
c - - compute horizontal gradient of the scalar array psi
c - - in spherical staggered coordinates.  This version computes the
c - - gradient, assuming psi is at cell-centers.  The gradient itself
c - - is computed on edges, with gradt on theta edges (TE), while gradp
c - - is computed on phi edges (PE).
c
c - - Input:  m,n - number of cell centers in theta, phi directions, resp.
c - - Input:  psi - array of the scalar potential at cell-centers, 
c             including ghost-zones.  Dimensions assume of size m+2,n+2
c - - Input:  rsun - units for the radius of the Sun.
c - - Input:  sinth_hlf, sin(theta) computed over the co-latitude range,
c             at cell centers, active zones only.  
c - - Input:  dtheta,dphi - the spacing of cells in theta,phi directions
c - - Output: gradt - 2d array of size m+1,n equal to gradient in
c             theta direction.  The output is at TE locations.
c - - Output: gradp - 2d array of size m,n+1, equal to gradient
c             in phi direction.  The output is at PE locations.
c - - NOTE:  Plate Carree grid spacing (dtheta, dphi are constants) assumed!
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
      integer :: m,n,mp1,np1,iph,imh,i,j,jph,jmh
      real*8 :: psi(m+2,n+2),gradt(m+1,n),gradp(m,n+1),sinth_hlf(m)
      real*8 :: rsun,dtheta,dphi,oneodt,oneodp,rsuninv
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
