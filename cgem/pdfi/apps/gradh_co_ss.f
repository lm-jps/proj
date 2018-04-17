      subroutine gradh_co_ss(m,n,psi,rsun,sinth,dtheta,dphi,gradt,
     1 gradp)
c
c+
c - - Purpose:  Compute the horizontal gradient of the function psi, where
c               psi is assumed to be located on corners (COE grid)
c
c - - Usage:  call gradh_co_ss(m,n,psi,rsun,sinth,dtheta,dphi,gradt,gradp)
c
c - - compute horizontal gradient of the scalar array psi
c - - in spherical staggered coordinates.  This version computes the
c - - gradient assuming psi is at cell-corners.  The output components
c - - of the gradient is computed on cell edges.
c - - Input:  m,n - number of cell centers in theta, phi directions, resp.
c - - Input:  psi - array of the scalar potential at cell-corners, 
c - -         including external corners.  Dimensions m+1,n+1
c - - Input:  rsun - units for the radius of the Sun.
c - - Input:  sinth, sin(theta) computed over the co-latitude range,
c - -         at theta cell edges. Dimension m+1.
c - - Input:  dtheta,dphi - the spacing of cells in theta,phi directions
c - - Output: gradt - 2d array of size m,n+1 equal to gradient in
c - -         theta direction.  This is computed on phi edges (PE)
c - - Output: gradp - 2d array of size m+1,n equal to gradient
c - -         in phi direction.  This is computed on theta edges (TE).
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
      integer :: m,n,mp1,np1,iph,i,j,jph,jp1,ip1
      real*8 :: rsun,dtheta,dphi,oneodt,oneodp,rsuninv
      real*8 :: psi(m+1,n+1),gradt(m,n+1),gradp(m+1,n),sinth(m+1)
c
      mp1=m+1
      np1=n+1
c
      oneodt=1.0d0/dtheta
      oneodp=1.0d0/dphi
      rsuninv=1.d0/rsun
c
c - - gradp loop:
c
      do i=2,m
         do jph=1,n
            j=jph
            jp1=jph+1
            gradp(i,jph)=oneodp*(psi(i,jp1)-psi(i,j))/sinth(i)
         enddo
      enddo
c
c - - logic for handling north and south poles:
c
      if(sinth(1) .le. 1.d-15) then
c - - north pole
        gradp(1,1:n)=0.d0
      else
        i=1
        do jph=1,n
          j=jph
          jp1=jph+1
          gradp(i,jph)=oneodp*(psi(i,jp1)-psi(i,j))/sinth(i)
        enddo
      endif
      if(sinth(m+1) .le. 1.d-15) then
c - - south pole
        gradp(m+1,1:n)=0.d0
      else
        i=m+1
        do jph=1,n
          j=jph
          jp1=jph+1
          gradp(i,jph)=oneodp*(psi(i,jp1)-psi(i,j))/sinth(i)
        enddo
      endif
c
c - - gradt loop:
c
      do j=1,np1
         do iph=1,m
            i=iph
            ip1=iph+1
            gradt(iph,j)=oneodt*(psi(ip1,j)-psi(i,j))
         enddo
      enddo
c
c - - finally, divide all results by rsun
c
      gradt(:,:)=gradt(:,:)*rsuninv
      gradp(:,:)=gradp(:,:)*rsuninv
c
      return
      end
