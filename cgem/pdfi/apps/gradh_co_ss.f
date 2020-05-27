      subroutine gradh_co_ss(m,n,psi,rsun,sinth,dtheta,dphi,gradt,
     1 gradp)
c
c+
c - - Purpose:  Compute the horizontal gradient of the function psi, where
c               psi is assumed to be located on corners including edge
c               corners (COE grid).  The theta component gradt is on the PE
c               grid, and the phi component gradp is on the TE grid.
c
c - - Usage:    call gradh_co_ss(m,n,psi,rsun,sinth,dtheta,dphi,gradt,gradp)
c
c - - Input:    m,n - integer number of cell centers in theta, phi directions, 
c               resp.
c
c - - Input:    psi(m+1,n+1) - real*8 array of the scalar potential at 
c               cell-corners, including external corners. [G km^2/sec]
c
c - - Input:    rsun - real*8 value of radius of the Sun. [km]. Normally 6.96d5
c
c - - Input:    sinth(m+1) -  real*8 array of sin(colatitude) computed over 
c               the co-latitude range, at theta cell edges.
c
c - - Input:    dtheta,dphi - real*8 values of the spacing of cells in 
c               theta,phi directions. [radians]
c
c - - Output:   gradt(m,n+1) - real*8 array equal to gradient in
c               theta direction.  This is computed on phi edges (PE grid)
c               [G km/sec]
c - - Output:   gradp(m+1,n) - real*8 array equal to gradient
c               in phi direction.  This is computed on theta edges (TE grid).
c               [G km/sec]
c
c - - NOTE:     Plate Carree grid spacing (dtheta, dphi are constants) assumed.
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
c - - input variables:
c
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: sinth(m+1)
      real*8 :: psi(m+1,n+1)
c
c - - output variables:
c
      real*8 :: gradt(m,n+1),gradp(m+1,n)
c
c - - local variables:
c
      integer :: mp1,np1,iph,i,j,jph,jp1,ip1
      real*8 :: oneodt,oneodp,rsuninv
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
