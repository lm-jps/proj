      subroutine gradh_sc(m,n,psi,rsun,sinth_gh,dtheta,dphi,gradt,
     1 gradp)
c
c+
c - - Purpose:  Compute gradient of psi, using centered grid formalism
c
c               NOTE: *DO NOT USE* for staggered grid variables in PDFI_SS.  
c               This subroutine *only* for use within the relax_psi_3d_ss 
c               subroutine, which assumes a centered grid formalism.
c
c - - Usage:    call gradh_sc(m,n,psi,rsun,sinth_gh,dtheta,dphi,gradt,gradp)
c
c - - Input:    m,n - number of cell centers in theta, phi directions, resp.
c
c               NOTE:  m,n are smaller by two than the values used for the
c               staggered grid description in PDFI_SS.
c
c - - Input:    psi(m+3,n+3) - real*8 array of the electric scalar potential, 
c               including ghost zones. [G km^2/sec]
c
c - - Input:    rsun - real*8 value for radius of the Sun [km]. Normally 6.96d5
c
c - - Input:    sinth_gh(m+3) - real*8 array of sin(colatitude) computed 
c               over the co-latitude range, with ghostzones included.
c
c - - Input:    dtheta,dphi - real*8 values of the spacing of cells in 
c               theta,phi directions [radians]
c
c - - Output:   gradt(m+1,n+1) - real*8 array equal to gradient of psi in theta 
c               direction. [G km/sec]
c
c - - Output:   gradp(m+1,n+1) - real*8 array equal to gradient of psi in phi 
c               direction. [G km/sec]
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
      integer :: m,n
c
c - - input variable declarations:
c
      real*8 :: rsun,dtheta,dphi
      real*8 :: psi(m+3,n+3),sinth_gh(m+3)
c
c - - output variable declarations:
c
      real*8 :: gradt(m+1,n+1),gradp(m+1,n+1)
c
c - - local variable declarations:
c
      real*8 :: gradttmp(m+3,n+3),gradptmp(m+3,n+3)
      integer :: ntmax,npmax,j
      real*8 :: halfodt,halfodp,rsuninv
c
      ntmax=m+3
      npmax=n+3
      rsuninv=1.d0/rsun
      halfodt=0.5d0/dtheta
      halfodp=0.5d0/dphi
c
      do j=2,npmax-1
         gradttmp(2:ntmax-1,j)=halfodt*(psi(3:ntmax,j)-psi(1:ntmax-2,j))
         gradptmp(2:ntmax-1,j)=halfodp*(psi(2:ntmax-1,j+1)-
     1   psi(2:ntmax-1,j-1))/sinth_gh(2:ntmax-1)
      enddo
c
      gradt(1:m+1,1:n+1)=gradttmp(2:m+2,2:n+2)*rsuninv
      gradp(1:m+1,1:n+1)=gradptmp(2:m+2,2:n+2)*rsuninv
c
      return
      end
