      subroutine curl_psi_rhat_co_ss(m,n,psi,rsun,sinth,dtheta,dphi,
     1 curlt, curlp)
c
c+
c - - Purpose:  To compute the curl of psi * rhat unit vector, with
c - - the psi variable located on the CO grid (interior corners).
c - - Output computed at cell edges.
c
c - - Usage:  call curl_psi_rhat_co_ss(m,n,psi,rsun,sinth,dtheta,
c    1  dphi,curlt,curlp)
c
c - - Input:  m,n - number of cell centers in theta, phi directions, resp.
c - - Input:  psi - array of the scalar potential at cell-corners, 
c - -         including ghost-zones.  Dimensions assume of size m+1,n+1
c - - Input:  rsun - units for the radius of the Sun.
c - - Input:  sinth, sin(theta) computed over the co-latitude range,
c - -         at theta cell edges.  Dimension m+1.
c - - Input:  dtheta,dphi - the spacing of cells in theta,phi directions
c - - Output: curlt - 2d array of size m+1,n equal to curl in
c - -         theta direction, evaluated at theta edges (TE).
c - - Output: curlp - 2d array of size m,n+1, equal to curl
c - -         in phi direction, evaluated at phi edges (PE)
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
      integer :: m,n
      real*8 :: rsun,dtheta,dphi
      real*8 :: psi(m+1,n+1),sinth(m+1),curlt(m+1,n),curlp(m,n+1)
      real*8 :: gradt(m,n+1),gradp(m+1,n)
c
      call gradh_co_ss(m,n,psi,rsun,sinth,dtheta,dphi,gradt,gradp)
c
c - - components of curl are related to the gradient in a simple way:
c
      curlt(:,:)=gradp(:,:)
      curlp(:,:)=-gradt(:,:)
      return
      end
