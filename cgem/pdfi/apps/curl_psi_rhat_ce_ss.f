      subroutine curl_psi_rhat_ce_ss(m,n,psi,rsun,sinth_hlf,dtheta,
     1 dphi,curlt,curlp)
c
c+
c - - Purpose:  To compute the curl of psi * rhat unit vector, for psi located
c - - on the CE (cell center) grid.  Output is computed on cell edges.
c
c - - Usage: call curl_psi_rhat_ce_ss(m,n,psi,rsun,sinth_hlf,dtheta,
c    1 dphi,curlt,curlp)
c
c - - Input:  m,n - numbers of cell-centers in theta, phi directions, resp
c - - Input:  psi - array of the scalar potential at cell-centers, 
c - -         including ghost-zones.  Dimensions assume of size m+2,n+2
c - - Input:  rsun - units for the radius of the Sun.
c - - Input:  sinth_hlf, sin(theta) computed over the co-latitude range,
c - -         at cell centers, active zones only.  
c - - Input:  dtheta,dphi - the spacing of cells in theta,phi directions
c - - Output: curlt - 2d array of size m,n+1 equal to curl in
c - -         theta direction, evaluated on phi edges (PE).
c - - Output: curlp - 2d array of size m+1,n equal to curl
c - -         in phi direction, evaluated on theta edges (TE).
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
      real*8 :: psi(m+2,n+2),curlt(m,n+1),curlp(m+1,n),gradt(m+1,n),
     1          gradp(m,n+1),sinth_hlf(m)
c
      call gradh_ce_ss(m,n,psi,rsun,sinth_hlf,dtheta,dphi,gradt,gradp)
c
c - - Components of curl are related to components of gradient in a simple way:
c
      curlt(:,:)=gradp(:,:)
      curlp(:,:)=-gradt(:,:)
c
      return
      end
