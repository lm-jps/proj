      subroutine ahpottp2ll_ss(m,n,p,atpepot,aptepot,alonpot,alatpot)
c
c+
c   Purpose:  To transpose potential-field 3D vector potential A_h data arrays 
c             from theta,phi,r order to lon,lat,r order, and flipping the sign 
c             to get alatpot.
c
c     Usage:  call ahpottp2ll_ss(m,n,p,atpepot,aptepot,alonpot,alatpot)
c
c     Input:  m,n,p - number of cell centers in the theta (lat), and phi (lon)
c             and radial (r) directions, respectively.
c
c     Input:  atpepot(m,n+1,p+1),aptepot(m+1,n,p+1) - real*8 arrays of the 
c             co-latitudinal and azimuthal components of the vector potential 
c             evaluated at PE grid and TE grid locations 
c             (phi and theta edges, resp.).  Values are on radial shells.
c             [G-km]
c
c    Output:  alonpot(n,m+1,p+1),alatpot(n+1,m,p+1) - real*8 arrays of 
c             longitudinal and latitudinal components of potential-field, 
c             vector potential, stored in lon,lat index order, TE and PE grid 
c             locations, respectively.  Values are on radial shells.  [G-km]
c
c    NOTE:    Because the 3d arrays can be huge, this subroutine was written
c             in such a way that it will still work properly if atpepot and
c             and alatpot actually occupy the same space in memory (similarly
c             with aptepot and alonpot).  This can
c             be accomplished with an equivalence statements in the calling
c             program, or more generally, by using both C and fortran pointers 
c             to define alatpot so that it points to the same memory as 
c             atpepot, but with a different shape (and a similar relationship
c             between aptepot and alonpot).  This can be done in 
c             the context of iso_c_binding in the fortran calling program.
c             More details can be found in the doc folder under the potential
c             field writeup.
c-
c   PDFI_SS Electric Field Inversion Software
c   http://cgem.ssl.berkeley.edu/cgi-bin/cgem/PDFI_SS/index
c   Copyright (C) 2015-2018 University of California
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
      integer :: m,n,p
      real*8 :: atpepot(m,n+1,p+1),aptepot(m+1,n,p+1)
c
c - - output variables:
c
      real*8 :: alonpot(n,m+1,p+1),alatpot(n+1,m,p+1)
c
c - - local variables:
c
      integer :: q
      real*8 :: sliceat(m,n+1),sliceap(m+1,n),slicealon(n,m+1),
     1 slicealat(n+1,m)
c
      do q=1,p+1
         sliceat(:,:)=atpepot(:,:,q)
         sliceap(:,:)=aptepot(:,:,q)
c
c - - use ehyeetp2ll_ss since vector potential and electric field on same grid
c
         call ehyeetp2ll_ss(m,n,sliceat,sliceap,slicealon,slicealat)
         alonpot(:,:,q)=slicealon(:,:)
         alatpot(:,:,q)=slicealat(:,:)
      enddo
c
c - - we're done
c
      return
      end
