      subroutine bhpottp2ll_ss(m,n,p,btpot,bppot,blonpot,blatpot)
c
c+
c - - Purpose: To transpose B_h potential-field arrays of the horizontal
c              components of the field from theta,phi,r order to lon,lat,r order
c              and flip sign to get B_lat.  Array values are assumed at
c              staggered Yee grid locations, and mid-way between radial
c              shells.
c
c - - Usage:   call bhpottp2ll_ss(m,n,p,btpot,bppot,blonpot,blatpot)
c
c - - Input:   m,n,p - number of cell centers in the theta (lat), phi (lon),
c              and radial directions, respectively.
c
c - - Input:   btpot (m+1,n,p),bppot(m,n+1,p) - real*8 arrays of the 
c              co-latitudinal and azimuthal components of the potential field 
c              evaluated at TE and PE locations (theta and phi edges, resp, 
c              midway between radial shells) [G]
c
c - - Output:  blonpot(n+1,m,p),blatpot(n,m+1,p) - real*8 arrays of 
c              longitudinal and latitudinal components of magnetic field, 
c              stored in lon,lat index order mid-way between radial shells,
c              PE and TE grid locations, respectively. [G]
c
c    NOTE:     Because the 3d arrays can be huge, this subroutine was written
c              in such a way that it will still work properly if btpot and
c              and blatpot actually occupy the same space in memory (similarly
c              with bppot and blonpot).  This can
c              be accomplished with equivalence statements in the calling
c              program, or more generally, by using both C and fortran pointers 
c              to define blatpot so that it points to the same memory as 
c              btpot, but with a different shape (and a similar relationship
c              between bppot and blonpot).  This can be done in 
c              the context of iso_c_binding in the fortran calling program.
c              More details can be found in the doc folder under the potential
c              field writeup.
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
      integer :: m,n,p
      real*8 :: btpot(m+1,n,p),bppot(m,n+1,p)
c
c - - output variables:
c
      real*8 :: blonpot(n+1,m,p),blatpot(n,m+1,p)
c
c - - local variables:
c
      integer :: q
      real*8 :: slicebt(m+1,n),slicebp(m,n+1)
      real*8 :: sliceblon(n+1,m),sliceblat(n,m+1)
c
      do q=1,p
         slicebt(:,:)=btpot(:,:,q)
         slicebp(:,:)=bppot(:,:,q)
         call bhyeetp2ll_ss(m,n,slicebt,slicebp,sliceblon,sliceblat)
         blonpot(:,:,q)=sliceblon(:,:)
         blatpot(:,:,q)=sliceblat(:,:)
      enddo
c
c - - we're done
c
      return
      end
