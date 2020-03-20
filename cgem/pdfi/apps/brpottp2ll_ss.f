      subroutine brpottp2ll_ss(m,n,p,brpottp,brpotll)
c
c+
c - - Purpose: To transpose B_r 3D array of potential field radial components
c              from theta,phi,r order to lon,lat,r order.
c
c - - Usage:   call brpottp2ll_ss(m,n,p,brpottp,brpotll)
c
c - - Input:   m,n,p - number of cell centers in the theta (lat), and phi (lon)
c              and radial (r) directions, respectively.
c
c - - Input:   brpottp(m,n,p+1) - real*8 array of the radial component of the
c              magnetic field evaluated at CE grid locations in theta,phi order
c              (cell-centers), for each radial shell. [G]
c
c - - Output:  brpotll(n,m,p+1) - real*8 array of the radial magnetic field 
c              component at CE grid locations stored in lon,lat,r index order. 
c              [G]
c
c - - NOTE:    Because the 3d arrays can be huge, this subroutine was written
c              in such a way that it will still work properly if brpottp
c              and brpotll actually occupy the same space in memory.  This can
c              be accomplished with an equivalence statement in the calling
c              program, or more generally, using both C and fortran pointers 
c              to define brpotll so that it points to the same memory as 
c              brpottp, but with a different shape.  This can be done in 
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
      real*8 :: brpottp(m,n,p+1)
c
c - - output variables:
c
      real*8 :: brpotll(n,m,p+1)
c
c - - local variables:
c
      integer :: q
      real*8 :: slicetp(m,n),slicell(n,m)
c
      do q=1,p+1
         slicetp(:,:)=brpottp(:,:,q)
         call bryeetp2ll_ss(m,n,slicetp,slicell)
         brpotll(:,:,q)=slicell(:,:)
      enddo
c
c - - we're done
c
      return
      end
