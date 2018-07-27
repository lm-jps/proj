      subroutine fix_mask_ss(mdim,ndim,flag,mask)
c
c+
c - - Purpose: For mask values between 0 and 1, sets them to 0 or 1
c              depending on the value of flag:  If flag=1, then intermediate
c              values set to 1; if flag=0, then intermediate values set to 0
c
c - - Usage:   call fix_mask(mdim,ndim,flag,mask)
c
c - - Input:   mdim,ndim - number of mask values in theta, phi directions, resp.
c
c - - Input:   flag - integer equal to 0 or 1, which determines whether
c              intermediate mask values are set to 0 or 1
c
c - - Input/Output: 
c              mask(mdim,ndim) - real*8 2d mask array
c
c-
c
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
c
      implicit none
c
c - - calling argument declarations:
c
      integer :: mdim,ndim,flag
      real*8 :: mask(mdim,ndim)
c
      if ((flag .ne. 0) .and. (flag .ne. 1)) then
         write(6,*) 'fix_mask_ss:  flag = ',flag,' is an illegal value'
         stop
      endif
      where((mask .gt. 0.d0) .and. (mask .lt. 1.d0))
          mask=dble(flag)
      endwhere
c
      return
      end
