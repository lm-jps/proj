       subroutine get_pils_ss(m,n,bmap,pilmap,thresh,dilation_param)
c
c+
c - -  Purpose: Given input "bmap" array (m-1,n-1) corresponding to magnetogram 
c               pixels, this returns "pilmap" an (m-1, n-1) bitmap (of integer 
c               type, with default precision) for all pixels closer to 
c               opposite polarity than dilation_param (integer) 
c - -  Usage:   call get_pils_ss(m,n,bmap,pilmap,thresh,dilation_param)
c - -  Input:   m,n - integer no. of cell centers in colat, lon, resp.
c - -  Input:   bmap(m-1,n-1) - real*8 magnetogram array
c - -  Output:  pilmap(m-1,n-1) - integer array of PIL locations, equal to 1
c               within PILs, 0 outside of PILs.
c - -  Input:   thresh - real*8 threshold (abs. value) for defining magnetic
c               regions.
c - -  Input:   dilation_param - integer value defining width of PILs.  Width
c               is 2*dilation_param + 1.
c-
c
c
c HISTORY: 2016/10/25, BT Welsch: started
c          2017/08/01 modified by GHF for compatibility with pdfi_ss library
c
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
      implicit none
c
c - - calling argument declarations (input):
c
      integer :: m,n
c
c - - max. dist. btwn. +/- to be a PIL
c
      integer :: dilation_param 
c
c  - - INPUT magnetogram, B_los (or B_r) 
c
      real*8 :: bmap(m-1,n-1)
c
c - - threshold in abs(Bmap) below which pixel values are ignored
c
      real*8 :: thresh 
c
c - - calling argument declarations (output):
c
c - - map of PIL pixels (2D integer array)
c
      integer :: pilmap(m-1,n-1) 
c
c - - local variable declarations:
c
      integer :: posmap(m-1,n-1), negmap(m-1,n-1) 
      integer :: dilpos(m-1,n-1),dilneg(m-1,n-1) 
c
c - - start doing the actual work:
c
      where (bmap .ge. thresh) 
         posmap = 1
      elsewhere
         posmap = 0
      end where
c
      where (bmap .le. -thresh) 
        negmap = 1
      elsewhere
        negmap = 0
      end where
c
      call dilate_ss(m,n,posmap,dilpos,dilation_param)
      call dilate_ss(m,n,negmap,dilneg,dilation_param)
c
      pilmap = dilpos*dilneg
c
c - - we're done
c
      return
      end
