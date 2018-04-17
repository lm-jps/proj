      subroutine find_mask_ss(m,n,bmag0,bmag1,bmag2,bthr,mask)
c
c+
c     Purpose: calculate mask for three consecutive images.  Mask is 1 if
c              signal above threshold in all 3 images, otherwise 0
c
c - - Usage:  call find_mask(m,n,bmag0,bmag1,bmag2,bthr,mask)
c - - Input:  m,n - number of cell centers in theta, phi directions, resp.
c - - Input:  bmag0(m+1,n+1) - |B| evaluated on interior cell corners at t=t0
c - - Input:  bmag1(m+1,n+1) - |B| evaluated on interior cell corners at t=t1
c - - Input:  bmag2(m+1,n+1) - |B| evaluated on interior cell corners at t=t2
c - - Input:  bthr -  mask threshold [Gauss]
c
c - - Output: mask(m+1,n+1) - mask evaluated on interior cell corners
c
c-
c - - Written: MKD 25 Jan 2016
c     Corrected: MKD 21 Nov 2016 - corrected a bug: before mask Mij=0 only 
c     if all three Bmag0, Bmag1, Bmag2 < bthr in ij
c
c     now Mij=0 if any of three Bmag0, Bmag1, Bmag2 < bthr in ij
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
      integer :: m,n
      real*8 :: bmag0(m+1,n+1)
      real*8 :: bmag1(m+1,n+1)
      real*8 :: bmag2(m+1,n+1)
      real*8 :: mask(m+1,n+1)
      real*8 :: bthr
      real*8 :: mask0(m+1,n+1),mask1(m+1,n+1)
      real*8 :: mask2(m+1,n+1)
c
c - - initialize mask0,mask1,mask2 to 0:
c
      mask0(:,:)=0.d0
      mask1(:,:)=0.d0
      mask2(:,:)=0.d0


      where(bmag0 .ge. bthr) 
          mask0=1.d0
      endwhere
      
      where(bmag1 .ge. bthr) 
          mask1=1.d0
      endwhere

      where(bmag2 .ge. bthr) 
          mask2=1.d0
      endwhere
      
      mask=mask0*mask1*mask2
            
      return
      end
