      subroutine find_mask_ss(m,n,bmag0,bmag1,bmag2,bthr,mask)
c
c+
c     Purpose: calculate mask for three consecutive images.  Mask is 1 if
c              signal above threshold in all 3 images, otherwise 0
c
c - - Usage:   call find_mask(m,n,bmag0,bmag1,bmag2,bthr,mask)
c
c - - Input:   m,n - number of cell centers in theta, phi directions, resp.
c
c - - Input:   bmag0(m+1,n+1) - real*8 array of |B| evaluated on COE grid 
c              corners at t=t0 [G]
c
c - - Input:   bmag1(m+1,n+1) - real*8 array of |B| evaluated on COE grid 
c              corners at t=t1 [G]
c
c - - Input:   bmag2(m+1,n+1) - real*8 array of |B| evaluated on COE grid 
c              corners at t=t2 [G]
c
c - - Input:   bthr -  real*8 value of mask threshold [G]
c
c - - Output:  mask(m+1,n+1) - real*8 array of mask values evaluated on 
c              COE grid corners
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
c - - Written: MKD 25 Jan 2016
c     Corrected: MKD 21 Nov 2016 - corrected a bug: before mask Mij=0 only 
c     if all three Bmag0, Bmag1, Bmag2 < bthr in ij
c
c     now Mij=0 if any of three Bmag0, Bmag1, Bmag2 < bthr in ij
c
c
      implicit none
c
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: bmag0(m+1,n+1)
      real*8 :: bmag1(m+1,n+1)
      real*8 :: bmag2(m+1,n+1)
      real*8 :: bthr
c
c - - output variable declarations:
c
      real*8 :: mask(m+1,n+1)
c
c - - local variable declarations:
c
      real*8 :: mask0(m+1,n+1),mask1(m+1,n+1)
      real*8 :: mask2(m+1,n+1)
c
c - - initialize mask0,mask1,mask2 to 0:
c
      mask0(:,:)=0.d0
      mask1(:,:)=0.d0
      mask2(:,:)=0.d0
c
      where(bmag0 .ge. bthr) 
          mask0=1.d0
      endwhere
      
      where(bmag1 .ge. bthr) 
          mask1=1.d0
      endwhere

      where(bmag2 .ge. bthr) 
          mask2=1.d0
      endwhere
c     
      mask=mask0*mask1*mask2
c           
      return
      end
