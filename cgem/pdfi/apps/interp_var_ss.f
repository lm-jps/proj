      subroutine interp_var_ss(m,n,btcoe,bpcoe,brcoe,
     1 a,b,bt,bp,br)
c+
c - -   Purpose: To interpolate a 3-d vector from COE grid to its proper 
c                locations on the TE, PE, and CE staggered grids.
c                Interpolation is done using straightforward linear
c                averages, motivated by simplicity and noise reduction.
c
c - -   Usage:   call interp_data_ss(m,n,btcoe,bpcoe,brcoe,a,b,
c                bt,bp,br)
c
c - -    Input:  m,n - integer number of cell centers in colat, lon directions,
c                respectively.
c
c - -    Input:  btcoe(m+1,n+1) - real*8 array of theta component of B located
c                on the COE grid (corners, including boundary corners) [G]
c
c - -    Input:  bpcoe(m+1,n+1) - real*8 array of phi component of B located
c                on the COE grid [G]
c
c - -    Input:  brcoe(m+1,n+1) - real*8 array of the radial component of B
c                located on the COE grid [G] 
c
c - -    Input:  a,b - real*8 values of the mininum and maximum colatitude on 
c                the COE grid [radians]
c
c - -   Output:  bt(m+1,n) - real*8 array of theta component of B located
c                at theta edge locations (TE grid) [G]
c
c - -   Output:  bp(m,n+1) - real*8 array of phi component of B located at 
c                phi edge locations (PE grid) [G]
c
c - -   Output:  br(m,n) - real*8 array of the radial component of B located
c                at cell center locations (CE grid) [G]
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
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: btcoe(m+1,n+1),bpcoe(m+1,n+1),brcoe(m+1,n+1),a,b
c
c - - output variable declarations:
c
      real*8 :: bt(m+1,n),bp(m,n+1),br(m,n)
c
c - - local variable declarations:
c
      integer :: i,j,jp1,jph,ip1,iph
      real*8 sinth(m+1),sinth_hlf(m)
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - interpolate br at cell centers from edges+corners using flux wtd avg
c
      do i=1,m
         iph=i
         ip1=i+1
         do j=1,n
            jph=j
            jp1=j+1
            br(iph,jph)=(sinth(i)*0.5d0*(brcoe(i,j)+brcoe(i,jp1))
     1              +sinth(ip1)*0.5d0*(brcoe(ip1,j)+brcoe(ip1,jp1)))
     2              /(sinth(i)+sinth(ip1))
         enddo
      enddo
c
c - - interpolate bt at theta edges:
c
      do i=1,m+1
         do j=1,n
            jph=j
            jp1=j+1
            bt(i,jph)=0.5d0*(btcoe(i,j)+btcoe(i,jp1))
         enddo
      enddo
c
c - - interpolate bp at phi edges:
c
      do i=1,m
         iph=i
         ip1=i+1
         do j=1,n+1
            bp(iph,j)=0.5d0*(bpcoe(i,j)+bpcoe(ip1,j))
         enddo
      enddo
c
      return
      end
