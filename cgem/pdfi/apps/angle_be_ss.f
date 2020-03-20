      subroutine angle_be_ss(m,n,et,ep,er,btcoe,bpcoe,brcoe,maskcoe,
     1           cosang)
c
c+
c - - Purpose: To compute the angle between the E and B vectors.
c
c - - Usage:  call angle_be_ss(m,n,et,ep,er,btcoe,bpcoe,brcoe,cosang)
c
c - - Input:  m,n - integer values for the number of cell centers in 
c             theta, phi directions, resp.
c
c - - Input:  et(m,n+1): real*8 array of theta component electric field 
c             [V/cm or G km/sec] on PE grid
c
c - - Input:  ep(m+1,n): real*8 array of phi component electric field 
c             [V/cm or G km/sec] on TE grid
c
c - - Input:  er(m+1,n+1): real*8 array of radial component of electric field 
c             [V/cm or G km/sec] on COE grid
c
c - - Input:  btcoe(m+1,n+1): real*8 array of theta component magnetic 
c             field [G] on COE grid
c
c - - Input:  bpcoe(m+1,n+1): real*8 array of phi component magnetic field [G] 
c             on COE grid
c
c - - Input:  brcoe(m+1,n+1): real*8 array of radial component magnetic field 
c             [G] on COE grid
c
c - - Input:  maskcoe(m+1,n+1): real*8 mask array for strong magnetic fields
c
c - - Output: cosang(m-1,n-1): real*8 array of the angle between E and B on 
c             CO grid expressed as cos(theta)
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
c - - Input variable declarations
c
      integer :: m,n
      real*8 :: et(m,n+1),ep(m+1,n),er(m+1,n+1)
      real*8 :: btcoe(m+1,n+1),bpcoe(m+1,n+1),brcoe(m+1,n+1)
      real*8 :: maskcoe(m+1,n+1)
c
c - - output variable declarations
c
      real*8 :: cosang(m-1,n-1)
c
c - - local variable declarations
c
      real*8 :: etc(m-1,n-1),epc(m-1,n-1),erc(m-1,n-1)
      real*8 :: btc(m-1,n-1),bpc(m-1,n-1),brc(m-1,n-1)
      real*8 :: maskco(m-1,n-1),bmag(m-1,n-1),emag(m-1,n-1)
c 
c - - interpolate e_horizontal to interior corners:
c
      call interp_eh_ss(m,n,et,ep,etc,epc)
c
      erc(1:m-1,1:n-1)=er(2:m,2:n)
      btc(1:m-1,1:n-1)=btcoe(2:m,2:n)
      bpc(1:m-1,1:n-1)=bpcoe(2:m,2:n)
      brc(1:m-1,1:n-1)=brcoe(2:m,2:n)
      maskco(1:m-1,1:n-1)=maskcoe(2:m,2:n)
c
      bmag(1:m-1,1:n-1)=sqrt(btc**2+bpc**2+brc**2)
      emag(1:m-1,1:n-1)=sqrt(etc**2+epc**2+erc**2)
c
      cosang(1:m-1,1:n-1)=0.d0
c
      where((maskco .eq. 1.d0) .and.(emag .gt. 0.))
        cosang=(etc*btc + epc*bpc + erc*brc)/(bmag*emag)
      endwhere
c
      return
      end
