      subroutine v_ideal_ss(m,n,btcoe,bpcoe,brcoe,et,ep,er,maskcoe,
     1 vtcoe,vpcoe,vrcoe)
c+
c - - Purpose:  To compute the velocity field v given E and B such that
c               cE is as close as possible to -v cross B.
c
c - - Method:   v = cE cross B / B^2 (where B is above threshold),
c               =0 otherwise.
c
c - - Usage:    call v_ideal_ss(m,n,btcoe,bpcoe,brcoe,et,ep,er,maskcoe,
c               vtco,vpco,vrco)
c
c - - Input:    m,n - integer values of the number of cell centers in 
c               co-latitude and longitude directions, respectively.
c
c - - Input:    btcoe(m+1,n+1),bpcoe(m+1,n+1),brcoe(m+1,n+1) - real*8 arrays
c               of the theta, phi, and radial components of the magnetic field
c               evaluated on the COE grid. [G]
c
c - - Input:    et(m,n+1),ep(m+1,n),er(m+1,n+1) - real*8 arrays of the theta,
c               phi, and radial components of electic field multiplied by c, 
c               on the PE, TE, and COE grids, respectively. [G km/s]
c
c - - Input:    maskcoe(m+1,n+1) - mask array on COE grid, set to 1 where the
c               field strength is above threshold, and to 0 otherwise.
c
c - - Output:   vtcoe(m+1,n+1) - real*8 array of the theta component of v. 
c               [km/s]
c
c - - Output:   vpcoe(m+1,n+1) - real*8 array of the phi component of v. [km/s]
c
c - - Output:   vrcoe(m+1,n+1) - real*8 array of the radial component of v. 
c               [km/s]
c
c - - Note:     vtcoe,vpcoe,vrcoe will almost certainly be 0 on the edges, 
c               because ercoe is 0 there, and et, ep are not interpolated
c               to the edges.
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
c - - Input variables:
c
      integer :: m,n
      real*8 :: btcoe(m+1,n+1),bpcoe(m+1,n+1),brcoe(m+1,n+1)
      real*8 :: maskcoe(m+1,n+1)
      real*8 :: et(m,n+1),ep(m+1,n),er(m+1,n+1)
c
c - - Output variables:
c
      real*8 :: vtcoe(m+1,n+1),vpcoe(m+1,n+1),vrcoe(m+1,n+1)
c
c - - Local variables:
c
      real*8 :: etco(m-1,n-1),epco(m-1,n-1),etcoe(m+1,n+1),
     1 epcoe(m+1,n+1),bmag2(m+1,n+1),bm2inv(m+1,n+1)
c
      etcoe(:,:)=0.d0
      epcoe(:,:)=0.d0
c
      call interp_eh_ss(m,n,et,ep,etco,epco)
c
      etcoe(2:m,2:n)=etco(1:m-1,1:n-1)
      epcoe(2:m,2:n)=epco(1:m-1,1:n-1)
c
      bmag2 = btcoe**2+bpcoe**2+brcoe**2
c
      where(maskcoe .gt. 0.d0)
        bm2inv=1.d0/bmag2
      elsewhere
        bm2inv=0.d0
      endwhere
c
c - - take cross product to get all 3 components of (E X B)/B^2:
c
      where(maskcoe .gt. 0.d0)
        vtcoe=(epcoe*brcoe-er*bpcoe)*bm2inv
        vpcoe=(er*btcoe-etcoe*brcoe)*bm2inv
        vrcoe=(etcoe*bpcoe-epcoe*btcoe)*bm2inv
      elsewhere
        vtcoe=0.d0
        vpcoe=0.d0
        vrcoe=0.d0
      endwhere
c
c - - we're done
c
      return
      end
