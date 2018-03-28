      subroutine stack_3d_ll(m,n,bloncoe0,blatcoe0,brllcoe0,bloncoe1,
     1 blatcoe1,brllcoe1,vloncoe0,vlatcoe0,vlosllcoe0,vloncoe1,
     2 vlatcoe1,vlosllcoe1,lloncoe0,llatcoe0,lrllcoe0,lloncoe1,
     3 llatcoe1, lrllcoe1, data_ll_3d)
c
c+
c - - Purpose:  Take the 18 2D input arrays needed for pdfi_wrapper4jsoc
c               and stack them up into a 3d array, with the 3rd dimension
c               equal to 18.  Arrays are all in lon,lat order on COE grid.
c
c - - Input:    m,n - integer number of cell-centers in latitude 
c               (or co-latitude) and longitude directions, respectively.
c - - Input:    bloncoe0,blatcoe0,brllcoe0 - real*8(n+1,m+1) arrays of magnetic
c               field components (lon,lat,radial) at time 0
c - - Input:    bloncoe1,blatcoe1,brllcoe1 - real*8(n+1,m+1) arrays of magnetic
c               field components (lon,lat,radial) at time 1
c - - Input:    vloncoe0,vlatcoe0 - real*8(n+1,m+1) arrays of horizontal vel.
c               components (lon,lat) computed from FLCT at time 0
c - - Input:    vlosllcoe0 - real*8(n+1,m+1) array of Doppler shift in m/s, 
c               positive sign is redshifted, at time 0.
c - - Input:    vloncoe1,vlatcoe1 - real*8(n+1,m+1) arrays of horizontal vel.
c               components (lon,lat) computed from FLCT at time 1
c - - Input:    vlosllcoe1 - real*8(n+1,m+1) array of Doppler shift in m/s, 
c               positive sign is redshifted, at time 1.
c - - Input:    lloncoe0,llatcoe0,lrllcoe0 - real*8(n+1,m+1) - the lon,lat,
c               and radial components of the LOS unit vector, at time 0.
c - - Input:    lloncoe1,llatcoe1,lrllcoe1 - real*8(n+1,m+1) - the lon,lat,
c               and radial components of the LOS unit vector, at time 1.
c - - Output:   data_ll_3d(n+1,m+1,18) - The 3D array consisting of the 18
c               2D arrays described above, in that order.
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
c - - Declaration of input variables:
c
      integer :: m,n
      real*8 :: bloncoe0(n+1,m+1),blatcoe0(n+1,m+1),brllcoe0(n+1,m+1)
      real*8 :: bloncoe1(n+1,m+1),blatcoe1(n+1,m+1),brllcoe1(n+1,m+1)
      real*8 :: vloncoe0(n+1,m+1),vlatcoe0(n+1,m+1),vlosllcoe0(n+1,m+1)
      real*8 :: vloncoe1(n+1,m+1),vlatcoe1(n+1,m+1),vlosllcoe1(n+1,m+1)
      real*8 :: lloncoe0(n+1,m+1),llatcoe0(n+1,m+1),lrllcoe0(n+1,m+1)
      real*8 :: lloncoe1(n+1,m+1),llatcoe1(n+1,m+1),lrllcoe1(n+1,m+1)
c
c - - Declaration of output variables:
c
      real*8 :: data_ll_3d(n+1,m+1,18)
c
c - - stack up the 2d arrays into the 3d array:
c
      data_ll_3d(:,:,1)=bloncoe0(:,:)
      data_ll_3d(:,:,2)=blatcoe0(:,:)
      data_ll_3d(:,:,3)=brllcoe0(:,:)
c
      data_ll_3d(:,:,4)=bloncoe1(:,:)
      data_ll_3d(:,:,5)=blatcoe1(:,:)
      data_ll_3d(:,:,6)=brllcoe1(:,:)
c
      data_ll_3d(:,:,7)=vloncoe0(:,:)
      data_ll_3d(:,:,8)=vlatcoe0(:,:)
      data_ll_3d(:,:,9)=vlosllcoe0(:,:)
c
      data_ll_3d(:,:,10)=vloncoe1(:,:)
      data_ll_3d(:,:,11)=vlatcoe1(:,:)
      data_ll_3d(:,:,12)=vlosllcoe1(:,:)
c
      data_ll_3d(:,:,13)=lloncoe0(:,:)
      data_ll_3d(:,:,14)=llatcoe0(:,:)
      data_ll_3d(:,:,15)=lrllcoe0(:,:)
c
      data_ll_3d(:,:,16)=lloncoe1(:,:)
      data_ll_3d(:,:,17)=llatcoe1(:,:)
      data_ll_3d(:,:,18)=lrllcoe1(:,:)
c
c - - we're done
c
      return
      end
