      subroutine interp_hmidata_3d_ll(m,n,data3d,m_int,n_int,data3d_int)
c
c+
c - - Purpose:  To interpolate all 18 2D input arrays from original resolution
c               (n+1,m+1) to new resolution (n_int+1,m_int+1).  On input
c               the 2d arrays are assumed to be stacked into a 3D array of
c               dimension (n+1,m+1,18), and on output the interpolated 2D
c               arrays are stacked into a 3d array of dimension
c               (n_new+1,m_new+1,18.
c
c - - Usage:    call interp_hmidata_3d_ll(m,n,data3d,m_int,n_int,data3d_int)
c 
c - - Input:    m,n - integers defining the number of cell-centers in the
c               latitude (or co-latitude) direction, respectively.
c
c - - Input:    data3d(n+1,m+1,18) - The 18 input arrays stacked into a 3d
c               real*8 array at COE grid locations, in lon,lat index order.  
c               The order of the 18 input arrays is, referring to
c               the arguments of pdfi_wrapper4jsoc_ss, are:
c               id=1: data3d(:,:,id) = bloncoe0(:,:)
c               id=2: data3d(:,:,id) = blatcoe0(:,:)
c               id=3: data3d(:,:,id) = brllcoe0(:,:)
c               id=4: data3d(:,:,id) = bloncoe1(:,:)
c               id=5: data3d(:,:,id) = blatcoe1(:,:)
c               id=6: data3d(:,:,id) = brllcoe1(:,:)
c               id=7: data3d(:,:,id) = vloncoe0(:,:)
c               id=8: data3d(:,:,id) = vlatcoe0(:,:)
c               id=9: data3d(:,:,id) = vlosllcoe0(:,:)
c               id=10: data3d(:,:,id) = vloncoe1(:,:)
c               id=11: data3d(:,:,id) = vlatcoe1(:,:)
c               id=12: data3d(:,:,id) = vlosllcoe1(:,:)
c               id=13: data3d(:,:,id) = lloncoe0(:,:)
c               id=14: data3d(:,:,id) = llatcoe0(:,:)
c               id=15: data3d(:,:,id) = lrllcoe0(:,:)
c               id=16: data3d(:,:,id) = lloncoe1(:,:)
c               id=17: data3d(:,:,id) = llatcoe1(:,:)
c               id=18: data3d(:,:,id) = lrllcoe1(:,:)
c
c - - Input:    m_int,n_int - integer values of of the number of cell-centers
c               in the latitude (co-latitude), and longitude directions, 
c               respectively, for the interpolated grid.
c
c - - Output:   data3d(n_int+1,m_int+1,18) -  real*8 array of the interpolated 
c               values for all 18 input arrays using the new values 
c               m_int, n_int. COE grid locations in each 2D slice, in lon,lat
c               index order.
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
c - - Input variable declarations:
c
      integer :: m,n,m_int,n_int
      real*8 :: data3d(n+1,m+1,18)
c
c - - Output variable declaration
c
      real*8 :: data3d_int(n_int+1,m_int+1,18)
c
c - - Local variable declarations:
c
      integer :: id, degree
      real*8 data2d(n+1,m+1),data2d_int(n_int+1,m_int+1)
c
c - - For now, assume default value of degree = 9:
c
      degree=9
c
      do id=1,18
        data2d(1:n+1,1:m+1)=data3d(1:n+1,1:m+1,id)
        call interp_hmidata_ll(m,n,data2d,m_int,n_int,data2d_int,degree)
        data3d_int(1:n_int+1,1:m_int+1,id) = data2d_int(1:n_int+1,
     1            1:m_int+1)
      enddo
c
      return
      end
