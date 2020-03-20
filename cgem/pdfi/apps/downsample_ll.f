      subroutine downsample_ll(m,n,elon,elat,mc,nc,elonc,elatc)
c
c+
c - - Purpose:  To downsample the 2 horizontal components of the electric 
c               field, plus array dimensions m,n to array dimensions mc,nc.  
c               The ratios m/mc and n/nc must be whole integers.
c
c - - Method:   Take line integral around the edges of the course voxels for
c               each of the faces of the layer of voxels, and make sure
c               Faraday's law is obeyed for each face.  For each coarse rail,
c               electric field on rail is straight arithmetic average of fine
c               grid cell values, since the length of each sub-rail segment
c               is equal.  This version assumes longitude,latitude array
c               orientation.
c
c - - Usage:    call downsample_ll(m,n,elon,elat,mc,nc,elonc,elatc)
c
c - - Input:    m,n - integers describing the number of cell
c               centers in the colatitude, and longitudinal directions, 
c               respectively.
c
c - - Input:    elon(n,m+1) - real*8 array of c E_lon at full resolution
c               [G km/sec]
c
c - - Input:    elat(n+1,m) - real*8 array of c E_lat at full resolution
c               [G km/sec]
c
c - - Input:    mc, nc - integer values of the coarse grid resolution in
c               the latitude, longitude directions, respectively
c
c - - Output:   elonc(nc,mc+1) - real*8 array of c E_lon at coarse resolution
c               [G km / sec]
c
c - - Output:   elatc(nc+1,mc) - real*8 array of c E_theta at coarse resolution
c               [G km / sec]
c
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
      integer :: mc,nc
      real*8 :: elon(n,m+1),elat(n+1,m)
c
c - - output variable declarations:
c
      real*8 :: elatc(nc+1,mc),elonc(nc,mc+1)
c
c - - local variable declarations:
c
      integer :: mrat,nrat
      real*8 :: et(m,n+1),ep(m+1,n)
      real*8 :: etc(mc,nc+1),epc(mc+1,nc)
c
      mrat=m/mc
      nrat=n/nc
c
c - - test for integer divisibility of m/mc, n/nc:
c
      if((mrat*mc .ne. m) .or. (nrat*nc .ne. n)) then
        write(6,*) 'downsample_ll: m or n not divisible by mc, nc'
        write(6,*) 'downsample_ll: m, n, mc, nc = ',m,n,mc,nc
        stop
      endif
c
c - - rotate elon,elat to et,ep:
c
      call ehyeell2tp_ss(m,n,elon,elat,et,ep)
c
c - - call downsample_ss:
c
      call downsample_ss(m,n,et,ep,mc,nc,etc,epc)
c
c - - rotate etc,epc to elonc,elatc:
c
      call ehyeetp2ll_ss(mc,nc,etc,epc,elonc,elatc)
c
c - - we're done
c
      return
      end
