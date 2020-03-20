      subroutine downsample_ss(m,n,et,ep,mc,nc,etc,epc)
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
c               is equal.
c
c - - Usage:    call downsample_ss(m,n,et,ep,mc,nc,etc,epc)
c
c - - Input:    m,n - integers describing the number of cell
c               centers in the colatitude, and longitudinal directions, 
c               respectively.
c
c - - Input:    et(m,n+1) - real*8 array of c E_theta at full resolution
c               [G km/sec]
c
c - - Input:    ep(m+1,n) - real*8 array of c E_phi at full resolution
c               [G km/sec]
c
c - - Input:    mc, nc - integer values of the coarse grid resolution in
c               the colatitude, longitude directions, respectively
c
c - - Output:   etc(mc,nc+1) - real*8 array of c E_theta at coarse resolution
c               [G km / sec]
c
c - - Output:   epc(mc+1,nc) - real*8 array of c E_phi at coarse resolution
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
      real*8 :: et(m,n+1),ep(m+1,n)
      integer :: mc,nc
c
c - - output variable declarations:
c
      real*8 :: etc(mc,nc+1),epc(mc+1,nc)
c
c - - local variable declarations:
c
      integer :: mrat,nrat,icph,jc,jcph,ic,i,j,iphlow,iphhigh
      integer :: iph,jph,jphlow,jphhigh
      real*8 :: ethbar,ephbar
c
      mrat=m/mc
      nrat=n/nc
c
c - - test for integer divisibility of m/mc, n/nc:
c
      if((mrat*mc .ne. m) .or. (nrat*nc .ne. n)) then
        write(6,*) 'downsample_ss: m or n not divisible by mc, nc'
        write(6,*) 'downsample_ss: m, n, mc, nc = ',m,n,mc,nc
        stop
      endif
c
c - - Phi edges (PE grid) (get coarse values for E_theta):
c
      do jc=1,nc+1
         do icph=1,mc
            j=nrat*(jc-1)+1
            iphlow=mrat*(icph-1)+1
            iphhigh=mrat*(icph-1)+mrat
            ethbar=0.d0
c
c - - find average of the high res values of et:
c
            do iph=iphlow,iphhigh
               ethbar=ethbar+et(iph,j)
            end do
            ethbar=ethbar/dble(mrat)
c
c - - set coarse values of E_theta equal to average value:
c
            etc(icph,jc)=ethbar
         end do
      end do
c
c - - Theta edges (TE grid) (Get coarse values of E_phi and d E_phi/dr):
c
      do ic=1,mc+1
         do jcph=1,nc
            i=mrat*(ic-1)+1
            jphlow=nrat*(jcph-1)+1
            jphhigh=nrat*(jcph-1)+nrat
            ephbar=0.d0
c
c - - Get average value of eph, dephdr from high res grid:
c
            do jph=jphlow,jphhigh
               ephbar=ephbar+ep(i,jph)
            end do
            ephbar=ephbar/dble(nrat)
c
c - - Set coarse grid values of eph to average value:
c
            epc(ic,jcph)=ephbar
         end do
      end do
c
c - - we're done
c
      return
      end
