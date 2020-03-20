      subroutine downsample3d_ss(m,n,et,ep,er,detdr,depdr,
     1           mc,nc,etc,epc,erc,detdrc,depdrc)
c
c+
c - - Purpose:  To downsample the 3 components of the electric field, plus
c               the radial derivative of the horizontal components, from
c               array dimensions m,n to array dimensions mc,nc.  The ratios
c               m/mc and n/nc must be whole integers.
c
c - - Method:   Take line integral around the edges of the course voxels for
c               each of the faces of the layer of voxels, and make sure
c               Faraday's law is obeyed for each face.  For horizontal E-fields,
c               and radial derivative of horizontal E-fields, each rail value
c               is a simple arithmetic average of the fine grid line segments
c               within each coarse rail length, since the fine-grid line
c               segments are of equal length and sum to coarse grid line
c               segment length.  For radial electric field, we simply sample
c               at the coarse grid vertices (COE grid).
c
c - - Usage:    call downsample3d_ss(m,n,et,ep,er,detdr,depdr,
c               mc,nc,etc,epc,erc,detdrc,depdrc
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
c - - Input:    er(m+1,n+1) - real*8 array of c E_r (inductive) at full 
c               resolution [G km/sec]
c
c - - Input:    detdr(m,n+1) - real*8 array of c dE_theta / dr at full
c               resolution [G / sec]
c
c - - Input:    depdr(m+1,n) - real*8 array of c dE_phi / dr at full
c               resolution [G / sec]
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
c - - Output:   erc(mc+1,nc+1) - real*8 array of c E_r (inductive) at coarse
c               resolution [G km / sec]
c
c - - Output:   detdrc(mc,nc+1) - real*8 array of c dE_theta / dr at coarse
c               resolution [G / sec]
c
c - - Output:   depdrc(mc+1,nc) - real*8 array of c dE_phi / dr at coarse
c               resolution [G / sec]
c
c-
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
c - - input variable declarations:
c
      integer :: m,n
      real*8 :: et(m,n+1),ep(m+1,n),er(m+1,n+1)
      real*8 :: detdr(m,n+1),depdr(m+1,n)
      integer :: mc,nc
c
c - - output variable declarations:
c
      real*8 :: detdrc(mc,nc+1)
      real*8 :: depdrc(mc+1,nc)
      real*8 :: etc(mc,nc+1),epc(mc+1,nc),erc(mc+1,nc+1)
c
c - - local variable declarations:
c
      integer :: mrat,nrat,icph,jc,jcph,ic,i,j,iphlow,iphhigh
      integer :: iph,jph,jphlow,jphhigh
      real*8 :: ethbar,dethbardr,ephbar,dephbardr
c
      mrat=m/mc
      nrat=n/nc
c
c - - test for divisibility of m/mc, n/nc:
c
      if((mrat*mc .ne. m) .or. (nrat*nc .ne. n)) then
        write(6,*) 'downsample3d_ss: m or n not divisible by mc, nc'
        write(6,*) 'downsample3d_ss: m,n,mc,nc = ',m,n,mc,nc
        stop
      endif
c
c - - Phi edges (PE grid) (get coarse values for E_theta, d E_theta/dr):
c
      do jc=1,nc+1
         do icph=1,mc
            j=nrat*(jc-1)+1
            iphlow=mrat*(icph-1)+1
            iphhigh=mrat*(icph-1)+mrat
            ethbar=0.d0
            dethbardr=0.d0
c
c - - find average of the high res values of et, detdr:
c
            do iph=iphlow,iphhigh
               ethbar=ethbar+et(iph,j)
               dethbardr=dethbardr+detdr(iph,j)
            end do
            ethbar=ethbar/dble(mrat)
            dethbardr=dethbardr/dble(mrat)
c
c - - set coarse values of E_theta, d E_theta/dr equal to average values:
c
            etc(icph,jc)=ethbar
            detdrc(icph,jc)=dethbardr
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
            dephbardr=0.d0
c
c - - Get average value of eph, dephdr from high res grid:
c
            do jph=jphlow,jphhigh
               ephbar=ephbar+ep(i,jph)
               dephbardr=dephbardr+depdr(i,jph)
            end do
            ephbar=ephbar/dble(nrat)
            dephbardr=dephbardr/dble(nrat)
c
c - - Set coarse grid values of eph, dephdr to average values:
c
            epc(ic,jcph)=ephbar
            depdrc(ic,jcph)=dephbardr
         end do
      end do
c
c - - Set coarse values of er (COE grid) (sampled at coarse grid coordinates 
c     from fine grid)
c
      do jc=1,nc+1
         do ic=1,mc+1
            i=mrat*(ic-1)+1
            j=nrat*(jc-1)+1
            erc(ic,jc)=er(i,j)
         end do
      end do
c
c - - we're done
c
      return
      end
