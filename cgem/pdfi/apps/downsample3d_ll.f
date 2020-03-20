      subroutine downsample3d_ll(m,n,elon,elat,erll,delondr,delatdr,
     1           mc,nc,elonc,elatc,erllc,delondrc,delatdrc)
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
c - - Usage:    call downsample3d_ll(m,n,elon,elat,erll,delondr,delatdr,
c               mc,nc,elonc,elatc,erllc,delondrc,delatdrc)
c
c - - Input:    m,n - integers describing the number of cell
c               centers in the latitude, and longitudinal directions, 
c               respectively.
c
c - - Input:    elon(n,m+1) - real*8 array of c E_lon at full resolution
c               [G km/sec]
c
c - - Input:    elat(n+1,m) - real*8 array of c E_lat at full resolution
c               [G km/sec]
c
c - - Input:    erll(n+1,m+1) - real*8 array of c E_r (inductive) at full 
c               resolution [G km/sec]
c
c - - Input:    delondr(n,m+1) - real*8 array of c dE_lon / dr at full
c               resolution [G / sec]
c
c - - Input:    delatdr(n+1,m) - real*8 array of c dE_lat / dr at full
c               resolution [G / sec]
c
c - - Input:    mc, nc - integer values of the coarse grid resolution in
c               the latitude, longitude directions, respectively
c
c - - Output:   elonc(nc,mc+1) - real*8 array of c E_lon at coarse resolution
c               [G km / sec]
c
c - - Output:   elatc(nc+1,mc) - real*8 array of c E_lat at coarse resolution
c               [G km / sec]
c
c - - Output:   erllc(nc+1,mc+1) - real*8 array of c E_r (inductive) at coarse
c               resolution [G km / sec]
c
c - - Output:   delondrc(nc,mc+1) - real*8 array of c dE_lon / dr at coarse
c               resolution [G / sec]
c
c - - Output:   delatdrc(nc+1,mc) - real*8 array of c dE_lat / dr at coarse
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
      integer :: mc,nc
      real*8 :: elon(n,m+1),elat(n+1,m),erll(n+1,m+1)
      real*8 :: delondr(n,m+1),delatdr(n+1,m)
c
c - - output variable declarations:
c
      real*8 :: elatc(nc+1,mc),elonc(nc,mc+1),erllc(nc+1,mc+1)
      real*8 :: delondrc(nc,mc+1),delatdrc(nc+1,mc)
c
c - - local variable declarations:
c
      integer :: mrat,nrat
      real*8 :: et(m,n+1),ep(m+1,n),er(m+1,n+1)
      real*8 :: detdr(m,n+1),depdr(m+1,n)
      real*8 :: detdrc(mc,nc+1)
      real*8 :: depdrc(mc+1,nc)
      real*8 :: etc(mc,nc+1),epc(mc+1,nc),erc(mc+1,nc+1)
c
      mrat=m/mc
      nrat=n/nc
c
c - - test for divisibility of m/mc, n/nc:
c
      if((mrat*mc .ne. m) .or. (nrat*nc .ne. n)) then
        write(6,*) 'downsample3d_ll: m or n not divisible by mc, nc'
        write(6,*) 'downsample3d_ll: m,n,mc,nc = ',m,n,mc,nc
        stop
      endif
c
c - - rotate input arrays from lon-lat order to colat-lon order:
c
      call ehyeell2tp_ss(m,n,elon,elat,et,ep)
      call ehyeell2tp_ss(m,n,delondr,delatdr,detdr,depdr)
      call eryeell2tp_ss(m,n,erll,er)
c
c - - call downsample3d_ss:
c
      call downsample3d_ss(m,n,et,ep,er,detdr,depdr,mc,nc,etc,epc,erc,
     1     detdrc,depdrc)
c
c - - rotate coarse array results back to lon,lat order from theta,phi order:
c
      call ehyeetp2ll_ss(mc,nc,etc,epc,elonc,elatc)
      call ehyeetp2ll_ss(mc,nc,detdrc,depdrc,delondrc,delatdrc)
      call eryeetp2ll_ss(mc,nc,erc,erllc)
c
c - - we're done
c
      return
      end
