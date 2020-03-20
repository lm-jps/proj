      subroutine pad_abcd_as_ss(m,n,mpadb,mpadt,npadl,npadr,a,b,c,d,
     1           mp,np,ap,bp,cp,dp)
c
c+
c - - Purpose:  To correct the values of a,b,c,d to reflect asymmetric buffers 
c               of zero-padding around a 2-d array containing non-zero values. 
c               Also compute corrected values of m,n (ie mp,np) to reflect the 
c               padding.
c
c - - Usage:    call pad_abcd_as_ss(m,n,mpadb,mpadt,npadl,npadr,a,b,c,d,
c               mp,np,ap,bp,cp,dp)
c
c - - Input:    m,n - integers describing the number of cell interiors
c               in latitude/colatitude, and in longitude, respectively, 
c               of the unpadded array.  
c
c - - Input:    mpadb,mpadt, - integers descibing the number of 0-padded
c               cells in latitude at the bottom and top, respectively, 
c               that will surround the unpadded array.
c
c - - Input:    npadl,npadr - integers descibing the number of 0-padded
c               cells in longitude at left and right, respectively, that will 
c               surround the unpadded array.
c
c - - Input:    a,b,c,d - real*8 variables describing the minimum and
c               maximum co-latitude, and minimum and maximum longitude, 
c               respectively, spanned by the unpadded array [radians]
c
c - - Output:   mp,np - integers describing the number of cell interiors in 
c               latitude/colatitude, and in longitude, respectively, of the
c               padded array.  mp=m+mpadb+mpadt, and np=n+npadl+npadr
c
c - - Output:   ap,bp,cp,dp - real*8 variables describing minimum and maximum
c               colatitude and longitude of the padded array. [radians]
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
c - - Input argument declarations:
c
      integer :: m,n,mpadb,mpadt,npadl,npadr
      real*8 :: a,b,c,d
c
c - - Output argument declarations:
c
      integer :: mp,np
      real*8 :: ap,bp,cp,dp
c
c - - local variable declarations:
c
      real*8 :: dum,pi,twopi,dtheta,dphi
c
c - - declaration of pimach function from FISHPACK/FFTPACK:
c
      real*8 :: pimach
c
      pi=pimach(dum)
      twopi=2.d0*pimach(dum)
c
      dtheta=(b-a)/m
      dphi=(d-c)/n
c
      ap=a-mpadt*dtheta
      bp=b+mpadb*dtheta
      cp=c-npadl*dphi
      dp=d+npadr*dphi
c
c - - error exit if ap < 0 or bp > pi or dp-cp > 2*pi
c
      if((ap .le. 0.d0) .or. (bp .gt. pi) .or. (dp-cp .gt. twopi)) then
        write(6,*) 'pad_abcd_as_ss: ',
     1  'ap, bp, or dp-cp out of range: ap,bp,dp-cp = ',ap,bp,dp-cp
        stop
      endif
      mp=m+mpadb+mpadt
      np=n+npadl+npadr
c
      return
      end
