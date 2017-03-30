      subroutine add_padding_ss(ms,ns,npadlat,npadlon,as,bs,cs,ds,arrs,
     1           m,n,a,b,c,d,arrpadded)
c
c+
c - - Purpose:  To add a buffer of zero-padding around a 2-d array containing
c               non-zero values.
c - - Usage:    call add_padding_ss(ms,ns,npadlat,npadlon,as,bs,cs,ds,arrs,
c               m,n,a,b,c,d,arrpadded)
c - - Input:    ms,ns - integers describing the number of cell interiors
c               in latitude/colatitude, and in longitude, respectively, 
c               of the unpadded array.  
c - - Input:    npadlat,npadlon - integers descibing the number of 0-padded
c               cells in latitude and longitude, respectively, that will 
c               surround the unpadded array.
c - - Input:    as,bs,cs,ds - real*8 variables describing the minimum and
c               maximum co-latitude, and minimum and maximum longitude, 
c               respectively, spanned by the unpadded array
c - - Input:    arrs(ns+1,ms+1) - real*8 array of the unpadded array values,
c               stored in lon,lat index order.
c - - Output:   m,n - integers describing the number of cell interiors in 
c               latitude/colatitude, and in longitude, respectively, of the
c               padded array.  m=ms+2*npadlat, and n=ns+2*npadlon
c - - Output:   a,b,c,d - real*8 variables describing the minimum and maximum
c               colatitude and longitude of the padded array.
c - - Output:   arpadded(ns+2*npadlon+1,ms+2*npadlat+1) - real*8 array of 
c               of the padded array values.  Note array is stored in lon,lat
c               index order, same as input array.
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
      implicit none
c
c - - Input argument declarations:
c
      integer :: ms,ns,npadlat,npadlon
      real*8 :: as,bs,cs,ds
      real*8 :: arrs(ns,ms)
c
c - - Output argument declarations:
c
      integer :: m,n
      real*8 :: a,b,c,d
      real*8 :: arrpadded(ns+2*npadlon+1,ms+2*npadlat+1)
c
c - - local variable declarations:
c
      real*8 :: dtheta,dphi
c
      dtheta=(bs-as)/ms
      dphi=(ds-cs)/ns
c
      a=as-npadlat*dtheta
      b=bs+npadlat*dtheta
      c=cs-npadlon*dphi
      d=ds+npadlon*dphi
c
      m=ms+2*npadlat
      n=ns+2*npadlon
c - - Initialize padded array to 0:
      arrpadded(:,:)=0.d0
c - - fill in non-zero values:
      arrpadded(1+npadlon:ns+1+npadlon,1+npadlat:ms+1+npadlat) =
     1     arrs(1:ns+1,1:ms+1)
c
      return
      end
