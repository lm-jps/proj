      subroutine add_padding_ss(m,n,npadlat,npadlon,arr,arrpadded)
c
c+
c - - Purpose:  To add a buffer of zero-padding around a 2-d array containing
c               non-zero values.  Padding is assumed symmetric between 
c               N and S edges of domain, and between left and right edges.
c
c - - Usage:    call add_padding_ss(m,n,npadlat,npadlon,arrs,arrpadded)
c
c - - Input:    m,n - integers describing the number of cell interiors
c               in latitude/colatitude, and in longitude, respectively, 
c               of the unpadded array.  
c
c - - Input:    npadlat,npadlon - integers descibing the number of 0-padded
c               cells in latitude and longitude, respectively, that will 
c               surround the unpadded array.
c
c - - Input:    arr(n+1,m+1) - real*8 array of the unpadded array values,
c               stored in lon,lat index order on the COE grid.
c
c - - Output:   arpadded(n+2*npadlon+1,m+2*npadlat+1) - real*8 array of 
c               of the padded array values.  Note array is stored in lon,lat
c               index order, same as input array, also on COE grid, but with
c               increase in assumed array size.
c
c - - Note:     This subroutine adds padding symmetrically on left and right
c               sides, and also on north and south sides.
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
      integer :: m,n,npadlat,npadlon
      real*8 :: arr(n+1,m+1)
c
c - - Output argument declarations:
c
      real*8 :: arrpadded(n+2*npadlon+1,m+2*npadlat+1)
c
c - - local variable declarations:
c
      integer :: mp,np
c
      mp=m+2*npadlat
      np=n+2*npadlon
c - - Initialize padded array to 0:
      arrpadded(1:np+1,1:mp+1)=0.d0
c - - fill in non-zero values:
      arrpadded(1+npadlon:n+1+npadlon,1+npadlat:m+1+npadlat) =
     1     arr(1:n+1,1:m+1)
c
      return
      end
