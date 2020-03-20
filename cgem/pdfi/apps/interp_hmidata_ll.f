      subroutine interp_hmidata_ll(m,n,data2d,m_new,n_new,
     1 data2d_new,degree)
c+
c - - Purpose: Interpolate a 2D Plate Carree HMI dataset from dimension
c              n+1,m+1 to a new 2D array of dimension n_new+1,m_new+1. COE grid
c              locations, lon,lat index order assumed for both input and output
c              arrays.
c
c - - Method:  Use bspline of order degree to interpolate from the original
c              array to the new array.  
c
c - - Usage:   call interp_hmidata_ll(m,n,data2d,m_new,n_new,data2d_new,
c              degree)
c
c - - Input:   m,n - integer values of the number of cell centers in the
c              latitude and longitude directions, respectively.
c
c - - Input:   data2d(n+1,m+1) - real*8 array of HMI data in longitude,
c              latitude order, on the COE grid
c
c - - Input:   m_new,n_new - integer values of cell centers in latitude and
c              longitude directions, respectively, for the new grid
c
c - - Output:  data2d_new(n_new+1,m_new+1) - real*8 array of interpolated
c              values at the new grid locations, also assumed COE grid,
c              in lon,lat index order
c
c - - Input:   degree - integer value of the degree of interpolation
c              (3 <= degree <=9), odd values only.
c
c - - Note:    degree=9 seems to work well.
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
c - - input argument declarations:
c
      integer :: m,n,m_new,n_new,degree
      real*8 :: data2d(n+1,m+1)
c
c - - output argument declatations:
c
      real*8 :: data2d_new(n_new+1,m_new+1)
c
c - - local variables:
c
      real*8 :: xnew(n_new+1), ynew(m_new+1)
      real*8 :: realnx,realny,realnxinterp,realnyinterp,one,reali,realj
      real*8 :: slopex,slopey
      integer :: nx,ny,nxinterp,nyinterp,i,j
c    
      if((degree/2)*2 .eq. degree) then
         write(6,*) 'interp_hmidata_ll:  degree must be odd; degree = ',
     1   degree
         stop
      endif
      if((degree .lt. 3) .or. (degree .gt. 9)) then
         write(6,*) 'interp_hmidata_ll:  degree < 3 or degree > 9; = ',
     1   degree
         stop
      endif
c
      nx=n+1
      ny=m+1
      nxinterp=n_new+1
      nyinterp=m_new+1
      realnx=nx
      realny=ny
      realnxinterp=nxinterp
      realnyinterp=nyinterp
      one=1
c
      slopex=(realnx-one)/(realnxinterp-one)
      slopey=(realny-one)/(realnyinterp-one)
      do i=1,nxinterp
        reali=i
        xnew(i) = one+slopex*(reali-one)
      enddo
      if(xnew(1) .lt. one) xnew(1) = one
      if(xnew(nxinterp) .gt. realnx) xnew(nxinterp)=realnx
c
      do j=1,nyinterp
         realj=j
         ynew(j)=one+slopey*(realj-one)
      enddo
      if(ynew(1) .lt. one) ynew(1) = one 
      if(ynew(nyinterp) .gt. realny) ynew(nyinterp) = realny
c
      call bspline_ss(data2d,nx,ny,xnew,nxinterp,ynew,nyinterp,
     1 data2d_new, degree)
c
c - - we're done
c
      return
      end
