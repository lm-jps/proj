      subroutine sinthta_ss(thmin,thmax,m,sinth,sinth_hlf)
c+
c - - Purpose: compute arrays of sin(theta) for the theta edge locations and 
c              at the half-cell locations.  
c
c - - Usage:   call sinthta_ss(thmin,thmax,m,sinth,sinth_hlf)
c
c - - Input:   thmin, thmax - real*8 minimum and maximum values of
c              co-latitude theta at theta edges of domain [radians]
c
c - - Input:   m -  integer number of cell interiors in the theta (colatitude) 
c              direction
c
c - - Output:  sinth(m+1) - real*8 array of cell edge values sin(colatitude)
c
c - - Output:  sinth_hlf(m) -  real*8 array of cell center values of 
c              sin(colatitude)
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
c
      implicit none
c
c - - input variables:
c
      integer :: m
      real*8 :: thmin,thmax
c
c - - output variables:
c
      real*8 :: sinth(m+1),sinth_hlf(m)
c
c - - local variables:
c
      integer :: i
      real*8 :: theta(m+1),theta_hlf(m)
      real*8 :: dth,dth_hlf
      real*8 :: pi,dum
c
c - - function declarations:
c
      real*8 :: pimach
c
c - - In production runs, use pimach() from fishpack to set value of pi
c
      pi=pimach(dum)
c
c
      if((thmin .lt. 0.d0) .or. (thmax .gt. pi)) then
        write(6,*) 'sinthta_ss: bad range, thmin, thmax = ',thmin,thmax
        stop
      endif
c
      dth=(thmax-thmin)/m
      dth_hlf=0.5d0*dth
      do i=1,m
         theta(i)=thmin + (i-1)*dth
         theta_hlf(i)=thmin+dth_hlf+(i-1)*dth
      enddo
      theta(m+1)=thmax
c
      sinth(:)=sin(theta(:))
      sinth_hlf(:)=sin(theta_hlf(:))
c
      return
      end
