      subroutine fluxbal_ss(m,n,a,b,c,d,rsun,brt,brtbal)
c
c+
c - - Purpose:  Given an array of cell-center radial magnetic field values,
c               this subroutine removes any flux imbalance by enhancing the
c               field strength values in the defficient flux pixels.  This
c               is done in a way that is proportional to the initial field
c               values, rather than by subtracting the extra flux uniformly from
c               the array.  The objective is to avoid changing the locations
c               of PILs as much as possible.
c
c - - Usage:    call fluxbal_ss(m,n,a,b,c,d,rsun,brt,brtbal)
c
c - - Input:    m,n - integer values of the number of cell centers in the
c               colatitude and longitude direction, respectively
c
c - - Input:    a,b,c,d - real*8 values of the colatitude limits (a,b)
c               and longitude limits (c,d) of the spherical domain [radians].  
c               Spherical wedge or global limits can be used, but if 
c               b .eq. pi, you must use the pimach function from FISHPACK 
c               to set b=pi.
c
c - - Input:    rsun - real*8 value of assumed radius for Sun [km].  Usual
c               value is 6.96e5
c
c - - Input:    brt(m,n) - real*8 array of radial magnetic field values on CE
c               grid [G].  brt(m,n) can have a flux imbalance.
c
c - - Output:   brtbal(m,n) - real*8 array of radial magnetic field values on
c               CE grid, with flux imbalance removed. [G]
c-
c   PDFI_SS Electric Field Inversion Software
c   http://cgem.ssl.berkeley.edu/cgi-bin/cgem/PDFI_SS/index
c   Copyright (C) 2015-2018 University of California
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
c - - Input variable declarations:
c
      integer :: m,n
      real*8 :: a,b,c,d,rsun
      real*8 :: brt(m,n)
c
c - - Output variable declarations:
c
      real*8 :: brtbal(m,n)
c
c - - Local variable declarations:
c
      integer :: i,j,mask(m,n)
      real*8 :: posflux,negflux,netflux,corfac
      real*8 :: rsun2,dtheta,dphi,oneodp,oneodt
      real*8 :: sinth(m+1),sinth_hlf(m)
c
c - - define dtheta, dphi, rsun2:
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
      oneodp=1.d0/dphi
      oneodt=1.d0/dtheta
      rsun2=rsun**2
c
c - - compute sinth, sinth_hlf arrays from call to sinthta_ss:
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - set mask to 1 where brt ge 0, and -1 where brt lt 0
c
      where(brt .ge. 0.d0)
           mask=1
      elsewhere
           mask=-1
      endwhere
c
c - - calculate netflux, positive flux, and negative flux:
c
      netflux=0.d0
      posflux=0.d0
      negflux=0.d0
      do i=1,m
         do j=1,n
            netflux=netflux+brt(i,j)*sinth_hlf(i)*dtheta*dphi*rsun2
            if(mask(i,j) .eq. 1) then
              posflux=posflux+brt(i,j)*sinth_hlf(i)*dtheta*dphi*rsun2
            else
              negflux=negflux+brt(i,j)*sinth_hlf(i)*dtheta*dphi*rsun2
            endif
         enddo
      enddo
c
c - - enhance the under-represented flux to flux-balance the two polarities:
c
      if(netflux .gt. 0.d0) then
         corfac=1.d0+abs(netflux)/(abs(negflux))
         where(mask .eq. -1)
           brtbal=brt*corfac
         elsewhere
           brtbal=brt
         endwhere
      else
         corfac=1.d0+abs(netflux)/(abs(posflux))
         where(mask .eq. 1)
           brtbal=brt*corfac
         elsewhere
           brtbal=brt
         endwhere
      endif
c
c - - we're done
c
      return
      end
