      subroutine fluxbal_ll(m,n,a,b,c,d,rsun,brtll,brtllbal)
c
c+
c - - Purpose:  Given an array of cell-center radial magnetic field values,
c               this subroutine removes any flux imbalance by enhancing the
c               field strength values in the defficient flux pixels.  This
c               is done in a way that is proportional to the initial field
c               values, rather than subtracting the extra flux uniformly from
c               the array.  The objective is to avoid changing the locations
c               of PILs as much as possible.  lon,lat array order assumed.
c - - Usage:    call fluxbal_ll(m,n,a,b,c,d,rsun,brtll,brtllbal)
c - - Input:    m,n - integer values of the number of cell centers in the
c               colatitude and longitude direction, respectively
c - - Input:    a,b,c,d - real*8 values of the colatitude limits (a,b)
c               and longitude limits (c,d) of the spherical domain [radians].  
c               Spherical wedge or global limits can be used, but if 
c               b .eq. pi, you must use the pimach function from FISHPACK 
c               to set b=pi.
c - - Input:    rsun - real*8 value of assumed radius for Sun [km].  Usual
c               value is 6.96e5
c - - Input:    brtll(n,m) - real*8 array of radial magnetic field values on CE
c               grid [G].  brtll(n,m) can have a flux imbalance.
c - - Output:   brtllbal(n,m) - real*8 array of radial magnetic field values on
c               CE grid, with flux imbalance removed. [G] lon,lat array order.
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
      real*8 :: brtll(n,m)
c
c - - Output variable declarations:
c
      real*8 :: brtllbal(n,m)
c
c - - Local variable declarations:
c
      real*8 :: brt(m,n),brtbal(m,n)
c
c - - transpose brtll array to brt array:
c
      call bryeell2tp_ss(m,n,brtll,brt)
c
c - - flux balance correction on brt to brtbal:
c
      call fluxbal_ss(m,n,a,b,c,d,rsun,brt,brtbal)
c
c - - transpose brtbal to lon,lat order (to brtllbal):
c
      call bryeetp2ll_ss(m,n,brtbal,brtllbal)
c
c - - we're done
c
      return
      end
