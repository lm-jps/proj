      subroutine delh2_co_ss(m,n,a,b,c,d,rsun,psi,delh2)
c
c+ 
c - -  Purpose: Compute horizontal laplacian of psi evaluated on the co grid.
c
c - -   Usage:  call delh2_co_ss (m,n,a,b,c,d,rsun,psi,delh2)
c
c - -   Input:  m,n - integer number of cell centers in theta and phi 
c               directions, resp.
c
c - -   Input:  a,b,c,d - real*8 values of the colatitude boundaries of the
c               domain (a < b) and longitude boundaries (c < d) [radians]
c
c - -   Input:  rsun - real*8 value of radius of Sun [km]. Normally 6.96d5
c
c - -   Input:  psi(m+1,n+1) - real*8 array of scalar potential (or its
c               time derivative) on COE grid [G km^2 or G km^2/s]
c
c - -  Output:  delh2(m-1,n-1) - real*8 array array of horizontal laplacian of 
c               psi at interior corners (CO grid). [G or G/(sec)]
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
      real*8 :: a,b,c,d,rsun
      real*8 :: psi(m+1,n+1)
c
c - - output variable declarations:
c
      real*8 :: delh2(m-1,n-1)
c
c - - local variable declarations:
c
      real*8 :: dtheta,dphi
      real*8 :: sinth_hlf(m),sinth(m+1)
      real*8 :: gradt(m,n+1),gradp(m+1,n)
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
      call gradh_co_ss(m,n,psi,rsun,sinth_hlf,dtheta,dphi,gradt,gradp)
      call divh_co_ss(m,n,gradt,gradp,rsun,sinth,sinth_hlf,dtheta,dphi,
     1     delh2)
c
      return
      end
