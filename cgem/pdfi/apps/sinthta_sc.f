      subroutine sinthta_sc(thmin,thmax,m,sinth,sinth_hlf,sinth_gh,
     1 sinth_hlf_gh)
c
c+
c - - Purpose: compute arrays of sin(theta) for the grid locations and at the
c              half-grid locations.  Compute arrays with and without ghost 
c              zones.
c
c              NOTE: this version uses a centered grid, not a staggered grid,
c              and is *NOT VALID* for the staggered grid variables.  It is
c              used *ONLY* within the context of the relax_psi_3d_ss subroutine,
c              which employs a centered grid formalism.
c
c - - Usage:   call(sinthta_sc,thmin,thmax,m,sinth,sinth_hlf,sinth_gh,
c              sinth_hlf_gh)
c
c - - Input:   thmin,thmax - real*8 values of the minimum and maximum of
c              co-latitude grid. [radians]
c
c              NOTE:  thmin,thmax are *NOT THE SAME* is the PDFI_SS limits a,b.
c              see source code details in relax_psi_3d_ss.
c
c - - Input:   m - integer number of grid interiors in the theta 
c              (colatitude) direction
c
c              NOTE:  The value is m here is two less than the value of m
c              in the staggered grid formalism.
c     
c - - Output:  sinth(m+1), a real*8 array of grid values of sin(theta)
c
c - - Output:  sinth_gh(m+3), a real*8 array of grid values of sin(theta)
c              including one ghost-zone at each end
c
c - - Output:  sinth_hlf(m), a real*8 array of half-grid values of sin(theta)
c
c - - Output:  sinth_hlf_gh(m+2), a real*8 array of half-grid values of
c              sin(theta) including one ghost-zone at each end
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
      integer :: m
      real*8 thmin,thmax
      real*8 sinth(m+1),sinth_gh(m+3),sinth_hlf(m),sinth_hlf_gh(m+2)
c
      integer :: i
      real*8 :: theta(m+1),theta_gh(m+3),theta_hlf(m),
     1 theta_hlf_gh(m+2)
c
      real*8 :: deltheta,deltheta_hlf
      deltheta=(thmax-thmin)/m
      deltheta_hlf=0.5d0*deltheta
      do i=1,m
         theta_hlf(i)=thmin+0.5d0*deltheta_hlf+(i-1)*deltheta
         theta_hlf_gh(i+1)=theta_hlf(i)
         theta(i)=thmin+(i-1)*deltheta
         theta_gh(i+1)=theta(i)
      enddo
      theta_hlf_gh(1)=theta_hlf(1)-deltheta
      theta(m+1)=thmax
      theta_gh(m+2)=thmax
      theta_gh(1)=thmin-deltheta
      theta_gh(m+3)=thmax+deltheta
      theta_hlf_gh(m+2)=thmax+0.5d0*deltheta
      sinth(1:m+1)=sin(theta(1:m+1))
      sinth_gh(1:m+3)=sin(theta_gh(1:m+3))
      sinth_hlf(1:m)=sin(theta_hlf(1:m))
      sinth_hlf_gh(1:m+2)=sin(theta_hlf_gh(1:m+2))
c
      return
      end
