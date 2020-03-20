      subroutine ehyeell2tp_ss(m,n,elon,elat,etpe,epte)
c
c+
c    Purpose: To transpose E_h data arrays from lon,lat to theta,phi order
c             and flip sign to get E_theta.  
c
c     Usage:  call ehyeell2tp_ss(m,n,elon,elat,etpe,epte)
c
c     Input:  m,n - integer number of cell centers in the theta (lat), and 
c             phi (lon) directions, respectively.
c
c     Input:  elon(n,m+1),elat(n+1,m) - real*8 arrays of longitudinal and
c             latitudinal components of electric field, at TE and PE grid
c             locations, respectively, and stored in lon,lat
c             index order. [G km/sec or V/cm]
c
c    Output:  etpe (m,n+1),epte(m+1,n) - real*8 arrays of the co-latitudinal and
c             azimuthal components of the electric field evaluated at
c             PE and TE locations (phi and theta edges, resp.)
c             [G km/sec or V/cm]
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
c - - input arguments
c
      integer :: m,n
      real*8 :: elat(n+1,m),elon(n,m+1)
c
c - - output arguments:
c
      real*8 :: etpe(m,n+1),epte(m+1,n)
c
c - - local variables:
c
      integer :: i,j
c
      do i=1,m+1
         do j=1,n
            epte(m+2-i,j)=elon(j,i)
         enddo
      enddo
c
c - - lat and colat unit vectors have opp. sign so must change sign:
c
      do i=1,m
         do j=1,n+1
            etpe(m+1-i,j)=-elat(j,i)
         enddo
      enddo
c
      return
      end
