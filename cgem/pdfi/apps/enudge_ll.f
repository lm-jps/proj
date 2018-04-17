      subroutine enudge_ll(m,n,a,b,c,d,rsun,brtll,elon,elat)
c+
c - -  Purpose:  Compute horizontal electric fields elon,elat from radial 
c      magnetic field time derivatives brt, using the simplest PTD electric 
c      field solution. Assumes for boundary conditions that the tangential 
c      electric field at the outer boundaries is zero.  
c
c - -  Usage:  call enudge_ll(m,n,a,b,c,d,rsun,brtll,elon,elat)
c
c - -  Input:  m,n - integers denoting the number of cell-centers in the
c              colatitude and longitude directions, respectively
c - -  Input:  a,b - the real*8 values of colatitude (theta) 
c              at the northern and southern edges of the problem boundary
c - -  Input:  c,d - the real*8 values of longitude edges
c - -  Input:  rsun - real*8 - units for the radius of the Sun
c - -  Input:  brtll(n,m) - real*8 array of cell-center radial magnetic field
c              time derivative values, arranged in lon,lat order
c - -  Output: elon(n,m+1) - real*8 array of longitudinal electric fields, 
c              multiplied by the speed of light, computed on theta edges 
c              [G km/s].
c - -  Output: elat(n+1,m) - real*8 array of latitudinal electric fields, 
c              multiplied by speed of light, computed on phi edges [G km/s].
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
c
      implicit none
c
c - - variable declarations:
c
c - - calling arguments:
c
      integer :: m,n
      real*8 :: brtll(n,m),elat(n+1,m),elon(n,m+1)
      real*8 :: rsun,a,b,c,d
c
c - - local subroutine variables:
c
      real*8 :: brt(m,n)
      real*8 :: et(m,n+1),ep(m+1,n)
c
c - - Transpose input brtll data from lon,lat to theta,phi order:
c
      call bryeell2tp_ss(m,n,brtll,brt)
c
c - - Now call spherical-polar coordinate version, enudge_ss:
c
      call enudge_ss(m,n,a,b,c,d,rsun,brt,et,ep)
c
c - - Transpose et,ep to elon,elat:
c
      call ehyeetp2ll_ss(m,n,et,ep,elon,elat)
c
c - - we're done
c
      return
      end
