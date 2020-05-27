      subroutine e_laplace_ll(m,n,a,b,c,d,rsun,en,es,el,er,elon,elat)
c+
c - -  Purpose: Compute horizontal electric fields elon,elat from solutions to
c               Laplace equation, given boundary conditions en,es,el,er 
c               on the 4 sides.  Output is converted to longitude, latitude 
c               format, including sign change from e_theta to elat.
c               Assumes for boundary conditions that the tangential electric 
c               field at the N,S boundaries is given by en and es, and at 
c               the left and right longitude boundaries is given by el, er.
c
c - -  Usage:   call e_laplace_ll(m,n,a,b,c,d,rsun,en,es,el,er,elon,elat)
c
c - -  Input:   m,n - integers denoting the number of cell-centers in the
c               colatitude and longitude directions, respectively
c
c - -  Input:   a,b - the real*8 values of colatitude (theta) 
c               at the northern and southern edges of the problem boundary
c               [radians]
c
c - -  Input:   c,d - the real*8 values of longitude edges of problem boundary
c               [radians]
c
c - -  Input:   rsun - real*8 value for the radius of the Sun [km] 
c               Normally 6.96d5
c
c - -  Input:   en(n),es(n),el(m),er(m) - real*8 arrays of electric field values
c               multiplied by c (the speed of light)
c               at north, south, left and right edges of domain [G km/s].
c
c - -  Output:  elon(n,m+1) - real*8 array of electric fields multiplied
c               by the speed of light computed on theta edges (TE grid)
c               [G km/s].
c
c - -  Output:  elat(n+1,m) - real*8 array of electric fields multiplied
c               by c (the speed of light) computed on phi edges (PE grid)
c               [G km/s]..
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
c - - input variables:
c
      integer :: m,n
      real*8 :: en(n),es(n),el(m),er(m)
      real*8 :: rsun,a,b,c,d
c
c - - output variables:
c
      real*8 :: elon(n,m+1),elat(n+1,m)
c
c - - local subroutine variables:
c
      real*8 :: ea(n),eb(n),ec(m),ed(m)
      real*8 :: et(m,n+1),ep(n,m+1)
      integer :: i
c
c - - First, convert en,es,el,er to spherical polar coordinates:
c
      do i=1,m
         ec(i)=-el(m+1-i)
         ed(i)=-er(m+1-i)
      end do
      ea(1:n)=en(1:n)
      eb(1:n)=es(1:n)
c
c - - Use spherical polar Laplace solution to return et,ep
c
      call e_laplace_ss(m,n,a,b,c,d,rsun,ea,eb,ec,ed,et,ep)
c
c - - rotate solutions from theta,phi orientation to lon,lat orientation:
c
      call ehyeetp2ll_ss(m,n,et,ep,elon,elat)
c
c - - we're done!
c
      return
      end
