      subroutine ptdsolve_ss(m,n,bt,bp,br,rsun,sinth,sinth_hlf,
     1 a,b,c,d,scrb,dscrbdr,scrj)
c
c - - documentation below can be seen using DOC_LIBRARY,'ptdsolve_ss'
c+
c - - Purpose:  This procedure solves the 2D PTD Poisson equations for 
c               scrb, dscrbdr, and scrj 
c               (or their time derivatives) from
c               the theta, phi, and radial components of the magnetic field 
c               (or the time derivatives of bt,bp, and br).
c
c - - Usage:  call ptdsolve_ss(m,n,bt,bp,br,rsun,sinth,sinth_hlf,a,b,c,d,
c             scrb, dscrbdr,scrj)
c
c - - Input:  m,n - number of cell-centers in the theta, phi directions, resp.
c - - Input:  bt - the component of the magnetic field pointing in the
c - - colatitude direction (positive is southward pointing).  The array
c - - is assumed to be indexed i,j, where i denotes the colatitude direction
c - - with colatitude increasing as i increases; j denotes the azimuthal
c - - (longitudinal) direction, with azimuthal angle increasing as j increases.
c - - bt dimensioned m+1,n and is defined on theta edges, midway between
c - - phi edges.
c
c - - Input: bp - array of the azimuthal component of the magnetic field, 
c - - indexed as described for bt, but bp dimensioned m,n+1, and is defined on
c - - phi edges, mid-way between theta edges.
c
c - - Input: br - array of the radial component of the magnetic field,
c - - indexed as described for bt, bp, but br is defined at cell centers.
c - - br is dimensioned m,n.
c
c - - Input: rsun - units for the radius of the Sun.
c
c - - Input: sinth, the value of sin(colatitude), computed for theta edge
c - - locations.  Is dimensioned m+1.
c 
c - - Input: sinth_hlf, the value of sin(colatitude) computed for cell-center
c - - values of theta.  Is dimensioned m.
c
c - - a,b:  Minimum and maximum values of co-latitude in radians
c - - corresponding to the
c - - range of the theta (colatitude) edge values.
c
c - - c,d:  Minimum and maximum values of longitude in radians corresponding
c - - to the range of longitude (azimuth) edge values.
c
c - - Output: scrb - the scribtb array, located at cell-centers,
c - - dimensioned m+2,n+2, ie 2 larger in the theta,phi
c - - directions than the input br array.  It includes  ghostzones.
c
c - - Output: dscrbdr - the dscriptb/dr array, dimensioned as for scrb.
c - - Also located at cell-centers, and also includes ghostzones,
c - - dimensioned m+2,n+2
c
c - - Output: scrj - the scriptj array, dimensioned m+1,n+1, and is located
c - - at cell corners.  Corner values on the problem boundary are set to
c - - 0 using Dirichlet boundary conditions.
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
c
      implicit none
c
      integer :: m,n,mp1,np1
c
c - - calling argument declarations:
c
      real*8 :: a,b,c,d,rsun
      real*8 :: bt(m+1,n),bp(m,n+1),br(m,n),sinth(m+1),sinth_hlf(m),
     1 scrb(m+2,n+2),dscrbdr(m+2,n+2),scrj(m+1,n+1)
c
c - - local allocatable work array:
c
      real*8, allocatable :: w(:)
c
c - - local variables to ptdsolve
c
      integer :: idimf,itmp,bcm,bcn,ierror
      real*8 :: dtheta,dphi,elm,pertrb
      real*8 :: curlbh(m-1,n-1),divbh(m,n),fce(m,n),fcoe(m+1,n+1),
     1 bdas(n),bdbs(n),bdcs(m),bdds(m),bdad(n),bdbd(n),
     2 bdcd(m),bddd(m),bdaj(n+1),bdbj(n+1),bdcj(m+1),bddj(m+1)
c
      elm=0.d0
      mp1=m+1
      np1=n+1
      dtheta=(b-a)/m
      dphi=(d-c)/n
c
c - - set boundary condition flags for Neumann:
c    
      bcm=3
      bcn=3
c
c - - compute curlbh on cell corners (m-1,n-1)
c
      call curlh_co_ss(m,n,bt,bp,rsun,sinth,sinth_hlf,dtheta,
     1 dphi,curlbh)
c
c - - compute divbh on cell centers (m,n)
c
      call divh_ce_ss(m,n,bt,bp,rsun,sinth,sinth_hlf,dtheta,
     1 dphi,divbh)
c
c - - set homogenous Neumann boundary conditions for scrb:
c
      bdas(1:n)=0.d0
      bdbs(1:n)=0.d0
      bdcs(1:m)=0.d0
      bdds(1:m)=0.d0
c
c - - Set RHS for Poisson equation for scrb:
c
      fce(1:m,1:n)=-br(1:m,1:n)*rsun**2
c
c - - compute size of work array, and allocate it
c
      itmp=(4*(n+2)+int(log(real(n+2))/log(2.0)+16)*(m+2))
      allocate(w(itmp))
      idimf=m
c  
c - - Solve Poisson equation for scrb on a staggered grid
c
      call hstssp(a,b,m,bcm,bdas,bdbs,c,d,n,bcn,bdcs,bdds,elm,fce,
     1           idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'ptdsolve_ss: scrb ierror .ne. 0; = ',ierror
         stop
      endif
c
c - - put solution points into scrb array:
c
      scrb(2:m+1,2:n+1)=fce(1:m,1:n)
c
c - - Now fill in ghost zone values using boundary condition arrays:
c
      scrb(1,2:n+1)=scrb(2,2:n+1)-1.d0*dtheta*bdas(1:n)
      scrb(m+2,2:n+1)=scrb(m+1,2:n+1)+1.d0*dtheta*bdbs(1:n)
      scrb(2:m+1,1)=scrb(2:m+1,2)-1.d0*dphi*bdcs(1:m)
      scrb(2:m+1,n+2)=scrb(2:m+1,n+1)+1.d0*dphi*bdds(1:m)
c
c - - fill in ghost zone corners (not necessary, but cosmetic)
c
      scrb(1,1)=0.5d0*(scrb(1,2)+scrb(2,1))
      scrb(m+2,1)=0.5d0*(scrb(m+2,2)+scrb(m+1,1))
      scrb(1,n+2)=0.5d0*(scrb(2,n+2)+scrb(1,n+1))
      scrb(m+2,n+2)=0.5d0*(scrb(m+2,n+1)+scrb(m+1,n+2))
c
c - - Move onto the solution for dscrbdr now:
c - - First compute Neumann boundary condtions at colat, lon edges:
c
c - - Longitude edges:
      bdcd(1:m)=rsun*sinth_hlf(1:m)*bp(1:m,1)
      bddd(1:m)=rsun*sinth_hlf(1:m)*bp(1:m,n+1)
c - - Co-latitude edges:
      bdad(1:n)=rsun*bt(1,1:n)
      bdbd(1:n)=rsun*bt(m+1,1:n)
c
c - - Set RHS for dscrbdr Poisson equation:
c
      fce(1:m,1:n)=divbh(1:m,1:n)*rsun**2
c
c - - Solve Poisson equation for dscrbdr on a staggered grid:
c
      call hstssp(a,b,m,bcm,bdad,bdbd,c,d,n,bcn,bdcd,bddd,elm,fce,
     1           idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'ptdsolve_ss: dscrbdr ierror .ne. 0; = ',ierror
         stop
      endif
c
c - - put solution points into dscrbdr array:
c
      dscrbdr(2:m+1,2:n+1)=fce(1:m,1:n)
c
c - - Now fill in ghost zone values using boundary condition arrays:
c
      dscrbdr(1,2:n+1)=dscrbdr(2,2:n+1)-1.d0*dtheta*bdad(1:n)
      dscrbdr(m+2,2:n+1)=dscrbdr(m+1,2:n+1)+1.d0*dtheta*bdbd(1:n)
      dscrbdr(2:m+1,1)=dscrbdr(2:m+1,2)-1.d0*dphi*bdcd(1:m)
      dscrbdr(2:m+1,n+2)=dscrbdr(2:m+1,n+1)+1.d0*dphi*bddd(1:m)
c
c - - fill in ghost zone corners (not necessary, but cosmetic)
c
      dscrbdr(1,1)=0.5d0*(dscrbdr(1,2)+dscrbdr(2,1))
      dscrbdr(m+2,1)=0.5d0*(dscrbdr(m+2,2)+dscrbdr(m+1,1))
      dscrbdr(1,n+2)=0.5d0*(dscrbdr(2,n+2)+dscrbdr(1,n+1))
      dscrbdr(m+2,n+2)=0.5d0*(dscrbdr(m+2,n+1)+dscrbdr(m+1,n+2))
c
c - - Poisson equation for scrj, using homog Dirichlet BC, non-staggered grid
c
c
c - - Boundary conditions flag for Dirichlet:
c
      bcm=1
      bcn=1
      idimf=mp1
c
c - - Set RHS for scriptj equation:
c
      fcoe(2:m,2:n)=-curlbh(1:m-1,1:n-1)*rsun**2
c
c - - Set up Homog. Dirichlet BC arrays:
c
      bdaj(1:n+1)=0.d0
      bdbj(1:n+1)=0.d0
      bdcj(1:m+1)=0.d0
      bddj(1:m+1)=0.d0
c
c - - set boundary values within the f array:
c
      fcoe(1,1:n+1)=bdaj(1:n+1)
      fcoe(m+1,1:n+1)=bdbj(1:n+1)
      fcoe(1:m+1,1)=bdcj(1:m+1)
      fcoe(1:m+1,n+1)=bddj(1:m+1)
c
c - - Solve Poisson equation for scrj, non-staggered grid: 
c
      call hwsssp(a,b,m,bcm,bdaj,bdbj,c,d,n,bcn,bdcj,bddj,elm,
     1 fcoe,idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'ptdsolve_ss: scrj ierror .ne. 0; = ',ierror
         write(6,*) 'ptdsolve_ss: c,d = ',c,d
         stop
      endif
c
c - - put solution into scrj array
c
      scrj(1:m+1,1:n+1)=fcoe(1:m+1,1:n+1)

c
c - - we're done:
c
      deallocate(w)
c
      return
      end
