      subroutine enudge_gl_ss(m,n,rsun,brt,et,ep)
c
c+
c - - Purpose: Compute horizontal electric field components for a global
c              spherical solar model on a staggered mesh.
c
c - - Usage:   call enudge_gl_ss(m,n,rsun,brt,et,ep)
c
c - - Input:   m,n - integer values of number of cell centers in latitude,
c              and longitude, respectively.
c
c - - Input:   rsun - real*8 value of the Sun's radius [km].  Normally 6.96d5.
c
c - - Input:   brt(m,n) - real*8 array of time derivatives of the radial
c              magnetic field component, evaluated at cell-centers (CE grid). 
c              brt is assumed to be ordered in colat,lon order on input.
c              [G/sec]
c
c - - Output:  et(m,n+1), ep(m+1,n) - real*8 arrays of the theta and phi
c              components of the electric field (multiplied by c), respectively,
c              on the edges surrounding each cell-center.  et is on the PE
c              grid, ep on the TE grid.  Arrays are in colat,lon order
c              (rather than lon,lat order). [G km/sec]
c
c - - Note:    FISHPACK's "global" boundary conditions assumed at the
c              N and S poles, and periodic boundary conditions assumed in 
c              longitude direction.
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
c
      implicit none
c
c - - variable declarations:
c
c - - calling arguments (input):
c
      integer :: m,n
      real*8 :: rsun
      real*8 :: brt(m,n)
c
c - - calling arguments (output):
c
      real*8 :: et(m,n+1),ep(m+1,n)
c
c - - local variables:
c
      integer :: i,iph,j,mp1,np1
      integer :: imh,jph,jmh
      real*8 :: f(m,n)
      real*8 :: sinth(m+1),sinth_hlf(m)
      real*8 :: bda(n), bdb(n), bdc(m), bdd(m) 
      real*8 :: ughost(m,n+2)
      real*8, allocatable :: w(:)       
      integer :: idimf, ierror, itmp
      integer :: bcn, bcm
      real*8 :: a,b,c,d,elm,pertrb,pi,twopi,dum
      real*8 :: dtheta,dphi,rsun2
      real*8 :: brtbal(m,n)
c
c - - pimach function declaration:
c
      real*8 :: pimach
c
c - - set initial defaults:
c
      pi=pimach(dum)
      twopi=2.*pimach(dum)
      rsun2=rsun**2
      mp1=m+1
      np1=n+1
c
c - - The value of idimf should be exactly equal to m (staggered)
c
      idimf=m
c
c - - Global case:  a=0.d0, b=pi; c=0.d0, d=2 pi.
c - - Set bcm=9, bcn=0
c
      a=0.d0
      b=pi
      dtheta=(b-a)/dble(m)
      c=0.d0
      d=twopi
      dphi=(d-c)/dble(n)
c
c - - remove monopole term from input brt array:
c
      call fluxbal_ss(m,n,a,b,c,d,rsun,brt,brtbal)
c
c - - For RHS of Poisson equation, add a minus sign to brtbal, mpy by rsun**2:
c
      f(:,:)=-brtbal(:,:)*rsun2
c
c - - compute sinth, sinth_hlf arrays from call to sinthta_ss:
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c 
c - - set bcm (boundary condition type at theta=a, theta=b) for Global:
c
      bcm=9
c
c - - set bcn (boundary condition at phi=c, phi=d) for periodic:
c
      bcn=0
c
c - - set bda,bdb to 0 (but don't think they are used):
c
      bda(1:n)=0.d0
      bdb(1:n)=0.d0
c
c - - set bdc,bdd to zero but not sure they are used.
c
      bdc(1:m)=0.d0
      bdd(1:m)=0.d0
c
c - - set elm (coefficient for Helmholtz term) to 0.
c
      elm=0.d0
c
c - - compute the dimension itmp for the work array w;
c - - m and n represent the number of cell interiors in t,p directions, resp.
c - - and are 1 less than mp1,np1, resp.
c
      itmp=(4*(n+1)+int(log(real(n+1))/log(2.0)+16)*(m+1))
c
c - - allocate work array w
c
      allocate(w(itmp))
c
c - - Finally, solve Poisson's equation! Use hwscrt subroutine from Fishpack
c - - for Cartesian case (or hstcrt for staggered case);
c - - Use hwsssp for spherical case (hstssp for staggered spherical case):
c
      call hstssp (a,b,m,bcm,bda,bdb,c,d,n,bcn,bdc,bdd,elm,f,
     1      idimf,pertrb,ierror,w)
c
c - - check for errors returned from hstssp:
c
      if(ierror .ne. 0) then
         write(6,*) 'enudge_gl_ss: error in hstssp, ierror = ',ierror
         stop
      endif
c
c - - set solution ughost to source term array (where hstssp puts the solution)
c
      ughost(1:m,2:n+1)=f(1:m,1:n)
c
c - - apply Periodic boundary values to fill ghost zones on phi edges: 
c
      ughost(1:m,1)=ughost(1:m,n+1)
      ughost(1:m,n+2)=ughost(1:m,2)
c
c - - compute et:
c
      do iph=1,m     
         do j=1,n+1
            jmh=j
            jph=j+1
            et(iph,j)=-(ughost(iph,jph)-ughost(iph,jmh))/
     1                 (rsun*sinth_hlf(iph)*dphi)
         end do
      end do
c
c - - compute ep:
c
c
c - - Note that i=1 is North Pole, and i=m+1 is South Pole.  These
c - - are special cases, treated separately from main loop.  By
c - - assuming that dbrdt is finite at the poles, and using Faraday's Law
c - - for a small flux loop around the pole, can show that E_phi must go to zero
c - - as the north or south pole is approached.  Therefore, set ep at the poles
c - - to 0.
c
      ep(1,1:n)=0.d0
      ep(m+1,1:n)=0.d0
c
      do i=2,m
         imh=i-1
         iph=i
         do jph=1,n
            ep(i,jph)=(ughost(iph,jph+1)-ughost(imh,jph+1))/
     1      (rsun*dtheta)
         enddo
      enddo
c
c - - deallocate the work array:
c
      deallocate(w)
c
c - - we're done!
c
      return
      end
