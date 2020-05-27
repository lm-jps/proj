      subroutine enudge_ss(m,n,a,b,c,d,rsun,brt,et,ep)
c+
c - -  Purpose:  Compute horizontal electric fields et,ep from radial magnetic 
c                field time derivatives brt, using the simplest PTD electric 
c                field solution. Assumes for boundary conditions that the 
c                the net flux time derivative is accounted for by a constant
c                amplitude electric field on all 4 boundaries.
c
c - -  Usage:    call enudge_ss(m,n,a,b,c,d,rsun,brt,et,ep)
c
c - -  Input:    m,n - integers denoting the number of cell-centers in the
c                colatitude and longitude directions, respectively
c
c - -  Input:    a,b - the real*8 values of colatitude (theta) 
c                at the northern and southern edges of the problem boundary
c                [radians]
c
c - -  Input:    c,d - the real*8 values of longitude edges [radians]
c
c - -  Input:    rsun - real*8 value of radius of Sun [km]. Normally 6.96d5.
c
c - -  Input:    brt(m,n) - real*8 array of cell-center (CE grid) 
c                radial magnetic field time derivative values [G/sec]
c
c - -  Output:   et(m,n+1) - real*8 array of electric fields, multiplied by 
c                the speed of light, computed on phi edges [G km/s].
c
c - -  Output:   ep(m+1,n) - real*8 array of electric fields, multiplied by 
c                the speed of light, computed on theta edges [G km/s].
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
c - - variable declarations:
c
c - - input arguments:
c
      integer :: m,n
      real*8 :: brt(m,n)
      real*8 :: rsun,a,b,c,d
c
c - - output arguments:
c
      real*8 :: et(m,n+1),ep(m+1,n)
c
c - - local subroutine variables:
c
      integer :: i,iph,j,mp1,np1
      integer :: imh,jph,jmh
      real*8 :: f(m,n) 
      real*8 :: sinth(m+1),sinth_hlf(m)
      real*8 :: bda(n), bdb(n), bdc(m), bdd(m) 
      real*8 :: ughost(m+2,n+2)
      real*8, allocatable :: w(:)       
      integer :: idimf, ierror, itmp
      integer :: bcn, bcm
      real*8 :: elm,pertrb,pi,twopi,dum
      real*8 :: dtheta,dphi,rsun2,cc,dd
      real*8 :: flux,lperim,eperim
c
c - - function for setting pi in fishpack:
c
      real*8 :: pimach
c
c - - set initial defaults:
c
      pi=pimach(dum)
      twopi=2.*pimach(dum)
      rsun2=rsun**2
c
      mp1=m+1
      np1=n+1
c
c - - Set RHS for Poisson solve:
c
      f(1:m,1:n)=-brt(1:m,1:n)*rsun2
c
c - - The value of idimf should be exactly equal to m (staggered), or mp1 
c - - (non-staggered)
c
      idimf=m
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
c
      cc=0.d0
      dd=dble(n)*dphi
c
c - - assign theta and sin(theta) arrays
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - Compute area integral of dB_r/dt: (here assumed = br).  Assumed non-zero:
c
      flux=0.d0
      do j=1,n
         do i=1,m
            flux=flux+brt(i,j)*sinth_hlf(i)*dtheta*dphi*rsun**2
         enddo
      enddo
c
c - - compute perimeter length of spherical wedge: (now ignoring l and r edges)
c
c     lperim=rsun*((d-c)*(sin(a)+sin(b))+2.d0*(b-a))
      lperim=rsun*((d-c)*(sin(a)+sin(b)))
c
c - - compute amplitude of uniform electric field on boundary that generates
c - - flux imbalance:
c
      eperim = -flux/lperim
c

c 
c - - set bcm (boundary condition type at theta=a, theta=b) for Neumann:
c
      bcm=3
c
c - - set bcn (boundary condition at phi=c, phi=d) for Neumann:
c
      bcn=3
c
c - - set Neumann boundary conditions for scrb based on E on perimeter:
c - - (recall E = -curl scrb rhat)
c
c - - E_phi at theta=a = -eperim
      bda(1:n)= -eperim*rsun
c - - E_phi at theta=b = eperim
      bdb(1:n)= eperim*rsun
c - - E_theta at phi=c = eperim (now zero)
c     bdc(1:m)= -eperim*rsun*sinth_hlf(1:m)
      bdc(1:m)= 0.d0
c - - E_theta at phi=d = -eperim (now zero)
c     bdd(1:m)=eperim*rsun*sinth_hlf(1:m)      
      bdd(1:m)= 0.d0
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
c - - Finally, solve Poisson's equation. Use hstssp subroutine from Fishpack
c - - for staggered case:
c
      call hstssp (a,b,m,bcm,bda,bdb,cc,dd,n,bcn,bdc,bdd,elm,f,
     1      idimf,pertrb,ierror,w)
c
c - - check for errors returned from hstssp:
c
      if(ierror .ne. 0) then
         write(6,*) 'enudge_ss: error in hstssp, ierror = ',ierror
         stop
      endif
c
c - - allocate the solution arrays u, ughost, and
c - - electric field arrays et and ep, and the angle arrays:
c
c - - set solution u to source term array (where hstssp puts the solution)
c
      ughost(2:m+1,2:n+1)=f(1:m,1:n)
c
c - - apply Neumann boundary values to fill ghost zones on edges: 
c
      ughost(1,2:n+1)=ughost(2,2:n+1)-1.d0*dtheta*bda(1:n)
      ughost(m+2,2:n+1)=ughost(m+1,2:n+1)+1.d0*dtheta*bdb(1:n)
      ughost(2:m+1,1)=ughost(2:m+1,2)-1.d0*dphi*bdc(1:m)
      ughost(2:m+1,n+2)=ughost(2:m+1,n+1)+1.d0*dphi*bdd(1:m)
c
c - - fill in ghost zone corners (not necessary, but cosmetic)
c
      ughost(1,1)=0.5d0*(ughost(1,2)+ughost(2,1))
      ughost(m+2,1)=0.5d0*(ughost(m+2,2)+ughost(m+1,1))
      ughost(1,n+2)=0.5d0*(ughost(2,n+2)+ughost(1,n+1))
      ughost(m+2,n+2)=0.5d0*(ughost(m+2,n+1)+ughost(m+1,n+2))
c
c - - compute et:
c
      do iph=1,m     
         do j=1,n+1
            jmh=j
            jph=j+1
            et(iph,j)=-(ughost(iph+1,jph)-ughost(iph+1,jmh))/
     1                 (rsun*sinth_hlf(iph)*dphi)
         end do
      end do
c
c - - compute ep:
c
      do i=1,m+1
         imh=i
         iph=i+1
         do jph=1,n
            ep(i,jph)=(ughost(iph,jph+1)-ughost(imh,jph+1))/
     1      (rsun*dtheta)
         enddo
      enddo
c
c - - deallocate the allocatable array w:
c
      deallocate(w)
c
c - - we're done!
c
      return
      end
