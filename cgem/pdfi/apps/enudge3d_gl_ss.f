      subroutine enudge3d_gl_ss(m,n,rsun,dr,btt,bpt,brt,ettop,etbot,
     1           eptop,epbot,er)
c
c+
c - - Purpose: Given time derivatives of the 3 magnetic field components
c              (btt,bpt,brt) through the corresponding voxel faces in a 
c              layer of voxels, computes the electric field on all edges 
c              of the spherical voxels.  Here, we assume global spherical 
c              geometry, rather than a spherical wedge domain.
c
c - - Usage:   call enudge3d_gl_ss(m,n,rsun,dr,btt,bpt,brt,ettop,etbot,eptop,
c              ebbot,er)
c
c - - Input:   m,n - integer values of number of cell centers in latitude,
c              and longitude, respectively.
c
c - - Input:   rsun - real*8 assumed value of the Sun's radius [km].
c
c - - Input:   dr - real*8 value of radial distance between top and bottom of
c              voxels [km].
c
c - - Input:   btt(m+1,n) - real*8 array of time derivatives of B_theta- TE grid
c              [G/sec]
c
c - - Input:   bpt(m,n+1) - real*8 array of time derivatives of B_phi  - PE grid
c              [G/sec]

c - - Input:   brt(m,n) - real*8 array of time derivatives of B_r - CE grid
c              [G/sec]
c
c - - Output:  ettop(m,n+1),etbot(m,n+1) - real*8 arrays of c*E_theta - PE grid 
c              along the top and bottom rails of the voxels, respectively.
c              [G km/sec]
c
c - - Output:  eptop(m+1,n),epbot(m+1,n) - real*8 arrays of c*E_phi - TE grid
c              along the top and bottom rails of the voxels, respectively.
c              [G km/sec]
c
c - - Output:  er(m+1,n+1) - real*8 array of c*E_r - COE grid - along the 
c              vertical rails of the voxels. [G km/sec]
c
c - - Note1:   The photospheric layer is assumed to lie halfway between the
c              top and bottom layer of voxels.  The top layer is at
c              R_S+0.5*dr, and the bottom layer is at R_s-0.5*dr.
c
c - - Note2:   FISHPACK's "global" boundary conditions assumed at the
c              N and S poles, and periodic boundary conditions assumed in 
c              longitude direction.
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
c - - variable declarations:
c
c - - calling arguments (input):
c
      integer :: m,n
      real*8 :: rsun,dr
      real*8 :: btt(m+1,n),bpt(m,n+1),brt(m,n)
c
c - - calling arguments (output):
c
      real*8 :: ettop(m,n+1),etbot(m,n+1),eptop(m+1,n),epbot(m+1,n)
      real*8 :: er(m+1,n+1)
c
c - - local variables:
c
      integer :: i,imh,iph,mp1,np1
      integer :: j,jph,jmh
      real*8 :: f(m,n),fj(m+1,n+1),divbh(m,n),curlbh(m-1,n-1)
      real*8 :: sinth(m+1),sinth_hlf(m),curl_lhs(m-1)
      real*8 :: bda(n), bdb(n), bdc(m), bdd(m) 
      real*8 :: bdaj(n+1), bdbj(n+1), bdcj(m+1), bddj(m+1) 
      real*8 :: ughost(m,n+2),uj(m+1,n+1)
      real*8 :: brtbal(m,n)
c - - Uncomment next statement if you want to create a version of uj with ghost
c - - zones
c     real*8 :: ujghost(m+1,n+3)
      real*8 :: et(m,n+1),etr(m,n+1),ep(m+1,n),epr(m+1,n)
      real*8, allocatable :: w(:)       
      integer :: idimf, ierror, itmp
      integer :: bcn, bcm
      real*8 :: a,b,c,d,elm,pertrb,pi,twopi,dum
      real*8 :: dtheta,dphi,rsun2
      real*8 :: sumbpnp,sumbpsp,curlbnp,curlbsp
      real*8 :: oneodp,oneodt
c
c - - sdf variables (for debugging):
c
c     integer*8 :: dims(20)
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
c - - Set a,b,c,d for global spherical geometry:
c
      a=0.d0
      b=pi
      c=0.d0
      d=twopi
c
c - - define dtheta, dphi:
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
      oneodp=1.d0/dphi
      oneodt=1.d0/dtheta
c
c - - compute sinth, sinth_hlf arrays from call to sinthta_ss:
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - Evaluate curlbh on interior corners:
c
      call curlh_co_ss(m,n,btt,bpt,rsun,sinth,sinth_hlf,dtheta,dphi,
     1     curlbh)
c
c - - Evaluate divbh on cell centers:
c
      call divh_ce_ss(m,n,btt,bpt,rsun,sinth,sinth_hlf,dtheta,dphi,
     1     divbh)
c
c - - for a global solution, need to make sure brt is flux balanced.
c - - call fluxbal_ss on brt input to provide brtbal (flux balanced) array:
c
      call fluxbal_ss(m,n,a,b,c,d,rsun,brt,brtbal)
c
c - - For RHS of Poisson equation for scrb, add a minus sign to brt and mpy 
c - - by rsun**2:
c
      f(:,:)=-brtbal(:,:)*rsun2
c
c - - The value of idimf should be exactly equal to m (staggered)
c
      idimf=m
c
c - - Global case:
c - - Set bcm=9, bcn=0, 
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
      itmp=(4*(n+2)+int(log(real(n+1))/log(2.0)+16)*(m+2))
c
c - - zero out ghost array:
c
      ughost(:,:)=0.d0
c
c - - allocate work array w
c
      allocate(w(itmp))
c
c - - Finally, solve Poisson's equation. Use hwscrt subroutine from Fishpack
c - - for Cartesian case (or hstcrt for staggered case);
c - - Use hwsssp for spherical case (hstssp for staggered spherical case):
c
      call hstssp (a,b,m,bcm,bda,bdb,c,d,n,bcn,bdc,bdd,elm,f,
     1      idimf,pertrb,ierror,w)
c
c - - check for errors returned from hstssp:
c
      if(ierror .ne. 0) then
         write(6,*) 'enudge3d_gl_ss(1): error in hstssp, ierror = ',
     1              ierror
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
c - - Now, we're going to do exactly the same set of operations, to find
c - - the curl of dscrbdr:
c
      f(:,:)=rsun2*divbh(:,:)
c
c - - zero out ghost array:
c
      ughost(:,:)=0.d0
c
c - - dscrbdr solution:
c - - Use hstssp for staggered spherical case:
c
      call hstssp (a,b,m,bcm,bda,bdb,c,d,n,bcn,bdc,bdd,elm,f,
     1      idimf,pertrb,ierror,w)
c
c - - check for errors returned from hstssp:
c
      if(ierror .ne. 0) then
         write(6,*) 'enudge3d_gl_ss(2): error in hstssp, ierror = ',
     1              ierror
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
c - - compute etr:
c
      do iph=1,m     
         do j=1,n+1
            jmh=j
            jph=j+1
            etr(iph,j)=-(ughost(iph,jph)-ughost(iph,jmh))/
     1                 (rsun*sinth_hlf(iph)*dphi)
         end do
      end do
c
c - - compute epr:
c
c
c - - Note that i=1 is North Pole, and i=m+1 is South Pole.  These
c - - are special cases, treated separately from main loop.  By
c - - assuming that dbrdt is finite at the poles, and using Faraday's Law
c - - for a small flux loop around the pole, can show that E_phi must go to zero
c - - as the north or south pole is approached.  Therefore, set ep at the poles
c - - to 0.  (This applies to radial derivatives of ep as well).
c
      epr(1,1:n)=0.d0
      epr(m+1,1:n)=0.d0
c
      do i=2,m
         imh=i-1
         iph=i
         do jph=1,n
            epr(i,jph)=(ughost(iph,jph+1)-ughost(imh,jph+1))/
     1      (rsun*dtheta)
         enddo
      enddo
c
c - - calculate ettop, etbot, eptop, epbot using modified code fragment 
c - - taken from evoxels3d_ss:
c
      ettop(1:m,1:n+1)=et(1:m,1:n+1)+0.5d0*dr*
     1     (etr(1:m,1:n+1)-et(1:m,1:n+1)/rsun)
      etbot(1:m,1:n+1)=et(1:m,1:n+1)-0.5d0*dr*
     1     (etr(1:m,1:n+1)-et(1:m,1:n+1)/rsun)
      eptop(1:m+1,1:n)=ep(1:m+1,1:n)+0.5d0*dr*
     1      (epr(1:m+1,1:n)-ep(1:m+1,1:n)/rsun)
      epbot(1:m+1,1:n)=ep(1:m+1,1:n)-0.5d0*dr*
     1      (epr(1:m+1,1:n)-ep(1:m+1,1:n)/rsun)
c
c - - Now solve the Poisson equation for scrj.  But first, need to
c - - find curl B_h at north and south poles.  First integrate bpt in longitude:
c
      sumbpnp=0.d0
      sumbpsp=0.d0
      do j=1,n
c       (integral of B_phi as a sum, using trapezoidal rule)
        sumbpnp=sumbpnp+0.5d0*(bpt(1,j)+bpt(1,j+1))
        sumbpsp=sumbpsp+0.5d0*(bpt(m,j)+bpt(m,j+1))
      enddo
c
c - - Now use these integrals to evaluate curl B_h at north and south poles:
c
      curlbnp=sumbpnp*dble(m)*4.d0/(dble(n)*pi*rsun)
      curlbsp=sumbpsp*dble(m)*4.d0/(dble(n)*pi*rsun)
c
c - - Evaluate curl at phi=c=0:  Following code fragment taken from
c - - subroutine curlh_co_ss, but evaluated at a single longitude index j=1
c - - (and corresponding jmh,jph) coinciding with the left boundary locations.
c
      j=1
      jph=1
c - - periodic boundary conditions: jmh=n
      jmh=n
      do i=1,m-1
        imh=i
        iph=i+1
        curl_lhs(i)=oneodt*(bpt(iph,j)*sinth_hlf(iph) - 
     1  bpt(imh,j)*sinth_hlf(imh))/sinth(i+1) -
     2  oneodp*(btt(i+1,jph)-btt(i+1,jmh))/sinth(i+1)
      enddo
c
c - - OK set F(:,:).  See HWSSSP documentation about setting RHS array for
c - - bcm=9, bcn=0.  First do interior active points, then boundary points:
c
      fj(2:m,2:n)=-curlbh(1:m-1,1:n-1)*rsun2
c
c - - curl at c=0, d=2*pi, (left and right boundaries) 
c     for i=2,m (and j=1 and j=n+1), (bcn=0):
c
      fj(2:m,1)=-curl_lhs(1:m-1)*rsun2
      fj(2:m,n+1)=-curl_lhs(1:m-1)*rsun2
c
c - - north and south poles, for j=1,n+1, and for i=1 and i=m+1, (bcm=9):
c
      fj(1,1:n+1)=-curlbnp*rsun2
      fj(m+1,1:n+1)=-curlbsp*rsun2
c
c - - boundary condition arrays should be ignored, but set them to
c - - 0 anyway:
c
      bdaj(1:n+1)=0.d0
      bdbj(1:n+1)=0.d0
      bdcj(1:m+1)=0.d0
      bddj(1:m+1)=0.d0
c
c - - for hwsssp, need to set idimf to m+1:
c
      idimf=m+1
c
c - - Poisson equation solution for scrj, non-staggered grid on COE gridpoints
c
      call hwsssp(a,b,m,bcm,bdaj,bdbj,c,d,n,bcn,bdcj,bddj,elm,
     1 fj,idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'enudge3d_gl_ss: error in hwsssp call, ierror = ',
     1              ierror
         stop
      endif
c
c - - Set uj to RHS array containing solution:
c
      uj(1:m+1,1:n+1)=fj(1:m+1,1:n+1)
      er(1:m+1,1:n+1)=-uj(1:m+1,1:n+1)
c
c - - If you wanted to set ujghost array to add ghost zones on left and right
c - - to conveniently set periodic
c - - boundary conditions, this is where you would do that
c
c     ujghost(1:m+1,2:n+2)=uj(1:m+1,1:n+1)
c     ujghost(1:m+1,1)=ughost(1:m+1,n+2)
c     ujghost(1:m+1,n+3)=ughost(1:m+1,2)
c
c - - deallocate the work array:
c
      deallocate(w)
c
c - - we're done!
c
      return
      end
