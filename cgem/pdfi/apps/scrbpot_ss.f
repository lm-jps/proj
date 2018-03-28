      subroutine scrbpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,brce,scrb3d)
c
c+ - Purpose: Compute 3-d potential field within a spherical wedge
c             domain, ie fixed co-latitude and azimuth boundaries, 
c             plus inner and outer radial shells at Rsun and
c             Rsun plus distance to source surface.  Solution is returned
c             as the poloidal potential scrb3d, from which all the magnetic
c             field components can be recovered, along with the vector
c             and scalar potentials.
c   - Usage:  call scrbpot_ss(m,n,p,a,b,c,d,rsun,rssmrs,brce,scrb3d)
c
c - - Input: m,n,p - number of cell interiors in the theta,phi, r directions
c - - Input: a,b - min, max values of theta (colatitude) at edges in 0,pi range 
c - - Input: c,d - min, max values of phi edges (longitude) in 0,2*pi range
c - - Input: rsun - radius of the Sun in assumed length units
c - - Input: rssmrs - distance from surface to source surface
c - - Input: brce - array of face-center interior values of B_r, 
c - -        at photosphere, dimensioned (m,n)
c - - Output: scrb3d - array of scriptb for potential field, 
c             dimensioned (m,n,p+1), at cell centers in theta,phi, 
c             and at edges in radius.
c - - Approach:  Solve the equation r^2 delh^2 scriptB + r^2 d^2/dr^2(scriptB)
c                =0 (Bercik's equation).  Convert finite difference in phi
c                within the horizontal Laplacian
c                to its Fourier representation, then solve equations in theta
c                and r for the amplitude of each Fourier mode, using the
c                blktri solution method in FISHPACK.
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
c - - Declarations of calling arguments:
c
      integer :: m,n,p
c
      real*8 :: a,b,c,d,rsun,rssmrs
      real*8 :: brce(m,n),scrb3d(m,n,p+1)
c
c - - Internal variables and arrays:
c
      integer :: i,j,q,iph,itmp,bcm,bcn,idimf,ierror
      integer :: kk,kl,mp,np,iflag
      real*8 :: dtheta,dphi,delr,dphifft,elm,pertrb,pi,dum
      real*8 :: brft(m,n),scrb(m,n),scrbft(m,n)
c     real*8 :: brtest(m,n),br(m,n)
      real*8 :: fce(m,n)
      real*8 :: k(n),fr(n)
      real*8 :: bdas(n),bdbs(n),bdcs(m),bdds(m)
      real*8 :: r(p+1),rce(p),sinth(m+1),sinth_hlf(m)
      real*8 :: am(m,n),bm(m,n),cm(m,n)
      real*8 :: an(p+1),bn(p+1),cn(p+1)
      real*8 :: y(m,p+1),scrbft3d(m,p+1,n)
      real*8, allocatable :: w(:),wb(:)
c
c - - function pimach in fishpack (fftpack) that computes pi
c
      real*8 :: pimach
c
c - - sdf variables (temporary)
c     integer*8 :: dims(20)
c
      elm=0.d0
      pi=pimach(dum)
c - - Define dtheta,dphi,dphifft,delr:
      dtheta=(b-a)/m
      dphi=(d-c)/n
c - - dphifft needed to FT in the phi direction
      dphifft=2.d0*pi/(n)
c - - radial spacing:
      delr=rssmrs/p
c
c - - define edge and cell-center radius arrays:
c
      r(1)=rsun
      do q=2,p+1
         r(q)=r(1)+(q-1)*delr
         rce(q-1)=rsun+0.5*delr+delr*(q-2)
      end do
c
c - - Compute sinth,sinth_hlf:
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - Compute the an,bn, and cn arrays needed in BLKTRI:
c - - These are the arrays that evaluate coefficients of r 2nd derivative
c - - terms in finite difference expansion
c
      do q=1,p+1
         if (q .eq. 1) then
c           At first radial point, no r-derivatives come into solution
c           (at q=1 BC set by r^2 del_h^2 scrb = -r^2 B_r(phot))
            an(q)=0.d0
            bn(q)=0.d0
            cn(q)=0.d0
         else if (q .eq. p+1) then
c           Homog Neumann BC applied half a zone above last radial point
c           (equivalent to source surface at 1/2 zone above last active point)
            an(q)=r(q)**2/delr**2
            bn(q)=-r(q)**2/delr**2
            cn(q)=0.d0
         else
c           Interior points in radius; regular 2nd derivative weightings
            an(q)=r(q)**2/delr**2
            bn(q)=-2.d0*r(q)**2/delr**2
            cn(q)=r(q)**2/delr**2
         endif
      end do   
c
c     Get fourier wave numbers stored into k array
      call kfft_ss(n,k)
c
c - - Compute the am,bm,cm arrays needed in BLKTRI:
c - - These are the coefficients of the finite difference expression
c - - for the horizontal Laplacian:
c
      do iph=1,m
         i=iph
         do j=1,n
            if(iph .eq. 1) then
c              Homog Neumann BC applied at small-theta edge
               am(iph,j)=0.d0
               bm(iph,j)=-sinth(i+1)/(sinth_hlf(iph)*dtheta**2)
               cm(iph,j)=sinth(i+1)/(sinth_hlf(iph)*dtheta**2)
            else if (iph .eq. m) then
c              Homog Neumann BC applied at large-theta edge
               am(iph,j)=sinth(i)/(sinth_hlf(iph)*dtheta**2)
               bm(iph,j)=-sinth(i)/(sinth_hlf(iph)*dtheta**2)
               cm(iph,j)=0.d0
            else
c              Interior points in theta
               am(iph,j)=sinth(i)/(sinth_hlf(iph)*dtheta**2)
               bm(iph,j)=-(sinth(i)+sinth(i+1))/
     1                  (sinth_hlf(iph)*dtheta**2)
               cm(iph,j)=sinth(i+1)/(sinth_hlf(iph)*dtheta**2)
            endif
c           add in contribution to Laplacian from phi terms in fourier space
c           to bm(iph,j)
            bm(iph,j)=bm(iph,j)-((2.*(1.-cos(k(j)*dphifft))
     1     /dphifft**2)*(((2.d0*pi)/(d-c))**2)
     1     /sinth_hlf(iph)**2)
         end do
      enddo
c
c - - Next order of business is to generate n Fourier modes in phi
c - - needed to represent the finite difference contribution to the
c - - horizontal Laplacian.  To do this, we will first solve the horizontal
c - - Poisson equation for scrb at the photosphere.
c - - Homogenous Neumann boundary conditions (theta edges)
c - - assumed, along with periodic boundary conditions (phi edges)
c - - in the Fishpack call,
c - - and then perform the horizonal Laplacian operator (in Fourier space)
c - - on scrb to generate the n Fourier modes that represent both B_r and
c - - scrb.
c
c - - compute size of work array, and allocate it
c
      itmp=(4*(n+2)+int(log(real(n+2))/log(2.0)+16)*(m+2))
      allocate(w(itmp))
c
c - - set boundary condition flags for Neumann (theta):
c    
      bcm=3
c
c - - boundary conditions periodic in phi
c
      bcn=0
c
c - - set homogenous Neumann boundary conditions for scrb:
c
      bdas(1:n)=0.d0
      bdbs(1:n)=0.d0
c - - bdcs and bdds will not actually be used in hstssp call because bcn=0
      bdcs(1:m)=0.d0
      bdds(1:m)=0.d0
c
c - - Set RHS for horiz. Poisson equation for scrb:
c
      fce(1:m,1:n)=-brce(1:m,1:n)*rsun**2
c
      idimf=m
c  
c - - Solve horiz. Poisson equation on photosphere for scrb on a staggered grid
c
      call hstssp(a,b,m,bcm,bdas,bdbs,c,d,n,bcn,bdcs,bdds,elm,fce,
     1           idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'scrbpot_ss: scrb ierror .ne. 0; = ',ierror
         stop
      endif
c
      deallocate(w)
c
c - - put solution points into scrb array:
c
      scrb(1:m,1:n)=fce(1:m,1:n)
c
c - - Now compute Fourier transform of scrb in phi direction:
c
c - - first allocate work array for Fourier transforms
c
      allocate(w(2*n+15))
c - - initialize fft work array
      call rffti(n,w)
c - transform scrb
      do i=1,m
         do j=1,n
           fr(j)=scrb(i,j)
         enddo 
         call rfftf(n,fr,w)
c        Divide by n since FFT in fftpack not normalized:
         scrbft(i,:)=fr(:)/(n)
      enddo
c
c - - Compute Fourier amplitude of -r^2 Br from r^2 grad_h^2 scrb = -r^2 B_r
c
      do j=1,n
         do i=1,m
c - - following if-test used to avoid out of range indexing
            if(i .eq. 1) then
               brft(i,j)=bm(i,j)*scrbft(i,j)+cm(i,j)*scrbft(i+1,j)
            else if (i .eq. m) then
               brft(i,j)=am(i,j)*scrbft(i-1,j)+bm(i,j)*scrbft(i,j)
            else
               brft(i,j)=am(i,j)*scrbft(i-1,j)+bm(i,j)*scrbft(i,j)+
     1         cm(i,j)*scrbft(i+1,j)
            endif
         end do
      end do
c - - Get Fourier transform of Br in phi direction as function of i:
c - - this code retained in comments but didn't work well.  Preserved
c - - in case it might be useful some other time.
c
c     do i=1,m
c        br(i,1:n)=-rsun**2*brce(i,1:n)
c        fr=br(i,1:n)
c        call rfftf(n,fr,w)
c        brft(i,1:n)=fr(1:n)/n
c     enddo
c
c
c - - Now, we're finally ready to do the solution to Bercik's equation
c - - in theta and r, using brft values at the photosphere for each Fourier
c - - mode:
c 
      kk=int(log(real(p+1)/log(real(2))))+1
      kl=2**(kk+1)
      itmp=(kk-2)*kl+5+max(2*(p+1),6*m)
c - - allocate work array for blktri:
      allocate(wb(itmp))
c     iflag=0 for blktri initialization
      iflag=0
c     np=1,mp=1 means non-periodic BC in radius and theta
      np=1
      mp=1
c
c - - initialization call for blktri
c
      y(1:m,1)=brft(1:m,1)
c RHS = 0 everwhere except the photosphere
      y(1:m,2:p+1)=0.d0
      call blktri(iflag,np,p+1,an,bn,cn,mp,m,am(1:m,1),bm(1:m,1),
     1 cm(1:m,1),m,y,ierror,wb)
      if(ierror .ne. 0) then
         write(6,*) 'scrbpot blktri init ierror .ne. 0, = ',ierror
         stop
      endif
c
c - - iflag = 1 means blktri has already been initialized
      iflag=1
c
c - - Outer loop over Fourier mode:
c
      do j=1,n
c        Set the RHS array for call to blktri:
         y(1:m,1)=brft(1:m,j)
c        RHS=0 except at photosphere
         y(1:m,2:p+1)=0.d0
         call blktri(iflag,np,p+1,an,bn,cn,mp,m,am(1:m,j),bm(1:m,j),
     1     cm(1:m,j),m,y,ierror,wb)
         if(ierror .ne. 0) then
            write(6,*) 'scrbpot blktri loop ierror .ne. 0, = ',ierror
            stop
         endif
c
c - - store results in 3d array
c
        scrbft3d(1:m,1:p+1,j)=y(1:m,1:p+1)
      end do
c
c - - end of blktri loop
c
      deallocate(wb)
c
c--------------------------------------------------------------------
c - - test values of Br in real space:
c - - this test retained as commented out text, in case future debugging
c - - needed
c
c     brtest(:,:)=brft(:,:)
c     do i=1,m
c        fr=brtest(i,1:n+2)
c        call rfftb(n+2,fr,w)
c        brtest(i,1:n+2)=fr
c     enddo
c     brtest(:,:)=-brtest(:,:)/rsun**2
c     dims(1)=m
c     dims(2)=n+2
c     call sdf_write_f77('scrbpot-debug.sdf','brtest','f',8,2,dims,
c    1 brtest)
c-------------------------------------------------------------------
c
c - - do inverse FT of scrbft3d into real space:
c      
      do q=1,p+1
         scrb3d(1:m,1:n,q)=scrbft3d(1:m,q,1:n)
      end do
      do q=1,p+1
         do i=1,m
            fr(1:n)=scrb3d(i,1:n,q)
            call rfftb(n,fr,w)
            scrb3d(i,1:n,q)=fr(1:n)
         end do
      end do
c
c - - we're done:
c
      deallocate(w)
      return
      end
