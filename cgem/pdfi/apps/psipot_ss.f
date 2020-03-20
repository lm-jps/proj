      subroutine psipot_ss(m,n,p,bcn,a,b,c,d,rsun,rssmrs,btte,bppe,
     1           psi3d)
c
c+ - Purpose: Compute 3-d potential field within a spherical wedge
c             domain, ie fixed co-latitude and azimuth boundaries, 
c             plus inner and outer radial shells at Rsun and
c             Rsun plus distance to source surface.  Solution is returned
c             as the scalar potential psi3d, from which all the magnetic
c             field components can be recovered. Solution returned along
c             radial faces.  At photosphere, solution matches the condition
c             r^2 delh^2 psi = -r^2 div_h (B_h).
c
c - - Method: Solve the equation r^2 delh^2 psi + d/dr(r^2 d psi/dr)=0
c             (Laplace equation).  Convert finite difference in phi
c             within the horizontal Laplacian
c             to its Fourier representation, then solve equations in theta
c             and r for the amplitude of each Fourier mode, using the
c             blktri solution method in FISHPACK.
c
c   - Usage:  call psipot_ss(m,n,p,bcn,a,b,c,d,rsun,rssmrs,btte,bppe,psi3d)
c
c - - Input:  m,n,p - integer number of cell interiors in the theta,phi,r 
c             directions, respectively.
c
c - - Input:  bcn - integer flag for boundary conditions in phi
c             (=0 for periodic BC, =3 for homogenous Neumann BC)
c
c - - Input:  a,b - real*8 min, max values of theta (colatitude) at edges 
c             in 0,pi range [radians]
c
c - - Input:  c,d - real*8 min, max values of phi edges (longitude) 
c             in 0,2*pi range [radians]
c
c - - Input:  rsun - real*8 value of radius of Sun [km].  Normally 6.96d5.
c
c - - Input:  rssmrs - real*8 distance from surface to source surface [km]
c
c - - Input:  btte(m+1,n) - real*8 array of theta-edge (TE grid) values of B_t, 
c             at photosphere, in colat,lon index order. [G]
c
c - - Input:  bppe(m,n+1) - real*8 array of phi-edge (PE grid) values of B_p, 
c             at photosphere, in colat,lon index order. [G]
c
c - - Output: psi3d(m,n,p+1) - real*8 3D array of psi (scalar potential)
c             for a potential magnetic field, evaluated at at cell centers in 
c             theta,phi, and at radial shells. [G km]
c
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
c - - input variables:
c
      integer :: m,n,p,bcn
      real*8 :: a,b,c,d,rsun,rssmrs
      real*8 :: btte(m+1,n),bppe(m,n+1)
c
c - - output variables:
c
      real*8 :: psi3d(m,n,p+1)
c
c - - Internal variables and arrays:
c
      integer :: i,j,q,iph,itmp,bcm,idimf,ierror
      integer :: kk,kl,mp,np,iflag
      real*8 :: dtheta,dphi,delr,dphifft,elm,pertrb,pi,dum
      real*8 :: rangecos,mflux
c     real*8 :: alpha
      real*8 :: rqmh2,rqph2
c     real*8 :: psitot,psibar
c     real*8 :: scrbtot,area
      real*8 :: brft(m,n),scrb(m,n),scrbft(m,n),divbh(m,n)
c     real*8 :: brtest(m,n),br(m,n)
      real*8 :: fce(m,n)
      real*8 :: k(n),fr(n)
      real*8 :: bdas(n),bdbs(n),bdcs(m),bdds(m)
      real*8 :: r(p+1),rce(p),sinth(m+1),sinth_hlf(m)
      real*8 :: am(m,n),bm(m,n),cm(m,n)
      real*8 :: an(p+1),bn(p+1),cn(p+1)
      real*8 :: y(m,p+1),psift3d(m,p+1,n)
      real*8, allocatable :: w(:),wb(:)
c
c - - function pimach in fishpack (fftpack) that computes pi
c
      real*8 :: pimach
c
c - - sdf variables (temporary)
c     integer*8 :: dims(20)
c
      if((bcn .ne. 0) .and. (bcn .ne. 3)) then
         write(6,*) 'psipot_ss: Illegal bcn = ',bcn, ' exiting'
         stop
      endif
c
      elm=0.d0
      pi=pimach(dum)
c - - Define dtheta,dphi,dphifft,delr:
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
c
c - - dphifft needed to FT in the phi direction
c
      if (bcn .eq. 3) then
         dphifft=pi/dble(n)
      else
         dphifft=2.d0*pi/dble(n)
      endif
c
      
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
c     CHECK THESE EQNS!!!!
c
      do q=1,p+1
         rqph2=(r(q)+0.5d0*delr)**2
         rqmh2=(r(q)-0.5d0*delr)**2
         if (q .eq. 1) then
c
c           At first radial point, (photosphere) no r-derivatives come 
c           into solution
c           (at q=1 BC set by r^2 del_h^2 scrb = -r^2 div B_h)
c
            an(q)=0.d0
            bn(q)=0.d0
            cn(q)=0.d0
         else if (q .eq. p+1) then
c
c - - at source surface, we demand psi-> 0, meaning psi(r_[p+2]=-psi(r_[p+1])
c - - which results in bn(q) being weighted by -3 instead of -2, and cn=0.
c
            an(q)=rqmh2/delr**2
            bn(q)=-(2.d0*rqph2+rqmh2)/delr**2
            cn(q)=0.d0
         else
c
c           Interior points in radius; 
c
            an(q)=rqmh2/delr**2
            bn(q)=-(rqmh2+rqph2)/delr**2
            cn(q)=rqph2/delr**2
         endif
      end do   
c
c     Get fourier wave numbers stored into k array
c
      if(bcn .eq. 3) then
        call kcost_ss(n,k)
      else
        call kfft_ss(n,k)
      endif
c
c - - Compute the am,bm,cm arrays needed in BLKTRI:
c - - These are the coefficients of the finite difference expression
c - - for the horizontal Laplacian:
c
      if(bcn .eq. 3) then
        rangecos=pi
      else
        rangecos=2.d0*pi
      endif
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
     1     /dphifft**2)*(((rangecos)/(d-c))**2)
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
c - - boundary conditions periodic in phi (= bcn) is read in as calling arg
c
c - - set homogenous Neumann boundary conditions for scrb (but if bcn=0
c - - these boundary conditions for phi=c and phi=d ignored):
c
      bdas(1:n)=0.d0
      bdbs(1:n)=0.d0
      bdcs(1:m)=0.d0
      bdds(1:m)=0.d0
c
c - - Set RHS for horiz. Poisson equation for psi:
c
      call divh_ce_ss(m,n,btte,bppe,rsun,sinth,sinth_hlf,dtheta,dphi,
     1     divbh)
      fce(1:m,1:n)=-divbh(1:m,1:n)*rsun**2
c
      idimf=m
c  
c - - Solve horiz. Poisson equation on photosphere for scrb on a staggered grid
c
      call hstssp(a,b,m,bcm,bdas,bdbs,c,d,n,bcn,bdcs,bdds,elm,fce,
     1           idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'psipot_ss: scrb ierror .ne. 0; = ',ierror
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
c - - (big enough for both regular fft and cost transforms)
c
      allocate(w(3*n+15))
c - - initialize fft work array
c
      if(bcn .eq. 3) then
        call costi(n,w)
      else
        call rffti(n,w)
      endif
c
c - transform scrb
      do i=1,m
         do j=1,n
           fr(j)=scrb(i,j)
         enddo 
c
         if (bcn .eq. 3) then
            call cost(n,fr,w)
         else
            call rfftf(n,fr,w)
         endif
c
c        Divide by n since FFT in fftpack not normalized:
c
         if(bcn .eq. 3) then
           scrbft(i,:)=fr(:)/dble(2*(n-1))
         else
           scrbft(i,:)=fr(:)/dble(n)
         endif
c
      enddo
c
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
c
c - - Now, we're finally ready to do the solution to Laplace equation
c - - in theta and r, using brft values at the photosphere for each Fourier
c - - mode:
c 
      kk=int(log(real(p+1))/log(real(2)))+2
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
      if(wb(1) .gt. dble(itmp)) then
         write(6,*) 'psipot_ss: itmp too small for wb, exiting'
         stop
      endif
      if(ierror .ne. 0) then
         write(6,*) 'psipot blktri init ierror .ne. 0, = ',ierror
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
            write(6,*) 'psipot blktri loop ierror .ne. 0, = ',ierror
            stop
         endif
c
c - - store results in 3d array
c
        psift3d(1:m,1:p+1,j)=y(1:m,1:p+1)
      end do
c
c - - end of blktri loop
c
      deallocate(wb)
c
c - - do inverse FT of scrbft3d into real space:
c      
      do q=1,p+1
         psi3d(1:m,1:n,q)=psift3d(1:m,q,1:n)
      end do
      do q=1,p+1
         do i=1,m
            fr(1:n)=psi3d(i,1:n,q)
c
c - - guts of inverse transform:
c
            if(bcn .eq. 3) then
              call cost(n,fr,w)
            else
              call rfftb(n,fr,w)
            endif
c
            psi3d(i,1:n,q)=fr(1:n)
         end do
      end do
c
c - - deallocate work array w:
c
      deallocate(w)
c
c - - As a final step, remove 1/r artifact from psi, created by conflict
c - - in boundary conditions psi=0 at R_SS and floating photospheric psi
c - - caused by homogenous Neumann BC in theta for photospheric solution 
c
      mflux=0.d0
      call psi_fix_ss(m,n,p,bcn,a,b,c,d,rsun,rssmrs,psi3d,mflux)
c
c - - we're done:
c
      return
      end
