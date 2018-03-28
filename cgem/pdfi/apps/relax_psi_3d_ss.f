      subroutine relax_psi_3d_ss(m,n,bvec,evec,rsun,a,b,c,d,
     1 max_iter,bthr,verbose,psi,dpsi_dr,etot)
c
c+
c - - Purpose: Relaxation code to compute an electric field contribution
c - - from the gradient of a scalar
c - - potential, which then results in the total electric field being (close to)
c - - perpendicular to the magnetic field.  Basic relaxation algorithm 
c - - was developed by Brian Welsch, and is 
c - - described in section 3.2 of Fisher et al. 2010, ApJ 715, 242, with 
c - - subsequent changes described in sections 2.2 and 4.1 of 
c - - Kazachenko et al., 2014, ApJ 795, 17.
c
c - - Usage:  call relax_psi_3d_ss(m,n,bvec,evec,rsun,a,b,c,d,
c    1 max_iter,bthr,verbose,psi,dpsi_dr,etot)
c
c - - Input:  m,n - number of cell centers in theta, phi directions, resp.
c - - Input:  bvec(m-1,n-1,3) - B_theta, B_phi, and B_r evaluated on interior
c             cell corners.
c - - Input:  evec(m-1,n-1,3) - input electric field, multiplied by the speed
c             of light, (E_theta, E_phi, E_r), evaluated on interior cell
c             corners [G km/s]
c - - Input:  rsun - radius of the Sun
c - - Input:  a,b,c,d - colatitude (a,b) and longitude (c,d) range of the
c             exterior boundary of the problem (exterior cell corners)
c - - Input:  max_iter - Maximum number of iterations for "iterative method"
c             of Welsch, (see section 3.2 of Fisher et al. [2010])
c - - Input:  bthr - magnetic field threshold for parallelizing E field
c - - Input:  verbose - integer set to nonzero value to print diagnostics
c - - Output: psi(m+1,n+1) - scalar potential values on COE grid corners
c - - Output: dpsi_dr(m+1,n+1) - radial derivative of psi on COE grid 
c             corners
c - - Output: etot(m-1,n-1,3) - (c*evec - grad psi) on interior corners,
c             i.e. etot is in fact electric field, multiplied by the 
c             speed of light [G km/s]
c - - Note:   Unlike the rest of PDFI_SS, all vector calculus operations
c             within this subroutine use centered grid formalism.  Only
c             the last step of setting boundary conditions for psi assumes
c             staggered grid.
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
c - - calling argument declarations:
c
      integer :: m,n,max_iter,verbose
c
      real*8 :: bvec(m-1,n-1,3),evec(m-1,n-1,3)
      real*8 :: rsun,a,b,c,d,bthr
      real*8 :: psi(m+1,n+1),dpsi_dr(m+1,n+1),etot(m-1,n-1,3)
c
c - - local variable declarations
c
      integer :: itmax,mc,nc,bcm,bcn,idimf,
     1          mcp1,ncp1,ierror,itmp,keep_going,iter
      real*8 :: ap,bp,cp,dp,dtheta,dphi,bmthr,elm,pertrb
      real*8 :: emin,emax,emean,eptot
      real*8 :: s1(m-1,n-1),s2(m-1,n-1),s3(m-1,n-1)
      real*8 :: bth(m-1,n-1),bph(m-1,n-1),br(m-1,n-1)
      real*8 :: bhatth(m-1,n-1),bhatph(m-1,n-1),bhatr(m-1,n-1)
      real*8 :: bmag(m-1,n-1),psi_source0(m-1,n-1),psi0(m+1,n+1)
      real*8 :: bhmag(m-1,n-1),dpsidth(m-1,n-1),dpsidph(m-1,n-1)
      real*8 :: dpsidth0(m-1,n-1),dpsidph0(m-1,n-1)
      real*8 :: dpsidrng(m-1,n-1),dpsidr0(m-1,n-1),epara(m-1,n-1)
      real*8 :: eparapsi(m-1,n-1),eparapsi0(m-1,n-1)
      real*8 :: sinth(m-1),sinth_gh(m+1),sinth_hlf(m-2),
     1          sinth_hlf_gh(m),psi_source(m-1,n-1),f(m-1,n-1)
      real*8 :: eth(m-1,n-1),eph(m-1,n-1),er(m-1,n-1)
      real*8 :: bca(n-1),bcb(n-1),bcc(m-1),bcd(m-1)
      real*8, allocatable :: w(:)
c
c - - some sdf variables (uncomment if needed):
c
c     integer*8 :: dims(20)
c     integer :: nbpw,ndim
c     character*100 :: lab
c
      itmax=max_iter
c - - don't change max_iter, just use local value and change it if necessary
      if(itmax .lt. 0) itmax=100
c - - don't change bthr, just use local value and change it if necessary
      bmthr=bthr
      if(bmthr .lt. 0.d0) bmthr=0.d0
c
c - - Set up boundary limits for interior corner problem
c
      dtheta=(b-a)/m
      dphi=(d-c)/n
      ap=a+dtheta
      bp=b-dtheta
      cp=c+dphi
      dp=d-dphi
c
c - - centered grid, interior values of m,n (for use in sinthta_sc,
c - - delh2_sc, divh_sc.  Note these are two smaller than m,n of full grid.
c
      mc=m-2
      nc=n-2
      mcp1=mc+1
      ncp1=nc+1
c
c - - Get sin(theta), etc for our grid using centered grid formalism
c
      call sinthta_sc(ap,bp,mc,sinth,sinth_hlf,sinth_gh,
     1 sinth_hlf_gh)
c - - components of B (B_theta, B_phi, B_r):
      bth(1:m-1,1:n-1)=bvec(1:m-1,1:n-1,1)
      bph(1:m-1,1:n-1)=bvec(1:m-1,1:n-1,2)
      br(1:m-1,1:n-1)=bvec(1:m-1,1:n-1,3)
c
      bmag=sqrt(bth*bth+bph*bph+br*br)
c - - Now set components of unit vector bhat, but set to 0 where B < bthr
      where(bmag .gt. bmthr) 
          bhatth=bth/bmag
          bhatph=bph/bmag
          bhatr=br/bmag
          bhmag=sqrt(bhatth**2+bhatph**2)
      elsewhere
          bhatth=0.d0
          bhatph=0.d0
          bhatr=0.d0
          bhmag=0.d0
      endwhere
c - - components of E (E_theta, E_phi, E_r):
      eth(1:m-1,1:n-1)=evec(1:m-1,1:n-1,1)
      eph(1:m-1,1:n-1)=evec(1:m-1,1:n-1,2)
      er(1:m-1,1:n-1)=evec(1:m-1,1:n-1,3)
c - - compute component of E parallel to B, set s1:
      epara=eth*bhatth + eph*bhatph + er*bhatr
      s1(1:m-1,1:n-1)=epara
      s2(1:m-1,1:n-1)=0.d0
      s3(1:m-1,1:n-1)=0.d0
c - - initial value of dpsidr:
      dpsidr0=epara*bhatr
c
c - - Initial guess:
c
c - - compute horizontal divergence of s1*bhat:
c     (Note call to centered grid subroutines uses mc,nc not m,n!)
      call divh_sc(mc,nc,s1*bhatth,s1*bhatph,rsun,sinth,dtheta,dphi,
     1 psi_source0)
c
c - - RHS to poisson eqn is psi_source0
c
      f(1:m-1,1:n-1)=psi_source0(1:m-1,1:n-1)*(rsun**2)
c
c - - compute size of work array and allocate it (size overestimate)
c
      itmp=(4*(n+2)+int(log(real(n+2))/log(2.0)+16)*(m+2))
      allocate(w(itmp))
      idimf=m-1
c
c - - Set homogenous Neumann boundary conditions: boundary condition flags
      bcm=3
      bcn=3
c - - theta edges:
      bca(1:ncp1)=0.d0
      bcb(1:ncp1)=0.d0
c - - phi edges:
      bcc(1:mcp1)=0.d0
      bcd(1:mcp1)=0.d0
c - - set lambda (Helmholtz eqn term) to 0.
      elm=0.d0
c - - Fishpack call
c - - centered grid version of fishpack call, with mcp1=m-1,ncp1=n-1
      call hwsssp(ap,bp,mc,bcm,bca,bcb,cp,dp,nc,bcn,bcc,bcd,elm,
     1     f,idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'relax_psi_3d_ss init ierror .ne. 0; =',ierror
         stop
      endif
c
c - - use mcp1,ncp1 notation here to make ghost zone indexing easier
c - - to compare with other code using the centered grid format
c
c - - get solution from f into the psi0 array:
      psi0(2:mcp1+1,2:ncp1+1)=f(1:mcp1,1:ncp1)
c
c - - fill in ghost zones at edges, using boundary conditions:
c
      psi0(1,2:ncp1+1)=psi0(3,2:ncp1+1)-2.d0*dtheta*bca(1:ncp1)
      psi0(mcp1+2,2:ncp1+1)=psi0(mcp1,2:ncp1+1)+2.d0*dtheta*
     1    bcb(1:ncp1)
      psi0(2:mcp1+1,1)=psi0(2:mcp1+1,3)-2.d0*dphi*bcc(1:mcp1)
      psi0(2:mcp1+1,ncp1+2)=psi0(2:mcp1+1,ncp1)+2.d0*dphi*
     1    bcd(1:mcp1)
c
c - - fill in ghost zone corners for display reasons
c
      psi0(1,1)=0.5d0*(psi0(1,2)+psi0(2,1))
      psi0(mcp1+2,1)=0.5d0*(psi0(mcp1+1,1)+psi0(mcp1+2,2))
      psi0(1,ncp1+2)=0.5d0*(psi0(1,ncp1+1)+psi0(2,ncp1+2))
      psi0(mcp1+2,ncp1+2)=0.5d0*(psi0(mcp1+2,ncp1+1)+
     1    psi0(mcp1+1,ncp1+2))
c
c - - initialize source for subsequent iteration:
c
      psi(1:mcp1+2,1:ncp1+2)=psi0(1:mcp1+2,1:ncp1+2)
c
c - - compute horizontal gradient of psi:
c
      call gradh_sc(mc,nc,psi,rsun,sinth_gh,dtheta,dphi,dpsidth0,
     1 dpsidph0)
c
c - - keep copies of initial psi derivatives:
c
      dpsidth(1:mcp1,1:ncp1)=dpsidth0(1:mcp1,1:ncp1)
      dpsidph(1:mcp1,1:ncp1)=dpsidph0(1:mcp1,1:ncp1)
      dpsidrng=dpsidr0
      eparapsi0=-(dpsidth*bhatth+dpsidph*bhatph+
     1 dpsidr0*bhatr)
      eparapsi=eparapsi0
c
c - - initial setup complete, now get ready to loop iteratively:
c
      if(itmax .gt. 0) then
         keep_going=1
      else
         keep_going=0
      endif
      iter=0
c
      do while (keep_going .eq. 1) 
c
c - - BEGINNING of iterative loop
c
      iter=iter+1
c - - possible convergence metrics or diagnostics:
      emin=minval(epara+eparapsi)
      emax=maxval(epara+eparapsi)
      emean=sum(abs(epara+eparapsi))/size(epara+eparapsi)
      eptot=sum(abs(epara+eparapsi))
      if(verbose .ne. 0) then 
        if((((iter/5)*5) .eq. iter) .or. (iter .eq. 1)) then
          write(6,*) 'iter,emean,eptot = ',
     1    iter,emean,eptot
        endif
      endif
c
c - - formulate a convergence test (NOT YET IMPLEMENTED)!
c
      continue
c
c - - compute s2,s3 where bhmag .ne. 0.):
c
      where (bhmag .ne. 0.)
         s2=(bhatth*dpsidph - bhatph*dpsidth)/bhmag**2
         s3=dpsidrng - (dpsidth*bhatth + dpsidph*bhatph)*
     1      bhatr/(bhmag**2)
      elsewhere
         s2=0.d0
         s3=0.d0
      endwhere
c
c - - take horizontal divergence:
c
      call divh_sc(mc,nc,s1*bhatth-s2*bhatph-s3*bhatr*bhatth,
     1 s1*bhatph+s2*bhatth-s3*bhatr*bhatph,rsun,sinth,dtheta,dphi,
     2 psi_source)
c
c - - Set rhs to psi_source:
c
      f(1:m-1,1:n-1)=psi_source(1:m-1,1:n-1)*(rsun**2)
c
c - - solve Poisson equation for psi in loop and get ghost zone values
c
c - - centered grid version of fishpack call, with mp1=m-1,np1=n-1
      call hwsssp(ap,bp,mc,bcm,bca,bcb,cp,dp,nc,bcn,bcc,bcd,elm,
     1     f,idimf,pertrb,ierror,w)
      if(ierror .ne. 0) then
         write(6,*) 'relax_psi_3d_ss loop ierror .ne. 0; =',ierror
         stop
      endif
c
c - - use mcp1,ncp1 notation here to make ghost zone indexing easier
c - - to compare with other code using the centered grid format
c
c - - get solution from f into the psi array:
      psi(2:mcp1+1,2:ncp1+1)=f(1:mcp1,1:ncp1)
c
c - - fill in ghost zones at edges, using boundary conditions:
c
      psi(1,2:ncp1+1)=psi(3,2:ncp1+1)-2.d0*dtheta*bca(1:ncp1)
      psi(mcp1+2,2:ncp1+1)=psi(mcp1,2:ncp1+1)+2.d0*dtheta*
     1    bcb(1:ncp1)
      psi(2:mcp1+1,1)=psi(2:mcp1+1,3)-2.d0*dphi*bcc(1:mcp1)
      psi(2:mcp1+1,ncp1+2)=psi(2:mcp1+1,ncp1)+2.d0*dphi*
     1    bcd(1:mcp1)
c
c - - fill in ghost zone corners for display reasons
c
      psi(1,1)=0.5d0*(psi(1,2)+psi(2,1))
      psi(mcp1+2,1)=0.5d0*(psi(mcp1+1,1)+psi(mcp1+2,2))
      psi(1,ncp1+2)=0.5d0*(psi(1,ncp1+1)+psi(2,ncp1+2))
      psi(mcp1+2,ncp1+2)=0.5d0*(psi(mcp1+2,ncp1+1)+
     1    psi(mcp1+1,ncp1+2))
c
c - - compute horizontal gradient of psi:
c
      call gradh_sc(mc,nc,psi,rsun,sinth_gh,dtheta,dphi,dpsidth,
     1 dpsidph)
c
      dpsidrng=s1*bhatr+s3*bhmag**2
c
      eparapsi=-(dpsidth*bhatth+dpsidph*bhatph+
     1 dpsidrng*bhatr)
c
      if (iter .ge. itmax) keep_going=0
c
c - - END of iterative loop
c
      enddo
c
c - - fill in psi boundary values using boundary conditions and
c - - assuming staggered grid before exiting
c
      psi(1,2:ncp1+1)=psi(2,2:ncp1+1)-1.d0*dtheta*bca(1:ncp1)
      psi(mcp1+2,2:ncp1+1)=psi(mcp1+1,2:ncp1+1)+1.d0*dtheta*
     1    bcb(1:ncp1)
      psi(2:mcp1+1,1)=psi(2:mcp1+1,2)-1.d0*dphi*bcc(1:mcp1)
      psi(2:mcp1+1,ncp1+2)=psi(2:mcp1+1,ncp1+1)+1.d0*dphi*
     1    bcd(1:mcp1)
c
c - - fill in psi ghost zone corners for display reasons
c
      psi(1,1)=0.5d0*(psi(1,2)+psi(2,1))
      psi(mcp1+2,1)=0.5d0*(psi(mcp1+1,1)+psi(mcp1+2,2))
      psi(1,ncp1+2)=0.5d0*(psi(1,ncp1+1)+psi(2,ncp1+2))
      psi(mcp1+2,ncp1+2)=0.5d0*(psi(mcp1+2,ncp1+1)+
     1    psi(mcp1+1,ncp1+2))
c
c - - set dpsi_dr interior values to dpsidrng, and set boundary
c - - values to 0
c
      dpsi_dr(1:mcp1+2,1:ncp1+2)=0.d0
      dpsi_dr(2:mcp1+1,2:ncp1+1)=dpsidrng(1:mcp1,1:ncp1)
c
c - - Finally, compute etot:
c
      etot(1:m-1,1:n-1,1)=evec(1:m-1,1:n-1,1)-dpsidth(1:m-1,1:n-1)
      etot(1:m-1,1:n-1,2)=evec(1:m-1,1:n-1,2)-dpsidph(1:m-1,1:n-1)
      etot(1:m-1,1:n-1,3)=evec(1:m-1,1:n-1,3)-dpsidrng(1:m-1,1:n-1)
c
c - - we're done
      deallocate(w)
      return
      end
