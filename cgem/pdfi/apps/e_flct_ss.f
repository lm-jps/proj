      subroutine e_flct_ss(m,n,brte,brpe,bhte,bhpe,
     1 vt,vp,rsun,sinth,sinth_hlf,a,b,c,d,ezetat,ezetap)
c
c - - documentation below can be seen using DOC_LIB,'e_flct_ss.f'
c+
c - - Purpose:  To compute non-inductive contribution to the electric field
c               from horizontal flows determined by FLCT (correlation tracking)
c               and radial magnetic fields.  See section 2.3.2 of Kazachenko
c               et al (2014).  Solves for ezetat, ezetap - 
c               non-inductive electric field components due to horizontal 
c               FLCT velocity field [vt,vp] 
c - - Usage:  
c - - call e_flct_ss(m,n,brte,brpe,bhte,bhpe,vt,vp,rsun,sinth,sinth_hlf,
c    1 a,b,c,d,ezetat,ezetap)
c
c - - Input: brte - array of the radial component of the magnetic field,
c - - dimensioned m+1,n and defined on theta edges, midway between
c - - phi edges.
c
c - - Input: brpe - array of the radial component of the magnetic field,
c - - dimensioned m,n+1 and defined on phi edges, midway between
c - - theta edges.
c
c - - Input: bhte - array of the horizontal amplitude of the magnetic field,
c - - dimensioned m+1,n and defined on theta edges
c
c - - Input: bhpe - array of the horizontal amplitude of the magnetic field,
c - - dimensioned m,n+1 and defined on phi edges
c
c - - Input:  vt - the component of the velocity field pointing in the
c - - colatitude direction (positive is southward pointing).  
c - - vt is dimensioned m+1,n and is defined on theta edges, midway between
c - - phi edges.
c
c - - Input:  vp - the component of the velocity field pointing in the
c - - logitude direction, is dimensioned m,n+1 and is defined on phi edges, 
c - - midway between theta edges.
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
c - - Input: sigma, value of the width of the PIL in pixels.
c
c - - Output: ezetat - theta component of the FLCT non-inductive electric 
c - - field contribution, dimensioned m,n+1 and defined on theta edges PE. 
c 
c - - Output: ezetap - phi component of the FLCT non-inductive electric 
c - - field contribution, dimensioned m+1,n and defined on theta edges TE. 
c-
c - - MKD, October 2015
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
      integer :: m,n
      integer :: m1,n1
      real*8 :: rsun
      real*8 :: sinth(m+1),sinth_hlf(m)
      real*8 :: vt(m+1,n),vp(m,n+1)
      real*8 :: brte(m+1,n),brpe(m,n+1)
      real*8 :: bhte(m+1,n),bhpe(m,n+1)
      real*8 :: a,b,c,d,sigma
      real*8 :: ezetat(m,n+1),ezetap(m+1,n)
c
c - - local allocatable work array:
c
      real*8, allocatable :: w(:)
c      
c - - local variables to e_flct_ss:
c
      real*8 :: dtheta,dphi,elm,pertrb
      integer :: idimf,itmp,bcm,bcn,ierror
      real*8 :: eflctpe(m,n+1),eflctte(m+1,n),rhspe(m,n+1),rhste(m+1,n)
      real*8 :: omega_te(m+1,n),omega_pe(m,n+1)
      real*8 :: bdas(n-1),bdbs(n-1),bdcs(m-1),bdds(m-1)
      real*8 :: ap,bpp,cp,dp
      real*8 :: divh_rhs(m-1,n-1),fce(m-1,n-1),zeta(m+1,n+1)
      real*8 :: gradt(m,n+1),gradp(m+1,n)

      sigma=1.0
      dtheta=(b-a)/m
      dphi=(d-c)/n  
      m1=m-1
      n1=n-1
  
c
c - - set boundary condition flags for Neumann:
c    
      bcm=3
      bcn=3
      elm=0.d0

c
c vh \cross br = [et,ep]=[vp br, -vt br]
c cE=-vh \cross br=[-vp br, vt br]
c We need to calculate omega at both PE and TE edges
c
c TE:  
c
      eflctte=brte*vt
c
      where(abs(bhte) .ne. 0.d0) 
          omega_te=exp(-(brte*brte)/(bhte*bhte*sigma*sigma))
      elsewhere
          omega_te=0.d0
      endwhere
c          
      rhste=eflctte-omega_te*eflctte
c PE:
      eflctpe=-brpe*vp
      where(abs(bhpe) .ne. 0.d0) 
          omega_pe=exp(-(brpe*brpe)/(bhpe*bhpe*sigma*sigma))
      elsewhere
          omega_pe=0.d0
      endwhere
c      
      rhspe=eflctpe-omega_pe*eflctpe
      call divh_co_ss(m,n,rhspe,rhste,rsun,sinth,sinth_hlf,dtheta,dphi,
     1 divh_rhs)
c
      ap=a+0.5d0*dtheta
      bpp=b-0.5d0*dtheta
      cp=c+0.5d0*dphi
      dp=d-0.5d0*dphi
c      
c - - compute size of work array, and allocate it
c  M*INT(LOG2(N))
      itmp=(4*(n+1)+int(log(real(n+1))/log(2.0)+16)*(m+1))
      allocate(w(itmp))
      idimf=m1
c
c
c - - Set RHS for Poisson equation for fce:
c
      fce(1:m-1,1:n-1)=-divh_rhs(1:m-1,1:n-1)*rsun**2
      
c - - set homogenous Neumann boundary conditions for fce:
       bdas(1:n-1)=0.d0
       bdbs(1:n-1)=0.d0
       bdcs(1:m-1)=0.d0
       bdds(1:m-1)=0.d0
c - - bdas,bdbs[n-1], bdcs,bdds[m-1],fce[m-1,n-1]  
      call hstssp(ap,bpp,m1,bcm,bdas,bdbs,cp,dp,n1,bcn,bdcs,bdds,elm,
     1 fce,idimf,pertrb,ierror,w)
       if(ierror .ne. 0) then
         write(6,*) 'e_flct_ss: zeta ierror .ne. 0; = ',ierror
         stop
      endif 
c
c - - put solution points into ezeta array psi=zeta[m+1,n+1]
c
      zeta(2:m,2:n)=fce(1:m-1,1:n-1)
c
c - - Now fill in ghost zone values using boundary condition arrays:
c  instead of n+1 -> n; m+1 -> m; n->n-1: m->m-1
      zeta(1,2:n)=zeta(2,2:n)-1.d0*dtheta*bdas(1:n-1)
      zeta(m+1,2:n)=zeta(m,2:n)+1.d0*dtheta*bdbs(1:n-1)
      zeta(2:m,1)=zeta(2:m,2)-1.d0*dphi*bdcs(1:m-1)
      zeta(2:m,n+1)=zeta(2:m,n)+1.d0*dphi*bdds(1:m-1)
c
c - - fill in ghost zone corners (not necessary, but cosmetic)
c
      zeta(1,1)=0.5d0*(zeta(1,2)+zeta(2,1))
      zeta(m+1,1)=0.5d0*(zeta(m+1,2)+zeta(m,1))
      zeta(1,n+1)=0.5d0*(zeta(2,n+1)+zeta(1,n))
      zeta(m+1,n+1)=0.5d0*(zeta(m+1,n)+zeta(m,n+1))
c    
c - - ezeta: gradt[m,np1], gradp[mp1,n]
c
      call gradh_co_ss(m,n,zeta,rsun,sinth,dtheta,dphi,gradt,
     1 gradp)
c
      ezetat(1:m,1:n+1)=-gradt(1:m,1:n+1)
      ezetap(1:m+1,1:n)=-gradp(1:m+1,1:n)
c
c - - deallocate w to avoid memory leak
c        
      deallocate(w)
c - - we're done
      return
      end
