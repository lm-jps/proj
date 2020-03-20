      subroutine ahpot_ss(m,n,p,bcn,a,b,c,d,rsun,rssmrs,scrb3d,mflux,
     1           atpot,appot)
c
c+ - - Purpose: compute potential field value of at,ap (vector potential
c               components) given the 3-d array scrb3d by taking curl of 
c               scriptB rhat. scrb3d is computed by subroutine scrbpot_ss.
c               For nonzero net flux, an additional term as added to appot
c               to account for the net flux contribution.
c
c - -  Usage:   call ahpot_ss(m,n,p,bcn,a,b,c,d,rsun,rssmrs,scrb3d,mflux,
c               atpot,appot)
c
c - -  Input:   m,n,p: integer values of numbers of cell centers in theta,
c               phi, and r directions
c
c - -  Input:   bcn: integer value of boundary condition type in phi direction.
c               bcn=0 => periodic BC; bcn=3 => homogenous Neumann BC
c
c - -  Input:   a,b,c,d:  real*8 values of min, max colatitude, min, max
c               values of longitude. [radians]
c
c - -  Input:   rsun,rssmrs: real*8 values of the radius of sun, and 
c               distance from phot to source surface. [km] Normally rsun=6.96d5.
c
c - -  Input:   scrb3d(m,n,p+1): real*8 array of poloidal potential scribtb
c               [G km^2]
c
c - -  Input:   mflux: real*8 value of the net radial magnetic flux [G km^2].
c               If you wanto assume zero net flux, set mflux=0.d0 on input.
c               mflux can be calculated with subroutine mflux_ss.
c
c - -  Output:  atpot(m,n+1,p+1): real*8 array of theta-comp of vector potential
c               [G km]
c
c - -  Output:  appot(m+1,n,p+1): real*8 array of phi-comp of vector potential
c               [G km]
c
c - -  Note:    atpot,appot are computed on phi and theta edges, and on
c               radial shells.
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
      real*8 :: a,b,c,d,rsun,rssmrs,mflux
      real*8 :: scrb3d(m,n,p+1)
c 
c - - output variables:
c
      real*8 :: atpot(m,n+1,p+1)
      real*8 :: appot(m+1,n,p+1)
c
c - - local variables:
c
      integer :: q,i,j
      real*8 :: scrb(m+2,n+2),bdas(n),bdbs(n),bdcs(m),bdds(m)
      real*8 :: sinth(m+1),sinth_hlf(m),r(p+1)
      real*8 :: costh(m+1),costh_hlf(m),theta(m+1),theta_hlf(m)
      real*8 :: curlt(m,n+1),curlp(m+1,n)
      real*8 :: dphi,dtheta,dtheta_hlf,delr,area,b0
c
c - - check illegal values of bcn:
c
      if((bcn .ne. 0) .and. (bcn .ne. 3)) then
         write(6,*) 'ahpot_ss: Illegal bcn = ',bcn,' exiting'
         stop
      endif
c
      bdas(1:n)=0.d0
      bdbs(1:n)=0.d0
      bdcs(1:m)=0.d0
      bdds(1:m)=0.d0
c
      delr=rssmrs/p
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - compute cos(theta) arrays
c
      dtheta_hlf=0.5d0*dtheta
      do i=1,m
         theta(i)=a + (i-1)*dtheta
         theta_hlf(i)=a+dtheta_hlf+(i-1)*dtheta
      enddo
      theta(m+1)=b
c
      costh(1:m+1)=cos(theta(1:m+1))
      costh_hlf(1:m)=cos(theta_hlf(1:m))
c
c - - calculate area of the spherical wedge at photosphere:
c
      area=(d-c)*(cos(a)-cos(b))*rsun**2
c
c - - get value of b0 by dividing mflux by area:
c
      b0=mflux/area
c
c
c - - define radial shell arrays:
c
      r(1)=rsun
      do q=2,p+1
         r(q)=r(1)+(q-1)*delr
      end do
c
      do q=1,p+1
c - - get script B on radial shells, set equal to scrb array:
         scrb(2:m+1,2:n+1)=scrb3d(1:m,1:n,q)
c - - fill in ghost zones in the shell surface:
         scrb(1,2:n+1)=scrb(2,2:n+1)-1.d0*dtheta*bdas(1:n)
         scrb(m+2,2:n+1)=scrb(m+1,2:n+1)+1.d0*dtheta*bdbs(1:n)
c
c - -    Neumann BC in phi:
         if(bcn .eq. 3) then
            if(q .eq. 1) then
              scrb(2:m+1,1)=scrb(2:m+1,2)-1.d0*dphi*bdcs(1:m)
              scrb(2:m+1,n+2)=scrb(2:m+1,n+1)+1.d0*dphi*bdds(1:m)
            else
              scrb(2:m+1,1)=scrb(2:m+1,3)-2.d0*dphi*bdcs(1:m)
              scrb(2:m+1,n+2)=scrb(2:m+1,n)+2.d0*dphi*bdds(1:m)
            endif
         else
c - -       Periodic BC in phi:
            scrb(2:m+1,1)=scrb(2:m+1,n+1)
            scrb(2:m+1,n+2)=scrb(2:m+1,2)
         endif
c
         call curl_psi_rhat_ce_ss(m,n,scrb,r(q),sinth_hlf,dtheta,
     1        dphi,curlt,curlp)
         atpot(1:m,1:n+1,q)=curlt(1:m,1:n+1)
         appot(1:m+1,1:n,q)=curlp(1:m+1,1:n)
c
c - - put in term in A_phi to account for non-zero net flux:
c
         do j=1,n
            appot(1:m+1,j,q)=appot(1:m+1,j,q)-b0*((rsun**2)/r(q))
     1                      *costh(1:m+1)/sinth(1:m+1)
         enddo
c
      enddo
c
      return
      end
