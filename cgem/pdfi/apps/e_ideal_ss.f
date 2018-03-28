      subroutine e_ideal_ss(m,n,btcoe,bpcoe,brcoe,
     1 et,ep,er,rsun,sinth,a,b,c,d,eti,epi,eri)
c
c+
c    Purpose:  To compute the non-inductive electric field which will make
c              the total electric field approximately perpendicular to B.
c              This subroutine finds eti,epi,eri - the non-inductive 
c              electric field components due to due the gradient of a scalar 
c              potential, which, when added to the input electric field 
c              et,ep,er minimizes the component of E parallel to B.
c              Uses Brian Welsch's relaxation method, implemented in 
c              relax_psi_3d_ss.f .
c
c    Usage: call e_ideal_ss(m,n,btcoe,bpcoe,brcoe,
c    1 et,ep,er,rsun,sinth,a,b,c,d,eti,epi,eri)
c
c    Input: btcoe [m+1,n+1] -  COE Bt from HMI data[Gauss]
c    Input: bpcoe [m+1,n+1] - COE Bp from HMI data [Gauss]
c    Input: brcoe [m+1,n+1] - COE Br from HMI data [Gauss]
c    Input: et [m,n+1] - input electric field multiplied by the speed of light 
c                        in theta (colat) direction [Gauss km/s]
c    Input: ep [m+1,n] - input electric field multiplied by the speed of light 
c                        in phi (azimuth) direction [Gauss km/s]
c    Input: er [m+1,n+1] - input electric field  multiplied by the speed of light 
c                        in radial direction [Gauss km/s]
c
c    Input: sinth, the value of sin(colatitude), computed for theta edge
c    locations.  Is dimensioned m+1.
c 
c    Input: a,b - Minimum and maximum values of co-latitude in radians
c    corresponding to the range of the theta (colatitude) edge values.
c
c    Input: c,d - Minimum and maximum values of longitude in radians 
c    corresponding to the range of longitude (azimuth) edge values.
c
c    Local variable: max_iter, number of iterations for relaxation code
c
c    Output: eti - theta component of the ideal non-inductive electric 
c    field contribution, multiplied by speed of light, dimensioned m,n+1 
c    and defined on phi edges (PE) [Gauss km/s]. 
c 
c    Output: epi - phi component of the ideal non-inductive electric 
c    field contribution, multiplied by speed of light, dimensioned m+1,n
c    and defined on theta edges (TE) [Gauss km/s]. 
c 
c    Output: eri [m+1,n+1] (COE) - radial component of the ideal 
c    non-inductive electric field contribution,multiplied by speed of light, 
c    dimensioned m+1,n+1 and defined in COE locations [Gauss km/s].
c    GHF, 2015
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
c     INPUT DATA
c
      integer :: m,n
      real*8 :: rsun,bthr
      real*8 :: sinth(m+1)
      real*8 :: btcoe(m+1,n+1),bpcoe(m+1,n+1),brcoe(m+1,n+1)
      real*8 :: et(m,n+1),ep(m+1,n),er(m+1,n+1)
      real*8 :: a,b,c,d
c
c - - output variables:
c
      real*8 :: eti(m,n+1),epi(m+1,n),eri(m+1,n+1)
c      
c - - local variables to e_ideal_ss:
c
      integer :: max_iter,verbose
      real*8 :: dtheta,dphi
      real*8 :: btt(m-1,n-1),btp(m-1,n-1),btr(m-1,n-1)
      real*8 :: edt(m-1,n-1),edp(m-1,n-1),edr(m-1,n-1)
      real*8 :: psi(m+1,n+1)
      real*8 :: gradt(m,n+1),gradp(m+1,n),dpsi_dr(m+1,n+1)
      real*8 :: etot(m-1,n-1,3),bvec(m-1,n-1,3),evec(m-1,n-1,3)
c      
c - - set max_iter,verbose:
c
      max_iter=25
c - - set verbose to 0 to eliminate printed diagnostics from relax_psi_3d_ss
      verbose=1
c
      bthr=0.d0
      dtheta=(b-a)/m
      dphi=(d-c)/n  
c     
c - - interpolate et, ep from edges to interior corners as edt,edp:
c
      call interp_eh_ss(m,n,et,ep,edt,edp)
c
c - - get interior corner values of er into edr:
      edr(1:m-1,1:n-1)=er(2:m,2:n)
c 
c - - Define evec array on interior corners
c
      evec(1:m-1,1:n-1,1)=edt
      evec(1:m-1,1:n-1,2)=edp
      evec(1:m-1,1:n-1,3)=edr
c
c    components of B on interior corners only
c
      btt(1:m-1,1:n-1)=btcoe(2:m,2:n)
      btp(1:m-1,1:n-1)=bpcoe(2:m,2:n)
      btr(1:m-1,1:n-1)=brcoe(2:m,2:n)
c
c - - Define bvec array on interior corners
c
      bvec(1:m-1,1:n-1,1)=btt
      bvec(1:m-1,1:n-1,2)=btp
      bvec(1:m-1,1:n-1,3)=btr
c
c - - call relaxation code to get psi, dpsi_dr:
c
c       output (with ghostzones): psi(m+1,n+1),dpsi_dr(m+1,n+1),etot(m-1,n-1,3)
      call relax_psi_3d_ss(m,n,bvec,evec,rsun,a,b,c,d,
     1 max_iter,bthr,verbose,psi,dpsi_dr,etot)
c
c - - Take gradient of psi in horizontal directions:
c
c  input:   psi,[m+1,n+1] sinth(m+1)
c  output: TE, PE gradt[m,n+1], gradp[m+1,n]  
      call gradh_co_ss(m,n,psi,rsun,sinth,dtheta,
     1 dphi,gradt,gradp)
c
c - - set ideal electric field as minus the gradient of psi:
c
      eti(1:m,1:n+1)=-gradt(1:m,1:n+1)
      epi(1:m+1,1:n)=-gradp(1:m+1,1:n)
      eri(1:m+1,1:n+1)=-dpsi_dr(1:m+1,1:n+1)
c
c - - we're done
c
      return
      end
