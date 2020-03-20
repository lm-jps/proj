      subroutine e_ideal_ss(m,n,btcoe,bpcoe,brcoe,
     1 et,ep,er,rsun,sinth,a,b,c,d,eti,epi,eri)
c
c+
c - -  Purpose:  To compute the non-inductive electric field which will make
c                the total electric field approximately perpendicular to B.
c                This subroutine finds eti,epi,eri - the non-inductive 
c                electric field components due to due the gradient of a scalar 
c                potential, which, when added to the input electric field 
c                et,ep,er minimizes the component of E parallel to B.
c                Uses Brian Welsch's relaxation method, implemented in 
c                relax_psi_3d_ss.f .
c
c - -  Usage:    call e_ideal_ss(m,n,btcoe,bpcoe,brcoe,
c                et,ep,er,rsun,sinth,a,b,c,d,eti,epi,eri)
c
c - -  Input:    btcoe(m+1,n+1) - real*8 array of theta component of B on
c                COE grid [G]
c
c - -  Input:    bpcoe(m+1,n+1) - real*8 array of phi component of B on
c                COE grid [G]
c
c - -  Input:    brcoe(m+1,n+1) - real*8 array of radial component of B on
c                COE grid [G]
c
c - -  Input:    et(m,n+1) - real*8 array of input electric field in theta 
c                (colat) direction multiplied by c, on PE grid [G km/s]
c
c - -  Input:    ep(m+1,n) - real*8 array of input electric field in phi 
c                (lon) direction multiplied by c, on TE grid [G km/s]
c
c - -  Input:    er(m+1,n+1) - real*8 array of input electric field in 
c                radial direction multiplied by c, on COE grid [G km/s]
c
c - -  Input:    sinth(m+1) -  real*8 array of the sin(colatitude), 
c                computed on theta edges.
c
c - -  Input:    a,b - real*8 values of minimum and maximum values of 
c                co-latitude corresponding to the range of the 
c                theta (colatitude) edge values. [radians]
c
c - -  Input:    c,d - real*8 values of minimum and maximum values of 
c                longitude corresponding to the range of longitude (azimuth) 
c                edge values. [radians]
c
c - - Output:    eti(m,n+1) - real*8 array of theta component of the ideal, 
c                non-inductive electric field contribution, multiplied by c,
c                and defined on phi edges (PE grid) [G km/s]. 
c
c - - Output:    epi(m+1,n) - real*8 array of phi component of the ideal, 
c                non-inductive electric field contribution, multiplied by c,
c                and defined on theta edges (TE grid) [G km/s].
c
c - - Output:    eri(m+1,n+1) - real*8 array of radial component of the ideal 
c                non-inductive electric field contribution, multiplied by 
c                c, and defined at COE grid locations [G km/s].
c
c    GHF, 2015
c
c    Local variable: max_iter, number of iterations for relaxation code
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
c     INPUT DATA
c
      integer :: m,n
      real*8 :: rsun
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
      real*8 :: dtheta,dphi,bthr
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
      verbose=0
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
c
      edr(1:m-1,1:n-1)=er(2:m,2:n)
c 
c - - Define evec array on interior corners
c
      evec(1:m-1,1:n-1,1)=edt
      evec(1:m-1,1:n-1,2)=edp
      evec(1:m-1,1:n-1,3)=edr
c
c - - components of B on interior corners only
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
