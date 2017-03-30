      subroutine e_doppler_ss(m,n,btcoe,bpcoe,brcoe,
     1 vlos,lt,lp,lr,rsun,sinth,a,b,c,d,gradt,gradp,dpsi_dr)
c
c+
c    Purpose: To compute non-inductive contribution to the electric field
c             from observed Doppler velocity and transverse magnetic field
c             components.  See sections 2.3.1 and 2.3.3 of 
c             Kazachenko et al. 2014
c    Usage: call e_doppler_ss(m,n,btcoe,bpcoe,brcoe,vlos,lt,lp,lr,rsun,sinth,
c   1 a,b,c,d,gradt,gradp,dpsi_dr)
c
c    This subroutine solves for gradt,gradp,dpsi_dr - the non-inductive 
c    electric field components due to due to Doppler VLOS; See ;
c    Input: btcoe [m+1,n+1] -  COE Bt from HMI data[Gauss]
c    Input: bpcoe [m+1,n+1] - COE Bp from HMI data [Gauss]
c    Input: brcoe [m+1,n+1] - COE Br from HMI data [Gauss]
c    Input: vlos [m+1,n+1] - COE Doppler velocity from HMI data [km/s] 
c    Input: lt [m+1,n+1] - COE co-latitudinal component of the LOS unit vector 
c    for HMI
c    Input: lp [m+1,n+1] - COE azimuthal component of the LOS  unit vector 
c    for HMI
c    Input: lr [m+1,n+1] - COE r-component of the LOS unit vector for HMI
c    Input: rsun - units for the radius of the Sun.
c    Input: sinth, the value of sin(colatitude), computed for theta edge
c    locations.  Is dimensioned m+1.
c 
c    Input: a,b - Minimum and maximum values of co-latitude in radians
c    corresponding to the range of the theta (colatitude) edge values.
c
c    Input: c,d - Minimum and maximum values of longitude in radians 
c    corresponding to the range of longitude (azimuth) edge values.
c
c    Local variable: sigma, value of the width of the PIL in pixels.
c    Local variable: max_iter, number of iterations for relaxation code
c
c    Output: gradt - theta component of the Doppler non-inductive electric 
c    field contribution, dimensioned m,n+1 and defined on phi edges (PE). 
c 
c    Output: gradp - phi component of the Doppler non-inductive electric 
c    field contribution, dimensioned m+1,n and defined on theta edges (TE). 
c 
c    Output: dpsi_dr [m+1,n+1] (COE) - radial component of the Doppler 
c    non-inductive electric field contribution, dimensioned m+1,n+1 and 
c    defined in COE locations.
c    MKD, 2015
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
      real*8 :: vlos(m+1,n+1)
      real*8 :: lt(m+1,n+1),lp(m+1,n+1),lr(m+1,n+1)
      real*8 :: a,b,c,d
c
c - - output variables:
c
      real*8 :: gradt(m,n+1),gradp(m+1,n),dpsi_dr(m+1,n+1)
c      
c - - local variables to e_doppler_ss:
c
      integer :: max_iter,verbose
      real*8 :: dtheta,dphi,sigma
      real*8 :: tt(m-1,n-1),tp(m-1,n-1),tr(m-1,n-1)
      real*8 :: blos(m-1,n-1)
      real*8 :: blost(m-1,n-1),blosp(m-1,n-1),blosr(m-1,n-1)
      real*8 :: btt(m-1,n-1),btp(m-1,n-1),btr(m-1,n-1)
      real*8 :: bdotl(m-1,n-1),btrn(m-1,n-1)
      real*8 :: omega(m-1,n-1)
      real*8 :: edt(m-1,n-1),edp(m-1,n-1),edr(m-1,n-1)
      real*8 :: edmag(m-1,n-1)
      real*8 :: tmp1(m-1,n-1),tmp2(m-1,n-1),tmp3(m-1,n-1)
      real*8 :: psi(m+1,n+1)
      real*8 :: etot(m-1,n-1,3),ed_h(m-1,n-1,3),q(m-1,n-1,3)
c      
c - - set sigma,max_iter,verbose:
c
      sigma=1.0
      max_iter=25
c - - set verbose to 0 to eliminate printed diagnostics from relax_psi_3d_ss
      verbose=1
c
      bthr=0.d0
      dtheta=(b-a)/m
      dphi=(d-c)/n  
c
      tt=vlos(2:m,2:n)*lt(2:m,2:n)
      tp=vlos(2:m,2:n)*lp(2:m,2:n)
      tr=vlos(2:m,2:n)*lr(2:m,2:n)
c     
c    e=[edt,edp,edr]=-vlos \times B. (interior corners, ie CO)   [m-1,n-1]
c     
      edt=-(tp*brcoe(2:m,2:n)-tr*bpcoe(2:m,2:n))
      edp=-(tr*btcoe(2:m,2:n)-tt*brcoe(2:m,2:n))
      edr=-(tt*bpcoe(2:m,2:n)-tp*btcoe(2:m,2:n))
c
c    components of B computed on interior corners only
c
      blos=sqrt(btcoe(2:m,2:n)*lt(2:m,2:n)*btcoe(2:m,2:n)*lt(2:m,2:n)+
     1 bpcoe(2:m,2:n)*lp(2:m,2:n)*bpcoe(2:m,2:n)*lp(2:m,2:n)+
     2 brcoe(2:m,2:n)*lr(2:m,2:n)*brcoe(2:m,2:n)*lr(2:m,2:n))
      bdotl=btcoe(2:m,2:n)*lt(2:m,2:n)+bpcoe(2:m,2:n)*lp(2:m,2:n)
     1 +brcoe(2:m,2:n)*lr(2:m,2:n)
      blost=bdotl*lt(2:m,2:n)
      blosp=bdotl*lp(2:m,2:n)
      blosr=bdotl*lr(2:m,2:n)
      btt=btcoe(2:m,2:n)-blost(1:m-1,1:n-1)
      btp=bpcoe(2:m,2:n)-blosp(1:m-1,1:n-1)
      btr=brcoe(2:m,2:n)-blosr(1:m-1,1:n-1)
      btrn=sqrt(btt*btt+btp*btp+btr*btr)
c 
      where(btrn .ne. 0.d0) 
          omega=exp(-blos*blos/(btrn*btrn*sigma*sigma))
      elsewhere
          omega=0.d0
      endwhere
c
c - - ed_h on interior corners only
c
      ed_h(1:m-1,1:n-1,1)=edt*omega
      ed_h(1:m-1,1:n-1,2)=edp*omega
      ed_h(1:m-1,1:n-1,3)=edr*omega
c
      edmag=sqrt(ed_h(1:m-1,1:n-1,1)*ed_h(1:m-1,1:n-1,1)+
     1 ed_h(1:m-1,1:n-1,2)*ed_h(1:m-1,1:n-1,2)+
     2 ed_h(1:m-1,1:n-1,3)*ed_h(1:m-1,1:n-1,3))
c    
      where(edmag .ne. 0.d0) 
        tmp1=edt*omega/edmag
        tmp2=edp*omega/edmag
        tmp3=edr*omega/edmag
      elsewhere
        tmp1=0.d0
        tmp2=0.d0
        tmp3=0.d0
      endwhere
c
c     q on interior corners
c
      q(1:m-1,1:n-1,1)=tmp1
      q(1:m-1,1:n-1,2)=tmp2
      q(1:m-1,1:n-1,3)=tmp3

c       input: q[m-1,n-1,3],ed_h(m-1,n-1,3)
c       output (with ghostzones): psi(m+1,n+1),dpsi_dr(m+1,n+1),etot(m-1,n-1,3)
      call relax_psi_3d_ss(m,n,q,ed_h,rsun,a,b,c,d,
     1 max_iter,bthr,verbose,psi,dpsi_dr,etot)
     
c  input:   psi,[m+1,n+1] sinth(m+1)
c  output: TE, PE gradt[m,n+1], gradp[m+1,n]  
      call gradh_co_ss(m,n,psi,rsun,sinth,dtheta,
     1 dphi,gradt,gradp)
c
c - - we're done; components of gradient of scalar function are
c - - returned as gradt, gradp      
c
      return
      end
