      subroutine curle3dphot_ss(m,n,a,b,c,d,rsun,et,ep,
     1           er,detdr,depdr,btt,bpt,brt)
c
c+
c - - Purpose:  Given cE and detdr and depdr at the photosphere,
c               compute the curl of cE, and express the results as btt,
c               bpt,brt, the negative of the curl of E in the theta,phi,
c               and radial directions.  The results can be viewed as the
c               time derivative of the theta,phi, and radial components of B
c               through Faraday's law.
c
c - - Usage:    call curle3dphot_ss(m,n,a,b,c,d,rsun,et,ep,er,
c               detdr,depdr,btt,bpt,brt)
c
c - - Input:    m,n - integer values of the number of cell centers in the
c               colatitude (latitude) direction, and longitude direction, resp.
c
c - - Input:    a,b,c,d - real*8 values of the limits in colatitude 
c               (a and b; 0 <= a < b <= pi) and longitude (c <=0 < d
c               <= 2*pi) of the domain in spherical coordinates. [radians]
c
c - - Input:    rsun - real*8 value of the radius of the Sun [km].  Usually
c               6.96d5.
c
c - - Input:    et(m,n+1) - real*8 array of the theta
c               component of cE [G km/sec], at the photosphere,
c               on the PE grid.
c
c - - Input:    ep(m+1,n) - real*8 array of the phi component
c               of cE [G km/sec], at the photosphere, and on the TE
c               grid.
c
c - - Input:    detdr(m,n+1) - real*8 array of radial derivative of
c               d cE_theta/dr [G/sec], on PE grid.
c
c - - Input:    depdr(m+1,n) - real*8 array of radial derivative of 
c               d cE_phi/dr [G/sec], on TE grid.
c
c - - Input:    er(m+1,n+1) - real*8 array of the radial component of cE
c               [G km/sec], at the photosphere.  er is on COE grid.
c
c - - Output:   btt(m+1,n),bpt(m,n+1),brt(m,n) - real*8 arrays of the time
c               derivative of the magnetic field components in the theta,phi,
c               and radial directions, respectively, by setting B dot = 
c               - curl cE. [G/sec]
c-
c   PDFI_SS Electric Field Inversion Software
c   http://cgem.ssl.berkeley.edu/cgi-bin/cgem/PDFI_SS/index
c   Copyright (C) 2015-2018 University of California
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
c - - input variables:
c
      integer :: m,n
      real*8 :: a,b,c,d,rsun
      real*8 :: et(m,n+1),ep(m+1,n),er(m+1,n+1)
      real*8 :: detdr(m,n+1),depdr(m+1,n)
c
c - - Output variables:
c
      real*8 :: btt(m+1,n),bpt(m,n+1),brt(m,n)
c
c - - Local variables:
c
      real*8 :: sinth_hlf(m),sinth(m+1)
      real*8 :: dethdr(m,n+1),dephdr(m+1,n)
      real*8 :: dtheta,dphi
      real*8 :: curlt(m+1,n),curlp(m,n+1)
      real*8 :: curlet(m+1,n),curlep(m,n+1),curler(m,n)
c
c - - compute dtheta,phi:
c
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
c
c - - Get values of sin(theta):
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - compute curl E contribution from E_r:
c
      call curl_psi_rhat_co_ss(m,n,er,rsun,sinth,dtheta,dphi,curlt,
     1 curlp)
c
c - - Get dephdr and dethdr from depdr and detdr:
c
      dethdr(:,:)=detdr(:,:)+et(:,:)/rsun
      dephdr(:,:)=depdr(:,:)+ep(:,:)/rsun
c
c - - compute theta and phi components of curl cE (note signs in 2nd terms):
c
      curlet(1:m+1,1:n)=curlt(1:m+1,1:n)-dephdr(1:m+1,1:n)
      curlep(1:m,1:n+1)=curlp(1:m,1:n+1)+dethdr(1:m,1:n+1)
c
c - - compute radial component of curl cE:
c
      call curlh_ce_ss(m,n,et,ep,rsun,sinth,sinth_hlf,dtheta,
     1     dphi,curler)
c
c - - btt,bpt,brt arrays are minus the curl of cE (Faraday's law):
c
      btt(1:m+1,1:n)=-curlet(1:m+1,1:n)
c
c - - If you wanted to do a special treatment of btt at the north and south
c - - poles, instead of letting btt be 0 there, this is where you'd probably
c - - do it.  One idea is to use btt(2:m,1:n) to interpolate values at
c - - btt(1,1:n) and btt(m+1,1:n); interpolating at e.g. the north
c - - pole by averaging
c - - values from the topmost ring of points, using points at 2,j 
c - - and 2,mod(j+n/2,n) (the point on the opposite side of the pole)
c - - to get values at 1,j.  Must remember if doing this that B_theta
c - - will point oppositely at mod(j+n/2,n) than from j. Similarly for 
c - - south pole.  Will also have to deal with interpolation between points 
c - - on opposite side of pole if n is odd. But this idea is not currently 
c - - implemented.
c
      bpt(1:m,1:n+1)=-curlep(1:m,1:n+1)
      brt(1:m,1:n)=-curler(1:m,1:n)
c
c - - we're done
c
      return
      end
