      subroutine interp_data_ss(m,n,btcoe,bpcoe,brcoe,vtcoe,vpcoe,
     1 a,b,bt,bp,br,brte,brpe,bhte,bhpe,vt,vp)
c+
c - -  Purpose:  Interpolate input data from COE grid locations to the staggered
c                grid locations for each physical variable.
c
c - -   Usage:   call interp_data_ss,m,n,btcoe,bpcoe,brcoe,vtcoe,vpcoe,a,b,
c                bt,bp,br,brte,brpe,bhte,bhpe,vt,vp)
c
c - -   Input:   btcoe(m+1,n+1) - real*8 array colatitude component of B on 
c                COE grid (corners including edge corners) [G]
c
c  - -  Input:   bpcoe(m+1,n+1) - real*8 array of longitude (phi) component of B
c                on COE grid [G]
c
c  - -  Input:   brcoe(m+1,n+1) - real*8 array of radial component of B 
c                on the COE grid [G]
c
c  - -  Input:   vtcoe(m+1,n+1) - real*8 array of colatitude component of FLCT 
c                velocity on the COE grid [km/sec]
c
c  - -  Input:   vpcoe(m+1,n+1) - real*8 array of phi component of FLCT
c                velocity on the COE grid [km/sec]
c
c  - -  Input:   a,b - real*8 values of mininum, maximum colatitude of COE grid
c                [radians]
c
c - -   Output:  bt(m+1,n) - real*8 array of colatitude (theta) component of B 
c                interpolated to theta edge locations (TE grid) [G]
c
c - -   Output:  bp(m,n+1) - real*8 array of phi component of B
c                interpolated to phi edge locations (PE grid) [G]
c
c - -   Output:  br(m,n) - real*8 array of radial component of B at cell 
c                center locations (CE grid) [G]
c
c - -   Output:  brte(m+1,n) - real*8 array of radial component of B
c                interpolated to theta edge locations (TE grid) [G]
c
c - -   Output:  brpe(m,n+1) - real*8 array of radial component of B 
c                interpolated to phi edge locations (PE grid) [G]
c
c - -   Output:  bhte(m+1,n) - real*8 array of |Bh| (magnitude of horizontal
c                magnetic field) at theta edge locations (TE grid) [G]
c
c - -   Output:  bhpe(m,n+1) - real*8 array of |Bh| 
c                interpolated to phi edge locations (PE grid) [G]
c
c - -   Output:  vt(m+1,n) - real*8 array of theta component of FLCT velocity,
c                interpolated to theta edge locations (TE grid) [km/sec]
c
c - -   Output:  vp(m,n+1) - real*8 array of phi component of FLCT velocity,
c                interpolated to phi edge locations (PE grid) [km/sec]
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
c - - input variables:
c
      integer :: m,n
      real*8 :: btcoe(m+1,n+1),bpcoe(m+1,n+1),brcoe(m+1,n+1),
     1          vtcoe(m+1,n+1),vpcoe(m+1,n+1),a,b
c
c - - output variables:
c
      real*8 :: bt(m+1,n),bp(m,n+1),br(m,n),brte(m+1,n),brpe(m,n+1),
     1          bhte(m+1,n),bhpe(m,n+1),vt(m+1,n),vp(m,n+1)
c
c - - local variables:
c
      integer :: i,j,jp1,jph,ip1,iph
      real*8 :: btte2,bpte2,btpe2,bppe2
      real*8 sinth(m+1),sinth_hlf(m)
c
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c - - interpolate br at cell centers from edges+corners using flux wtd avg
c
      do i=1,m
         iph=i
         ip1=i+1
         do j=1,n
            jph=j
            jp1=j+1
            br(iph,jph)=(sinth(i)*0.5d0*(brcoe(i,j)+brcoe(i,jp1))
     1              +sinth(ip1)*0.5d0*(brcoe(ip1,j)+brcoe(ip1,jp1)))
     2              /(sinth(i)+sinth(ip1))
         enddo
      enddo
c
c - - interpolate bt,vt,brte,bhte at theta edges:
c
      do i=1,m+1
         do j=1,n
            jph=j
            jp1=j+1
            bt(i,jph)=0.5d0*(btcoe(i,j)+btcoe(i,jp1))
            vt(i,jph)=0.5d0*(vtcoe(i,j)+vtcoe(i,jp1))
            brte(i,jph)=0.5d0*(brcoe(i,j)+brcoe(i,jp1))
            btte2=(0.5d0*(btcoe(i,j)+btcoe(i,jp1)))**2
            bpte2=(0.5d0*(bpcoe(i,j)+bpcoe(i,jp1)))**2
            bhte(i,jph)=sqrt(btte2+bpte2)
         enddo
      enddo
c
c - - interpolate bp,vp,brpe,bhpe at phi edges:
c
      do i=1,m
         iph=i
         ip1=i+1
         do j=1,n+1
            bp(iph,j)=0.5d0*(bpcoe(i,j)+bpcoe(ip1,j))
            vp(iph,j)=0.5d0*(vpcoe(i,j)+vpcoe(ip1,j))
            brpe(iph,j)=0.5d0*(brcoe(i,j)+brcoe(ip1,j))
            btpe2=(0.5d0*(btcoe(i,j)+btcoe(ip1,j)))**2
            bppe2=(0.5d0*(bpcoe(i,j)+bpcoe(ip1,j)))**2
            bhpe(iph,j)=sqrt(btpe2+bppe2)
         enddo
      enddo
c
      return
      end
