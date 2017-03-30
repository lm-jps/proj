      subroutine interp_var_ss(m,n,btcoe,bpcoe,brcoe,
     1 a,b,bt,bp,br)
c+
c  Purpose: To interpolate a 3-d vector from COE grid to its proper locations
c           on the staggered grid.
c
c  Usage: call interp_data_ss,m,n,btcoe,bpcoe,brcoe,vtcoe,thmin,thmax,vpcoe,
c      bt,bp,br,brte,brpe,bhte,bhpe,vt,vp)
c  Input:  btcoe(m+1,n+1) - Bt at the corners/edges 
c  Input:  bpcoe(m+1,n+1) - Bp at the corners/edges 
c  Input:  brcoe(m+1,n+1) - Br at the corners/edges 
c  Input:  a - mininum colatitude at Brian's interp data points
c  Input:  b - maximum colatitude at Brian's interp data points
c Output:  bt(m+1,n) - Bt at theta edge locations
c Output:  bp(m,n+1) - Bp at phi edge locations
c Output:  br(m,n) - Br at cell center locations
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
      integer :: m,n,i,j,jp1,jph,ip1,iph
c
      real*8 :: btcoe(m+1,n+1),bpcoe(m+1,n+1),brcoe(m+1,n+1),
     1          a,b,
     2          bt(m+1,n),bp(m,n+1),br(m,n)
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
c - - interpolate bt at theta edges:
c
      do i=1,m+1
         do j=1,n
            jph=j
            jp1=j+1
            bt(i,jph)=0.5d0*(btcoe(i,j)+btcoe(i,jp1))
         enddo
      enddo
c
c - - interpolate bp at phi edges:
c
      do i=1,m
         iph=i
         ip1=i+1
         do j=1,n+1
            bp(iph,j)=0.5d0*(bpcoe(i,j)+bpcoe(ip1,j))
         enddo
      enddo
c
      return
      end
