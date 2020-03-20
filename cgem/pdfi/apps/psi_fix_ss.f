      subroutine psi_fix_ss(m,n,p,bcn,a,b,c,d,rsun,rssmrs,psi3d,mflux)
c
c+ - Purpose: Remove 1/r artifact from the psi3d array, and return
c             corrected psi3d array in place.
c
c - - Method: Fit the 1/r dependence by evaluating uniform Br contribution
c             at photosphere, and then subtract off the artifact term from
c             the array.  Re-uses much of the bpot_psi_ss code.
c
c   - Usage:  call psi_fix_ss(m,n,p,bcn,a,b,c,d,rsun,rssmrs,psi3d,mflux)
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
c - - Input:  psi3d(m,n,p+1) - real*8 3D array of psi (scalar potential)
c             for a potential magnetic field, evaluated at at cell centers in 
c             theta,phi, and at radial shells. [G km]
c
c - - Input:  mflux - real*8 value of the net radial magnetic flux computed
c             from subroutine mflux_ss [G km^2]
c
c - - Output: psi3d(m,n,p+1) - real*8 array of psi (scalar potential) for
c             a potential magnetic field, with 1/r term removed. [G km]
c             The input array is replaced by the output array.
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
      real*8 :: mflux,psi3d(m,n,p+1)
c
c - - output variables:
c
c     psi3d(m,n,p+1) (already declared in the input array section)     
c
c - - Internal variables and arrays:
c
      integer :: q,qph,iph,jph
      real*8 :: dtheta,dphi,delr,pi,dum,badmflux,brbad
      real*8 :: wtop,wbot,wtot,wtotinv,b0,area
      real*8 :: darea,areass,psissbar
      real*8 :: r(p+1),rce(p),sinth(m+1),sinth_hlf(m)
      real*8 :: psi2d1(m+2,n+2)
      real*8 :: bt0(m+1,n),bt1(m+1,n),bp0(m,n+1),bp1(m,n+1)
      real*8 :: br(m,n),brtop(m,n),brbot(m,n)
      real*8 :: bdas(n),bdbs(n),bdcs(m),bdds(m)
      real*8 :: gradt(m+1,n),gradp(m,n+1)
      real*8 :: btpot(m+1,n),bppot(m,n+1)
      real*8 :: psifix(m,n)
c
c - - function pimach in fishpack (fftpack) that computes pi
c
      real*8 :: pimach
c
c - - sdf variables (temporary)
c     integer*8 :: dims(20)
c
      if((bcn .ne. 0) .and. (bcn .ne. 3)) then
         write(6,*) 'psi_fix_ss: Illegal bcn = ',bcn, ' exiting'
         stop
      endif
c
      pi=pimach(dum)
c - - Define dtheta,dphi,dphifft,delr:
      dtheta=(b-a)/dble(m)
      dphi=(d-c)/dble(n)
c
c - - set brbad to 0 initially. (brbad will be net radial flux
c - - from possible artifact in psi solution later on).  
c
      brbad=0.d0
c
c - - calculate area of the spherical wedge at photosphere:
c
      area=(d-c)*(cos(a)-cos(b))*rsun**2
c
c - - get value of b0 by dividing mflux by area:
c
      b0=mflux/area
c
c - - radial spacing:
c
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
c - - area at source-surface:
c
      areass=(r(p+1)**2/rsun**2)*area
c
c - - compute mean of psi at source-surface:
c
      psissbar=0.d0
      do iph=1,m
         darea=(r(p+1)**2)*dphi*dtheta*sinth_hlf(iph)
         do jph=1,n
            psissbar=psissbar+psi3d(iph,jph,p+1)*darea
         enddo
      enddo
c
      psissbar=psissbar/areass
c
c - - Set homogenous boundary condition arrays:
c
      bdas(1:n)=0.d0
      bdbs(1:n)=0.d0
      bdcs(1:m)=0.d0
      bdds(1:m)=0.d0
c
c - - Next get psi at photosphere, and add ghost zones:
c 
      q=1
      psi2d1(2:m+1,2:n+1)=psi3d(:,:,q)
      psi2d1(1,2:n+1)=psi2d1(2,2:n+1)-1.d0*dtheta*bdas(1:n)
      psi2d1(m+2,2:n+1)=psi2d1(m+1,2:n+1)+1.d0*dtheta*bdbs(1:n)
c
c - -    Neumann BC in phi:
c
      if(bcn .eq. 3) then
         if(q .eq. 1) then
           psi2d1(2:m+1,1)=psi2d1(2:m+1,2)-1.d0*dphi*bdcs(1:m)
           psi2d1(2:m+1,n+2)=psi2d1(2:m+1,n+1)+1.d0*dphi*bdds(1:m)
         else
           psi2d1(2:m+1,1)=psi2d1(2:m+1,3)-2.d0*dphi*bdcs(1:m)
           psi2d1(2:m+1,n+2)=psi2d1(2:m+1,n)+2.d0*dphi*bdds(1:m)
         endif
      else
c
c - -       Periodic BC in phi:
c
         psi2d1(2:m+1,1)=psi2d1(2:m+1,n+1)
         psi2d1(2:m+1,n+2)=psi2d1(2:m+1,2)
      endif
c
c - - Take minus gradient of psi at photosphere:
c
      call gradh_ce_ss(m,n,psi2d1,r(q),sinth_hlf,dtheta,
     1        dphi,gradt,gradp)
      bt1(1:m+1,1:n)=-gradt(1:m+1,1:n)
      bp1(1:m,1:n+1)=-gradp(1:m,1:n+1)
c
c - - Now, ready for big loop over r index q from 2,p
c
      do qph=1,p
         q=qph+1
c
c - - set bottom layer values of bt,bp equal to former top layers of bt, bp
c
         bt0(:,:)=bt1(:,:)
         bp0(:,:)=bp1(:,:)
c
c - - calculate B_r at mid-voxel locations by taking minus r deriv of psi:
c - - (note that here in psi_fix_ss, we only need to evaluate br for q .eq. 2)
c
         if (q .eq. 2) then
            br(:,:) = -(psi3d(:,:,q)-psi3d(:,:,q-1))/delr
c
c
c - - calculate ghost zones for top layer of psi, calculate gradient of psi
c - - (repeat of code done at photosphere, except that q > 1)
c
            psi2d1(2:m+1,2:n+1)=psi3d(:,:,q)
            psi2d1(1,2:n+1)=psi2d1(2,2:n+1)-1.d0*dtheta*bdas(1:n)
            psi2d1(m+2,2:n+1)=psi2d1(m+1,2:n+1)+1.d0*dtheta*bdbs(1:n)
c
c - -    Neumann BC in phi:
c
            if(bcn .eq. 3) then
               if(q .eq. 1) then
                 psi2d1(2:m+1,1)=psi2d1(2:m+1,2)-1.d0*dphi*bdcs(1:m)
                 psi2d1(2:m+1,n+2)=psi2d1(2:m+1,n+1)+1.d0*dphi*bdds(1:m)
               else
                 psi2d1(2:m+1,1)=psi2d1(2:m+1,3)-2.d0*dphi*bdcs(1:m)
                 psi2d1(2:m+1,n+2)=psi2d1(2:m+1,n)+2.d0*dphi*bdds(1:m)
               endif
            else
c
c - -       Periodic BC in phi:
c
               psi2d1(2:m+1,1)=psi2d1(2:m+1,n+1)
               psi2d1(2:m+1,n+2)=psi2d1(2:m+1,2)
            endif
c
c - - Take minus horizontal gradient of psi:
c
            call gradh_ce_ss(m,n,psi2d1,r(q),sinth_hlf,dtheta,
     1           dphi,gradt,gradp)
            bt1(1:m+1,1:n)=-gradt(1:m+1,1:n)
            bp1(1:m,1:n+1)=-gradp(1:m,1:n+1)
c
c - - now calculate weights for interpolation of layer1 and layer0 horizontal
c - - fields to mid-voxel locations:
c
            wtop = r(q)**2-rce(qph)**2
            wbot = rce(qph)**2-r(q-1)**2
            wtot = r(q)**2 - r(q-1)**2
            wtotinv=1.d0/wtot
            btpot(1:m+1,1:n)=wtotinv*(wtop*bt1(:,:)+wbot*bt0(:,:))
            bppot(1:m,1:n+1)=wtotinv*(wtop*bp1(:,:)+wbot*bp0(:,:))
c
c - - Now call br_voxels3d_ss to get B_r on bottom and top layers:
c - - (apart from last layer, the brtop values will be ignored)
c - - Since we're not computing Br anywhere except at phot, only need
c - - to call br_voxels3d_ss once, when q .eq. 2, here in psi_fix_ss:
c
c
c - - if spurious radial net flux is present, find out what it is:
c
            call br_voxels3d_ss(m,n,rce(qph),sinth,sinth_hlf,dtheta,
     1        dphi,btpot(:,:),bppot(:,:),br,delr,brtop,brbot)
            call mflux_ss(m,n,a,b,c,d,r(q-1),brbot,badmflux)
            brbad=badmflux/area
         endif
c
c - - Now compute corrected psi in this layer (=psifix):
c - - (don't forget to add net flux term, and to remove spurious flux)
c
         psifix(1:m,1:n)=psi3d(1:m,1:n,q-1)+(b0-brbad)*
     1             (rsun+rssmrs-r(q-1))*rsun**2/((rsun+rssmrs)*r(q-1))
     2             -psissbar
c
c - - as last step, replace psi3d at q-1:
c
         psi3d(1:m,1:n,q-1)=psifix(1:m,1:n)
      enddo         
c
c - - Loop over qph is done.  We still need to set psi3d at topmost layer:
c
      q=p+1
      psifix(1:m,1:n)=psi3d(1:m,1:n,q)+(b0-brbad)*
     1       (rsun+rssmrs-r(q))*rsun**2/((rsun+rssmrs)*r(q))
     2       -psissbar
      psi3d(1:m,1:n,q)=psifix(1:m,1:n)
c
c - - we're done:
c
      return
      end
