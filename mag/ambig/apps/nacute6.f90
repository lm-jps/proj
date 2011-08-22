!***********************************************************************
! subroutine nacute6                                                   *
!    In pixels where the magnitude of the transverse field is below a  *
! given threshold, resolve the 180 degree ambiguity by picking the     *
! direction closest to the average direction of the field at           *
! neighboring pixels. The solution is iterated as the average may      *
! change. Based on the approach in UHIM code.                          *
! This version starts with the pixels which have the most above-       *
! theshold neightbors and works down from there, and allows all pixels *
! not yet visited to change.                                           *
!***********************************************************************

subroutine nacute6(bthresh0,bthresh1)

!-----------------------------------------------------------------------
   use bobs
   use maskvec
   use mgram_data
   use sizes
   use spherical_position
   use verbose

   implicit none

   integer :: i,j,i1d,j1d,nflip,nvisit,inn,jnn,non0
   integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,j1,iflipt
   integer,dimension(:),allocatable :: indx
   integer, dimension(9) :: iexp,ishift,jshift
   integer,dimension(:,:),allocatable :: nnn,nnsav,wt
   real :: Bxs,Bys,bdot,bdotb,bthresh0,bthresh1,bdotm
   real,dimension(:),allocatable :: bt1d
!-----------------------------------------------------------------------
! --> Allocate memory.
   allocate(bt1d(nx*ny),indx(nx*ny))
   allocate(nnn(nx,ny),nnsav(nx,ny),wt(nx,ny))

! --> Determine the number of above threshold neighbors for each pixel. 
! --> Also flag which pixels are themselves above threshold. 

   nnn=0
   do i=2,nx-1
      do j=2,ny-1
         if(bt(i,j).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(j))*cos(t(i,j)))) then
            wt(i,j)=0
            nnn(i-1,j-1)=nnn(i-1,j-1)+1
            nnn(i-1,j)=nnn(i-1,j)+1
            nnn(i-1,j+1)=nnn(i-1,j+1)+1
            nnn(i,j-1)=nnn(i,j-1)+1
            nnn(i,j+1)=nnn(i,j+1)+1
            nnn(i+1,j-1)=nnn(i+1,j-1)+1
            nnn(i+1,j)=nnn(i+1,j)+1
            nnn(i+1,j+1)=nnn(i+1,j+1)+1
         else
            wt(i,j)=1
         endif
      enddo
   enddo

! --> Edge pixels need to be treated separately. 
   do i=2,nx-1
      if(bt(i,1).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(1))*cos(t(i,1)))) then
         wt(i,1)=0
         nnn(i-1,1)=nnn(i-1,1)+1
         nnn(i-1,2)=nnn(i-1,2)+1
         nnn(i,2)=nnn(i,2)+1
         nnn(i+1,1)=nnn(i+1,1)+1
         nnn(i+1,2)=nnn(i+1,2)+1
      else
         wt(i,1)=1
      endif

      if(bt(i,ny).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(ny))*cos(t(i,ny)))) then
         wt(i,ny)=0
         nnn(i-1,ny-1)=nnn(i-1,ny-1)+1
         nnn(i-1,ny)=nnn(i-1,ny)+1
         nnn(i,ny-1)=nnn(i,ny-1)+1
         nnn(i+1,ny-1)=nnn(i+1,ny-1)+1
         nnn(i+1,ny)=nnn(i+1,ny)+1
      else
         wt(i,ny)=1
      endif
   enddo

   do j=2,ny-1
      if(bt(1,j).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(j))*cos(t(1,j)))) then
         wt(1,j)=0
         nnn(1,j-1)=nnn(1,j-1)+1
         nnn(1,j+1)=nnn(1,j+1)+1
         nnn(2,j-1)=nnn(2,j-1)+1
         nnn(2,j)=nnn(2,j)+1
         nnn(2,j+1)=nnn(2,j+1)+1
      else
         wt(1,j)=1
      endif

      if(bt(nx,j).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(j))*cos(t(nx,j)))) then
         wt(nx,j)=0
         nnn(nx,j-1)=nnn(nx,j-1)+1
         nnn(nx,j+1)=nnn(nx,j+1)+1
         nnn(nx-1,j-1)=nnn(nx-1,j-1)+1
         nnn(nx-1,j)=nnn(nx-1,j)+1
         nnn(nx-1,j+1)=nnn(nx-1,j+1)+1
      else
         wt(nx,j)=1
      endif
   enddo

! --> Finally corner pixels. 
   if(bt(1,1).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(1))*cos(t(1,1)))) then
      wt(1,1)=0
      nnn(1,2)=nnn(1,2)+1
      nnn(2,1)=nnn(2,1)+1
      nnn(2,2)=nnn(2,2)+1
   else
      wt(1,1)=1
   endif

   if(bt(1,ny).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(ny))*cos(t(1,ny)))) then
      wt(1,ny)=0
      nnn(1,ny-1)=nnn(1,ny-1)+1
      nnn(2,ny)=nnn(2,ny)+1
      nnn(2,ny-1)=nnn(2,ny-1)+1
   else
      wt(1,ny)=1
   endif

   if(bt(nx,1).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(1))*cos(t(nx,1)))) then
      wt(nx,1)=0
      nnn(nx,2)=nnn(nx,2)+1
      nnn(nx-1,1)=nnn(nx-1,1)+1
      nnn(nx-1,2)=nnn(nx-1,2)+1
   else
      wt(nx,1)=1
   endif

   if(bt(nx,ny).gt.(bthresh1+(bthresh0-bthresh1)*sin(p(ny))*cos(t(nx,ny)))) then
      wt(nx,ny)=0
      nnn(nx,ny-1)=nnn(nx,ny-1)+1
      nnn(nx-1,ny)=nnn(nx-1,ny)+1
      nnn(nx-1,ny-1)=nnn(nx-1,ny-1)+1
   else
      wt(nx,ny)=1
   endif

! --> Create 1-d array of transverse field, rescaled based on variation
! --> of threshold field, and zeroed for pixels outside mask.
   do i=1,nx
      do j=1,ny
         nnsav(i,j)=nnn(i,j)
         i1d=i+(j-1)*nx
         bt1d(i1d)=bt(i,j)*mask(i,j)*(bthresh1+(bthresh0-bthresh1)*sin(p(j))*cos(t(i,j)))
      enddo
   enddo

! --> Arrays containing shifts about present pixel.
   ishift(1)=-1
   ishift(2)=0
   ishift(3)=1
   ishift(4)=-1
   ishift(5)=0
   ishift(6)=1
   ishift(7)=-1
   ishift(8)=0
   ishift(9)=1
   jshift(1)=-1
   jshift(2)=-1
   jshift(3)=-1
   jshift(4)=0
   jshift(5)=0
   jshift(6)=0
   jshift(7)=1
   jshift(8)=1
   jshift(9)=1

! --> Create an index array for pixels sorted on transverse field.

   call sortrx(nx*ny,-bt1d,indx)

! --> How many pixels have a finite transverse field?
   if(bt1d(indx(nx*ny)).gt.0.) then
      non0=nx*ny
   else
      non0=1
      do while(bt1d(indx(non0)).gt.0.)
         non0=non0+1
      enddo
   endif

   jnn=8
   nflip=0
   do while(jnn.gt.0)

! --> Start with pixels with most neighbors above threshold. 
      inn=8

      nvisit=1

! --> Keep iterating as long as any pixels remain with a given value of
! --> nnn.
      do while(inn.ge.jnn)
         nvisit=0

! --> This is a local test, so loop over pixels and apply at each one. 

         do i1d=1,non0
            j1d=indx(i1d)

! --> 2-d indices for this pixel.
            j=1+(j1d-1)/nx
            i=j1d-(j-1)*nx

! --> Only apply to points below threshold and not already visited.
! --> Also must be in mask.

            if(wt(i,j).gt.0.5.and.nnsav(i,j).ge.inn.and.mask(i,j).eq.1) then

! --> Compute dot product of transverse field at this pixel with
! --> transverse field at neighboring pixels. 
! --> Also update neighbors to indicate this pixel has been visited.

! --> If edge pixel, has fewer neighbors.
               if(i.eq.1) then
                  if(j.eq.1) then
                     wt(1,1)=0
                     nnn(1,2)=nnn(1,2)+1
                     nnn(2,1)=nnn(2,1)+1
                     nnn(2,2)=nnn(2,2)+1
                     Bxs=Bx(i+1,j+1)+Bx(i,j+1)+Bx(i+1,j)
                     Bys=By(i+1,j+1)+By(i,j+1)+By(i+1,j)
                  elseif(j.eq.ny) then
                     wt(1,ny)=0
                     nnn(1,ny-1)=nnn(1,ny-1)+1
                     nnn(2,ny)=nnn(2,ny)+1
                     nnn(2,ny-1)=nnn(2,ny-1)+1
                     Bxs=Bx(i+1,j-1)+Bx(i,j-1)+Bx(i+1,j)
                     Bys=By(i+1,j-1)+By(i,j-1)+By(i+1,j)
                  else
                     wt(1,j)=0
                     nnn(1,j-1)=nnn(1,j-1)+1
                     nnn(2,j-1)=nnn(2,j-1)+1
                     nnn(2,j)=nnn(2,j)+1
                     nnn(2,j+1)=nnn(2,j+1)+1
                     nnn(1,j+1)=nnn(1,j+1)+1
                     Bxs=Bx(i,j-1)+Bx(i+1,j-1)+Bx(i+1,j)+Bx(i,j+1)+Bx(i+1,j+1)
                     Bys=By(i,j-1)+By(i+1,j-1)+By(i+1,j)+By(i,j+1)+By(i+1,j+1)
                  endif
                  bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
                  iflipt=0
                  if(bdotb.lt.0.) iflipt=1
               elseif(i.eq.nx) then
                  if(j.eq.1) then
                     wt(nx,1)=0
                     nnn(nx,2)=nnn(nx,2)+1
                     nnn(nx-1,1)=nnn(nx-1,1)+1
                     nnn(nx-1,2)=nnn(nx-2,1)+1
                     Bxs=Bx(i-1,j+1)+Bx(i,j+1)+Bx(i-1,j)
                     Bys=By(i-1,j+1)+By(i,j+1)+By(i-1,j)
                  elseif(j.eq.ny) then
                     wt(nx,ny)=0
                     nnn(nx,ny-1)=nnn(nx,ny-1)+1
                     nnn(nx-1,ny)=nnn(nx-1,ny)+1
                     nnn(nx-1,ny-1)=nnn(nx-1,ny-1)+1
                     Bxs=Bx(i-1,j-1)+Bx(i,j-1)+Bx(i-1,j)
                     Bys=By(i-1,j-1)+By(i,j-1)+By(i-1,j)
                  else
                     wt(nx,j)=0
                     nnn(nx,j-1)=nnn(nx,j-1)+1
                     nnn(nx-1,j-1)=nnn(nx-1,j-1)+1
                     nnn(nx-1,j)=nnn(nx-1,j)+1
                     nnn(nx-1,j+1)=nnn(nx-1,j+1)+1
                     nnn(nx,j+1)=nnn(nx,j+1)+1
                     Bxs=Bx(i,j-1)+Bx(i-1,j-1)+Bx(i-1,j)+Bx(i,j+1)+Bx(i-1,j+1)
                     Bys=By(i,j-1)+By(i-1,j-1)+By(i-1,j)+By(i,j+1)+By(i-1,j+1)
                  endif
                  bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
                  iflipt=0
                  if(bdotb.lt.0.) iflipt=1
               else
                  if(j.eq.1) then
                     wt(i,1)=0
                     nnn(i-1,1)=nnn(i-1,1)+1
                     nnn(i-1,2)=nnn(i,2-1)+1
                     nnn(i,2)=nnn(i,2)+1
                     nnn(i+1,2)=nnn(i,2+1)+1
                     nnn(i+1,1)=nnn(i+1,1)+1
                     Bxs=Bx(i-1,j+1)+Bx(i-1,j)+Bx(i,j+1)+Bx(i+1,j+1)+Bx(i+1,j)
                     Bys=By(i-1,j+1)+By(i-1,j)+By(i,j+1)+By(i+1,j+1)+By(i+1,j)
                     bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
                     iflipt=0
                     if(bdotb.lt.0.) iflipt=1
                  elseif(j.eq.ny) then
                     wt(i,ny)=0
                     nnn(i-1,ny)=nnn(i-1,ny)+1
                     nnn(i-1,ny-1)=nnn(i-1,ny-1)+1
                     nnn(i,ny-1)=nnn(i,ny-1)+1
                     nnn(i+1,ny-1)=nnn(i+1,ny-1)+1
                     nnn(i+1,ny)=nnn(i+1,ny)+1
                     Bxs=Bx(i-1,j-1)+Bx(i-1,j)+Bx(i,j-1)+Bx(i+1,j-1)+Bx(i+1,j)
                     Bys=By(i-1,j-1)+By(i-1,j)+By(i,j-1)+By(i+1,j-1)+By(i+1,j)
                     bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
                     iflipt=0
                     if(bdotb.lt.0.) iflipt=1
                  else
                     nnn(i-1,j-1)=nnn(i-1,j-1)+1
                     nnn(i-1,j)=nnn(i-1,j)+1
                     nnn(i-1,j+1)=nnn(i-1,j+1)+1
                     nnn(i,j-1)=nnn(i,j-1)+1
                     nnn(i,j+1)=nnn(i,j+1)+1
                     nnn(i+1,j-1)=nnn(i+1,j-1)+1
                     nnn(i+1,j)=nnn(i+1,j)+1
                     nnn(i+1,j+1)=nnn(i+1,j+1)+1
! --> Loop over possible directions of transverse field at each pixel. 
! --> If pixel has already been visited or is above threshold, do not
! --> allow it to flip.
                     bdotm=-1.e9
                     do i1=0,wt(i-1,j-1)
                       iexp(1)=i1
                       do i2=0,wt(i,j-1)
                         iexp(2)=i2
                         do i3=0,wt(i+1,j-1)
                           iexp(3)=i3
                           do i4=0,wt(i-1,j)
                             iexp(4)=i4
                             do i5=0,wt(i,j)
                               iexp(5)=i5
                               do i6=0,wt(i+1,j)
                                 iexp(6)=i6
                                 do i7=0,wt(i-1,j+1)
                                   iexp(7)=i7
                                   do i8=0,wt(i,j+1)
                                     iexp(8)=i8
                                     do i9=0,wt(i+1,j+1)
                                       iexp(9)=i9
                                       bdot=0.
                                       do j1=1,9
                                         Bxs=Bx(i+ishift(j1),j+jshift(j1))*(-1.)**iexp(j1)
                                         Bys=By(i+ishift(j1),j+jshift(j1))*(-1.)**iexp(j1)
                                         if(j1.ne.5) then
                                           bdot=bdot+(Bxs*Bx(i,j)+Bys*By(i,j))*(-1.)**iexp(5)
                                         endif
                                       enddo
 

                                       if(bdot.gt.bdotm) then
                                         bdotm=bdot
                                         iflipt=i5
                                       endif
                                     enddo
                                   enddo
                                 enddo
                               enddo
                             enddo
                           enddo
                         enddo
                       enddo
                     enddo
                     wt(i,j)=0

                  endif
               endif

               nvisit=nvisit+1

! --> Flip direction of transverse field if it improves dot product.
               if(iflipt.gt.0.5) then
                  nflip=nflip+1
                  Bx(i,j)=-Bx(i,j)
                  By(i,j)=-By(i,j)
               endif

            endif
         enddo
         if(iverb.eq.2) write(*,*) 'iteration',nflip,inn,nvisit
         if(nvisit.eq.0) then
            inn=inn-1
         else
            inn=8
            do i=1,nx
               do j=1,ny
                  nnsav(i,j)=nnn(i,j)
               enddo
            enddo
         endif
      enddo

      jnn=jnn-1
   enddo
   if(iverb.eq.2) then
      write(*,*) 'total number of flips',nflip
      write(*,*) 'fraction flipped',float(nflip)/float(nx*ny)
   endif

! --> Deallocate memory.
   deallocate(bt1d,indx)
   deallocate(nnn,nnsav,wt)

end subroutine nacute6
