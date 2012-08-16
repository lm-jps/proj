!***********************************************************************
! subroutine nacute6p                                                  *
!    In pixels where the magnitude of the transverse field is below a  *
! given threshold, resolve the 180 degree ambiguity by picking the     *
! direction closest to the average direction of the field at           *
! neighboring pixels. The solution is iterated as the average may      *
! change. Based on the approach in UHIM code.                          *
! This version starts with the pixels which have the most above-       *
! theshold neightbors and works down from there, and allows all pixels *
! not yet visited to change.                                           *
!***********************************************************************

subroutine nacute6p(bthresh)

!-----------------------------------------------------------------------
   use bobs
   use maskvec
   use mgram_data
   use sizes
   use verbose

   implicit none

   integer :: i,j,i1d,j1d,nflip,nvisit,inn,jnn
   integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,j1,iflipt
   integer,dimension(:),allocatable :: indxdo1,indxdo2
   real,dimension(:),allocatable :: btdo
   !integer,dimension(:,:),allocatable :: nnn,wt
! Not clear if the following helps, but it does reduce cache use.
! Assumes that nnn never exceeds 127 (test gave max of 12)
   integer, parameter :: byte = selected_int_kind (2)
   integer(kind=byte),dimension(:,:),allocatable :: nnn,wt
   real :: Bxs,Bys,bdot,bdotb,bthresh,bdotm0,bdotm1
   real :: Bxs0,Bxs1,Bxs2,Bxs3,Bxs4,Bxs6,Bxs7,Bxs8,Bxs9
   real :: Bys0,Bys1,Bys2,Bys3,Bys4,Bys6,Bys7,Bys8,Bys9

!-----------------------------------------------------------------------


! --> Allocate memory.
   allocate(btdo(nx*ny),indxdo1(nx*ny),indxdo2(nx*ny))
   allocate(nnn(nx,ny),wt(nx,ny))

!   open(8,file='nearest.dat')
! --> Determine the number of above threshold neighbors for each pixel. 
! --> Also flag which pixels are themselves above threshold. 

   do j=2,ny-1
      do i=2,nx-1
         if (bt(i,j).gt.bthresh.and.mask(i,j).eq.1) then
            wt(i,j)=0
         else
            wt(i,j)=1
         endif
         nnn(i,j)=0
         if (bt(i-1,j-1).gt.bthresh.and.mask(i-1,j-1).eq.1) nnn(i,j)=nnn(i,j)+1
         if (bt(i,j-1).gt.bthresh.and.mask(i,j-1).eq.1) nnn(i,j)=nnn(i,j)+1
         if (bt(i+1,j-1).gt.bthresh.and.mask(i+1,j-1).eq.1) nnn(i,j)=nnn(i,j)+1
         if (bt(i-1,j).gt.bthresh.and.mask(i-1,j).eq.1) nnn(i,j)=nnn(i,j)+1
         if (bt(i+1,j).gt.bthresh.and.mask(i+1,j).eq.1) nnn(i,j)=nnn(i,j)+1
         if (bt(i-1,j+1).gt.bthresh.and.mask(i-1,j+1).eq.1) nnn(i,j)=nnn(i,j)+1
         if (bt(i,j+1).gt.bthresh.and.mask(i,j+1).eq.1) nnn(i,j)=nnn(i,j)+1
         if (bt(i+1,j+1).gt.bthresh.and.mask(i+1,j+1).eq.1) nnn(i,j)=nnn(i,j)+1
      enddo
   enddo


! --> Edge pixels need to be treated separately. 
   do i=2,nx-1
      if (bt(i,1).gt.bthresh.and.mask(i,1).eq.1) then
         wt(i,1)=0
      else
         wt(i,1)=1
      endif
      nnn(i,1)=0
      if (bt(i-1,1).gt.bthresh.and.mask(i-1,1).eq.1) nnn(i,1)=nnn(i,1)+1
      if (bt(i+1,1).gt.bthresh.and.mask(i+1,1).eq.1) nnn(i,1)=nnn(i,1)+1
      if (bt(i-1,2).gt.bthresh.and.mask(i-1,2).eq.1) nnn(i,1)=nnn(i,1)+1
      if (bt(i,2).gt.bthresh.and.mask(i,2).eq.1) nnn(i,1)=nnn(i,1)+1
      if (bt(i+1,2).gt.bthresh.and.mask(i+1,2).eq.1) nnn(i,1)=nnn(i,1)+1

      if (bt(i,ny).gt.bthresh.and.mask(i,ny).eq.1) then
         wt(i,ny)=0
      else
         wt(i,ny)=1
      endif
      nnn(i,ny)=0
      if (bt(i-1,ny-1).gt.bthresh.and.mask(i-1,ny-1).eq.1) nnn(i,ny)=nnn(i,ny)+1
      if (bt(i,ny-1).gt.bthresh.and.mask(i,ny-1).eq.1) nnn(i,ny)=nnn(i,ny)+1
      if (bt(i+1,ny-1).gt.bthresh.and.mask(i+1,ny-1).eq.1) nnn(i,ny)=nnn(i,ny)+1
      if (bt(i-1,ny).gt.bthresh.and.mask(i-1,ny).eq.1) nnn(i,ny)=nnn(i,ny)+1
      if (bt(i+1,ny).gt.bthresh.and.mask(i+1,ny).eq.1) nnn(i,ny)=nnn(i,ny)+1
   enddo

   do j=2,ny-1
      if (bt(1,j).gt.bthresh.and.mask(1,j).eq.1) then
         wt(1,j)=0
      else
         wt(1,j)=1
      endif
      nnn(1,j)=0
      if (bt(1,j-1).gt.bthresh.and.mask(1,j-1).eq.1) nnn(1,j)=nnn(1,j)+1
      if (bt(1,j+1).gt.bthresh.and.mask(1,j+1).eq.1) nnn(1,j)=nnn(1,j)+1
      if (bt(2,j-1).gt.bthresh.and.mask(2,j-1).eq.1) nnn(1,j)=nnn(1,j)+1
      if (bt(2,j).gt.bthresh.and.mask(2,j).eq.1) nnn(1,j)=nnn(1,j)+1
      if (bt(2,j+1).gt.bthresh.and.mask(2,j+1).eq.1) nnn(1,j)=nnn(1,j)+1

      if (bt(nx,j).gt.bthresh.and.mask(nx,j).eq.1) then
         wt(nx,j)=0
      else
         wt(nx,j)=1
      endif
      nnn(nx,j)=0
      if (bt(nx-1,j-1).gt.bthresh.and.mask(nx-1,j-1).eq.1) nnn(nx,j)=nnn(nx,j)+1
      if (bt(nx-1,j).gt.bthresh.and.mask(nx-1,j).eq.1) nnn(nx,j)=nnn(nx,j)+1
      if (bt(nx-1,j+1).gt.bthresh.and.mask(nx-1,j+1).eq.1) nnn(nx,j)=nnn(nx,j)+1
      if (bt(nx,j-1).gt.bthresh.and.mask(nx,j-1).eq.1) nnn(nx,j)=nnn(nx,j)+1
      if (bt(nx,j+1).gt.bthresh.and.mask(nx,j+1).eq.1) nnn(nx,j)=nnn(nx,j)+1
   enddo

! --> Finally corner pixels. 
   if (bt(1,1).gt.bthresh.and.mask(1,1).eq.1) then
      wt(1,1)=0
   else
      wt(1,1)=1
   endif
   nnn(1,1)=0
   if (bt(1,2).gt.bthresh.and.mask(1,2).eq.1) nnn(1,1)=nnn(1,1)+1
   if (bt(2,1).gt.bthresh.and.mask(2,1).eq.1) nnn(1,1)=nnn(1,1)+1
   if (bt(2,2).gt.bthresh.and.mask(2,2).eq.1) nnn(1,1)=nnn(1,1)+1

   if (bt(1,ny).gt.bthresh.and.mask(1,ny).eq.1) then
      wt(1,ny)=0
   else
      wt(1,ny)=1
   endif
   nnn(1,ny)=0
   if (bt(1,ny-1).gt.bthresh.and.mask(1,ny-1).eq.1) nnn(1,ny)=nnn(1,ny)+1
   if (bt(2,ny).gt.bthresh.and.mask(2,ny).eq.1) nnn(1,ny)=nnn(1,ny)+1
   if (bt(2,ny-1).gt.bthresh.and.mask(2,ny-1).eq.1) nnn(1,ny)=nnn(1,ny)+1

   if (bt(nx,1).gt.bthresh.and.mask(nx,1).eq.1) then
      wt(nx,1)=0
   else
      wt(nx,1)=1
   endif
   nnn(nx,1)=0
   if (bt(nx-1,1).gt.bthresh.and.mask(nx-1,1).eq.1) nnn(nx,1)=nnn(nx,1)+1
   if (bt(nx,2).gt.bthresh.and.mask(nx,2).eq.1) nnn(nx,1)=nnn(nx,1)+1
   if (bt(nx-1,2).gt.bthresh.and.mask(nx-1,2).eq.1) nnn(nx,1)=nnn(nx,1)+1

   if (bt(nx,ny).gt.bthresh.and.mask(nx,ny).eq.1) then
      wt(nx,ny)=0
   else
      wt(nx,ny)=1
   endif
   nnn(nx,ny)=0
   if (bt(nx-1,ny).gt.bthresh.and.mask(nx-1,ny).eq.1) nnn(nx,ny)=nnn(nx,ny)+1
   if (bt(nx,ny-1).gt.bthresh.and.mask(nx,ny-1).eq.1) nnn(nx,ny)=nnn(nx,ny)+1
   if (bt(nx-1,ny-1).gt.bthresh.and.mask(nx-1,ny-1).eq.1) nnn(nx,ny)=nnn(nx,ny)+1

! --> write out number of neighbors
!   do j=1,ny
!      write(8,*) (nnn(i,j),i=1,nx)
!   enddo


   jnn=8
   nflip=0
   do while(jnn.gt.0)

! --> Start with pixels with most neighbors above threshold. 
      inn=8

      nvisit=1

! --> Keep iterating as long as any pixels remain with a given value of
! --> nnn.
      do while(inn.ge.jnn)

! At this point the pixels to be visited is known, is based only on wt and nnsav
! and nnsav=nnn. Order is also known based on bt.
! So let's find those pixels and sort them
! Get rid of nnsav (which used to keep a copy of nnn)
         nvisit=0
         do j=1,ny
            do i=1,nx
! --> Only apply to points below threshold and not already visited.
! Split if statement in two to avoid cache misses on nnn(i,j) access.
               if (wt(i,j).ne.0) then
                 if (nnn(i,j).ge.inn) then
                   nvisit=nvisit+1
                   indxdo1(nvisit)=i+(j-1)*nx
                   btdo(nvisit)=-bt(i,j)
                 endif
               endif
            enddo
         enddo
         call sortrx(nvisit,btdo,indxdo2)

         do i1d=1,nvisit
            j1d=indxdo1(indxdo2(i1d))


! --> 2-d indices for this pixel.
            j=1+(j1d-1)/nx
            i=j1d-(j-1)*nx

!write(8,*) i,j
! --> Compute dot product of transverse field at this pixel with
! --> transverse field at neighboring pixels. 
! --> Also update neighbors to indicate this pixel has been visited.

! --> If edge pixel, has fewer neighbors.
            if (i.eq.1) then
               if (j.eq.1) then
                  nnn(1,2)=nnn(1,2)+1
                  nnn(2,1)=nnn(2,1)+1
                  nnn(2,2)=nnn(2,2)+1
                  Bxs=Bx(i+1,j+1)+Bx(i,j+1)+Bx(i+1,j)
                  Bys=By(i+1,j+1)+By(i,j+1)+By(i+1,j)
               elseif (j.eq.ny) then
                  nnn(1,ny-1)=nnn(1,ny-1)+1
                  nnn(2,ny)=nnn(2,ny)+1
                  nnn(2,ny-1)=nnn(2,ny-1)+1
                  Bxs=Bx(i+1,j-1)+Bx(i,j-1)+Bx(i+1,j)
                  Bys=By(i+1,j-1)+By(i,j-1)+By(i+1,j)
               else
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
               if (bdotb.lt.0.) iflipt=1
            elseif (i.eq.nx) then
               if (j.eq.1) then
                  nnn(nx,2)=nnn(nx,2)+1
                  nnn(nx-1,1)=nnn(nx-1,1)+1
                  nnn(nx-1,2)=nnn(nx-2,1)+1
                  Bxs=Bx(i-1,j+1)+Bx(i,j+1)+Bx(i-1,j)
                  Bys=By(i-1,j+1)+By(i,j+1)+By(i-1,j)
               elseif (j.eq.ny) then
                  nnn(nx,ny-1)=nnn(nx,ny-1)+1
                  nnn(nx-1,ny)=nnn(nx-1,ny)+1
                  nnn(nx-1,ny-1)=nnn(nx-1,ny-1)+1
                  Bxs=Bx(i-1,j-1)+Bx(i,j-1)+Bx(i-1,j)
                  Bys=By(i-1,j-1)+By(i,j-1)+By(i-1,j)
               else
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
               if (bdotb.lt.0.) iflipt=1
            else
               if (j.eq.1) then
                  nnn(i-1,1)=nnn(i-1,1)+1
                  nnn(i-1,2)=nnn(i,2-1)+1
                  nnn(i,2)=nnn(i,2)+1
                  nnn(i+1,2)=nnn(i,2+1)+1
                  nnn(i+1,1)=nnn(i+1,1)+1
                  Bxs=Bx(i-1,j+1)+Bx(i-1,j)+Bx(i,j+1)+Bx(i+1,j+1)+Bx(i+1,j)
                  Bys=By(i-1,j+1)+By(i-1,j)+By(i,j+1)+By(i+1,j+1)+By(i+1,j)
                  bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
                  iflipt=0
                  if (bdotb.lt.0.) iflipt=1
               elseif (j.eq.ny) then
                  nnn(i-1,ny)=nnn(i-1,ny)+1
                  nnn(i-1,ny-1)=nnn(i-1,ny-1)+1
                  nnn(i,ny-1)=nnn(i,ny-1)+1
                  nnn(i+1,ny-1)=nnn(i+1,ny-1)+1
                  nnn(i+1,ny)=nnn(i+1,ny)+1
                  Bxs=Bx(i-1,j-1)+Bx(i-1,j)+Bx(i,j-1)+Bx(i+1,j-1)+Bx(i+1,j)
                  Bys=By(i-1,j-1)+By(i-1,j)+By(i,j-1)+By(i+1,j-1)+By(i+1,j)
                  bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
                  iflipt=0
                  if (bdotb.lt.0.) iflipt=1
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
                  bdotm0=-1.e9
                  bdotm1=1.e9
                  do i1=0,wt(i-1,j-1)
                    if (i1.eq.0) then
                      Bxs1=Bx(i-1,j-1)
                      Bys1=By(i-1,j-1)
                    else
                      Bxs1=-Bx(i-1,j-1)
                      Bys1=-By(i-1,j-1)
                    endif
                    do i2=0,wt(i,j-1)
                      if (i2.eq.0) then
                        Bxs2=Bxs1+Bx(i,j-1)
                        Bys2=Bys1+By(i,j-1)
                      else
                        Bxs2=Bxs1-Bx(i,j-1)
                        Bys2=Bys1-By(i,j-1)
                      endif
                      do i3=0,wt(i+1,j-1)
                        if (i3.eq.0) then
                          Bxs3=Bxs2+Bx(i+1,j-1)
                          Bys3=Bys2+By(i+1,j-1)
                        else
                          Bxs3=Bxs2-Bx(i+1,j-1)
                          Bys3=Bys2-By(i+1,j-1)
                        endif
                        do i4=0,wt(i-1,j)
                          if (i4.eq.0) then
                            Bxs4=Bxs3+Bx(i-1,j)
                            Bys4=Bys3+By(i-1,j)
                          else
                            Bxs4=Bxs3-Bx(i-1,j)
                            Bys4=Bys3-By(i-1,j)
                          endif
                          do i6=0,wt(i+1,j)
                            if (i6.eq.0) then
                              Bxs6=Bxs4+Bx(i+1,j)
                              Bys6=Bys4+By(i+1,j)
                            else
                              Bxs6=Bxs4-Bx(i+1,j)
                              Bys6=Bys4-By(i+1,j)
                            endif
                            do i7=0,wt(i-1,j+1)
                              if (i7.eq.0) then
                                Bxs7=Bxs6+Bx(i-1,j+1)
                                Bys7=Bys6+By(i-1,j+1)
                              else
                                Bxs7=Bxs6-Bx(i-1,j+1)
                                Bys7=Bys6-By(i-1,j+1)
                              endif
                              do i8=0,wt(i,j+1)
                                if (i8.eq.0) then
                                  Bxs8=Bxs7+Bx(i,j+1)
                                  Bys8=Bys7+By(i,j+1)
                                else
                                  Bxs8=Bxs7-Bx(i,j+1)
                                  Bys8=Bys7-By(i,j+1)
                                endif
                                do i9=0,wt(i+1,j+1)
                                  if (i9.eq.0) then
                                    Bxs9=Bxs8+Bx(i+1,j+1)
                                    Bys9=Bys8+By(i+1,j+1)
                                  else
                                    Bxs9=Bxs8-Bx(i+1,j+1)
                                    Bys9=Bys8-By(i+1,j+1)
                                  endif

                                  bdot=Bxs9*Bx(i,j)+Bys9*By(i,j)
                                  if (bdot.gt.bdotm0) then
                                    bdotm0=bdot ! Most positive so far
                                  endif
                                  if (bdot.lt.bdotm1) then
                                    bdotm1=bdot ! Most negative so far
                                  endif
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                  iflipt=0
                  if (abs(bdotm1).gt.bdotm0) iflipt=1

               endif
            endif

! --> Flip direction of transverse field if it improves dot product.
            if (iflipt.ne.0) then
               nflip=nflip+1
               Bx(i,j)=-Bx(i,j)
               By(i,j)=-By(i,j)
            endif
            wt(i,j)=0
         enddo
         if (iverb.eq.2) write(*,*) 'iteration',nflip,inn,nvisit
         if (nvisit.eq.0) then
            inn=inn-1
         else
            inn=8
         endif
      enddo

      jnn=jnn-1

   enddo

   if (iverb.eq.2) then
      write(*,*) 'total number of flips',nflip
      write(*,*) 'fraction flipped',float(nflip)/float(nx*ny)
   endif

! --> Deallocate memory.
   deallocate(btdo,indxdo1,indxdo2)
   deallocate(nnn,wt)


end subroutine nacute6p
