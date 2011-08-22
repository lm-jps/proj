!***********************************************************************
! subroutine nacute5                                                   *
!    In pixels where the magnitude of the transverse field is below a  *
! given threshold, resolve the 180 degree ambiguity by picking the     *
! direction closest to the average direction of the field at           *
! neighboring pixels. The solution is iterated as the average may      *
! change. Based on the approach in UHIM code.                          *
! This version starts with the pixel with the strongest field and      *
! works radially out from there.                                       *
!***********************************************************************

subroutine nacute5(bthresh)

!-----------------------------------------------------------------------
   use sizes
   use bobs
   use maskvec
   use mgram_data
   use verbose

   implicit none

   integer :: i,j,i1d,j1d,iter,niter,nflip,nflip1,imax,jmax
   integer,dimension(:),allocatable :: indx
   real :: Bxs,Bys,bdotb,bthresh,bmax,bmag,dmax
   real,dimension(:),allocatable :: bt1d,dist1d
!-----------------------------------------------------------------------
!   open(3,file='nacute.dat')
! --> Set up equivalent 1-d arrays.

   allocate(bt1d(nx*ny),indx(nx*ny),dist1d(nx*ny))

! --> Determine the pixel with the strongest field. 
! --> Along the way, also convert 2-d arrays to 1-d arrays.

   bmax=-1.
   do i=1,nx
      do j=1,ny
         bmag=Bx(i,j)**2+By(i,j)**2+Bz(i,j)**2
         if(bmag.gt.bmax) then
           bmax=bmag
           imax=i
           jmax=j
         endif
         i1d=i+(j-1)*nx
         bt1d(i1d)=bt(i,j)
      enddo
   enddo
   if(iverb.eq.2) write(*,*) 'max field at',imax,jmax,bmax

! --> Calculate distance from each pixel to strongest field pixel.
   dmax=float(nx)**2+float(ny)**2
   do i=1,nx
      do j=1,ny
         i1d=i+(j-1)*nx
         dist1d(i1d)=sqrt(float(i-imax)**2+float(j-jmax)**2)
      enddo
   enddo

! --> Create an index array for pixels sorted on distance.

   call sortrx(nx*ny,dist1d,indx)

! --> Maximum number of iterations. 
   niter=50

   nflip=0
   nflip1=1
   iter=1
   do while(nflip1.gt.0.and.iter.lt.niter)

! --> This is a local test, so loop over pixels and apply at each one. 

      nflip1=0
      do i1d=1,nx*ny
         j1d=indx(i1d)

! --> 2-d indices for this pixel.
         j=1+(j1d-1)/nx
         i=j1d-(j-1)*nx

! --> Only apply to points below threshold or not in the mask.

         if(bt1d(j1d).le.bthresh.or.mask(i,j).eq.0) then

! --> Compute dot product of transverse field at this pixel with
! --> transverse field at neighboring pixels. 

! --> If edge pixel, has fewer neighbors.
            if(i.eq.1) then
               if(j.eq.1) then
                  Bxs=Bx(i+1,j+1)+Bx(i,j+1)+Bx(i+1,j)
                  Bys=By(i+1,j+1)+By(i,j+1)+By(i+1,j)
               elseif(j.eq.ny) then
                  Bxs=Bx(i+1,j-1)+Bx(i,j-1)+Bx(i+1,j)
                  Bys=By(i+1,j-1)+By(i,j-1)+By(i+1,j)
               else
                  Bxs=Bx(i,j-1)+Bx(i+1,j-1)+Bx(i+1,j)+Bx(i,j+1)+Bx(i+1,j+1)
                  Bys=By(i,j-1)+By(i+1,j-1)+By(i+1,j)+By(i,j+1)+By(i+1,j+1)
               endif
            elseif(i.eq.nx) then
               if(j.eq.1) then
                  Bxs=Bx(i-1,j+1)+Bx(i,j+1)+Bx(i-1,j)
                  Bys=By(i-1,j+1)+By(i,j+1)+By(i-1,j)
               elseif(j.eq.ny) then
                  Bxs=Bx(i-1,j-1)+Bx(i,j-1)+Bx(i-1,j)
                  Bys=By(i-1,j-1)+By(i,j-1)+By(i-1,j)
               else
                  Bxs=Bx(i,j-1)+Bx(i-1,j-1)+Bx(i-1,j)+Bx(i,j+1)+Bx(i-1,j+1)
                  Bys=By(i,j-1)+By(i-1,j-1)+By(i-1,j)+By(i,j+1)+By(i-1,j+1)
               endif
            else
               if(j.eq.1) then
                  Bxs=Bx(i-1,j+1)+Bx(i-1,j)+Bx(i,j+1)+Bx(i+1,j+1)+Bx(i+1,j)
                  Bys=By(i-1,j+1)+By(i-1,j)+By(i,j+1)+By(i+1,j+1)+By(i+1,j)
               elseif(j.eq.ny) then
                  Bxs=Bx(i-1,j-1)+Bx(i-1,j)+Bx(i,j-1)+Bx(i+1,j-1)+Bx(i+1,j)
                  Bys=By(i-1,j-1)+By(i-1,j)+By(i,j-1)+By(i+1,j-1)+By(i+1,j)
               else
                  Bxs=Bx(i-1,j-1)+Bx(i,j-1)+Bx(i+1,j-1)+Bx(i-1,j+1)+Bx(i,j+1)+Bx(i+1,j+1)+Bx(i-1,j)+Bx(i+1,j)
                  Bys=By(i-1,j-1)+By(i,j-1)+By(i+1,j-1)+By(i-1,j+1)+By(i,j+1)+By(i+1,j+1)+By(i-1,j)+By(i+1,j)
               endif
            endif

            bdotb=Bxs*Bx(i,j)+Bys*By(i,j)

! --> Flip direction of transverse field if dot product is negative.
            if(bdotb.lt.0.) then
               nflip1=nflip1+1
               Bx(i,j)=-Bx(i,j)
               By(i,j)=-By(i,j)
            endif

         endif

      enddo
      iter=iter+1
      nflip=nflip+nflip1
      if(iverb.eq.2) write(*,*) 'iteration',iter,nflip1,nflip
   enddo
   if(iverb.eq.2) write(*,*) 'total number of flips',nflip
   if(iverb.eq.2) write(*,*) 'fraction flipped',float(nflip)/float(nx*ny)

! --> Deallocate memory.
   deallocate(bt1d,indx,dist1d)

end subroutine nacute5
