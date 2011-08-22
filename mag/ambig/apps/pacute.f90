!***********************************************************************
! subroutine pacute                                                    *
!    In pixels where the magnitude of the transverse field is below a  *
! given threshold, resolve the 180 degree ambiguity by picking the     *
! direction closest to the potential field.                            *
!***********************************************************************

subroutine pacute(bthresh)

!-----------------------------------------------------------------------
   use sizes
   use bobs
   use mgram_data
   use pad
   use pot_field
   use trnsfrm
   use verbose

   implicit none

   integer :: i,j,nflip
   real :: bdotb,bthresh
   real,parameter :: pi=3.1415926535897932
!-----------------------------------------------------------------------
! --> This is a local test, so loop over pixels and apply at each one. 

   nflip=0
   do i=1,nx
      do j=1,ny

! --> Only apply to points below threshold.
         if(bt(i,j).le.bthresh) then

! --> Compute dot product of observed transverse field with potential field.
            bdotb=Bpix(i,j)*Bx(i,j)+Bpiy(i,j)*By(i,j)

! --> Flip direction of transverse field if dot product is negative. 
            if(bdotb.lt.0.) then
               nflip=nflip+1
               Bx(i,j)=-Bx(i,j)
               By(i,j)=-By(i,j)
            endif

         endif

      enddo
   enddo
   if(iverb.eq.2) then
      write(*,*) 'number of flips',nflip
      write(*,*) 'fraction flipped',float(nflip)/float(nx*ny)
   endif

end subroutine pacute
