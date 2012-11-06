!***********************************************************************************************************************************
subroutine get_ran_pix
!===================================================================================================================================
! This subroutine determines the coordinates for a sequence of pixels with a random start point
!
! 2009,2010, Ashley Crouch, ash@cora.nwra.com
!===================================================================================================================================
   use sizes
   use bounds
   use ran_pix
   use ranseed
   use maskvec

   implicit none

   integer,parameter :: max_tries=50

   integer :: i,j,ii,jj,ntries,get_pix,nxg,nyg
   real :: ran3
   integer,dimension(:),allocatable :: ivec1,jvec1
!
! Allocate memory for temporary arrays.
!
   allocate(ivec1((nx/jump)+1),jvec1((ny/jump)+1))

   if ((nxjump.eq.1).and.(nyjump.eq.1)) then
!
! Select a single random pixel
!
      i=int(nxny*ran3(seed))
      if (i.ge.nxny) i=nxny-1
      j=i/dnx
      i=i-j*dnx+nxa
      j=j+nya

      ng=1
      ivec(ng)=i
      jvec(ng)=j
   else
!
! Select a random start point and avoid the grid from the previous sequence
!
      ntries=0
      get_pix=1
      do while ((get_pix.eq.1).and.(ntries.le.max_tries))
         ia=int(ran3(seed)*jumpsq)
         if (ia.eq.jumpsq) ia=jumpsq-1
         ja=ia/jump
         ia=ia-ja*jump+nxa
         ja=ja+nya
         if ((ia.ne.ia_prev).and.(ja.ne.ja_prev)) get_pix=0
         ntries=ntries+1
      enddo
      ia_prev=ia
      ja_prev=ja
!
! Construct the sequence of pixels with random start point and separation: jump
!
      nxg=0
      ii=1
      get_pix=1
      do while (get_pix.eq.1)
         i=ia+(ii-1)*jump
         if (i.le.nxb) then
            nxg=nxg+1
            ivec1(nxg)=i
         else
            get_pix=0
         endif
         ii=ii+1
      enddo
      nyg=0
      jj=1
      get_pix=1
      do while (get_pix.eq.1)
         j=ja+(jj-1)*jump
         if (j.le.nyb) then
            nyg=nyg+1
            jvec1(nyg)=j
         else
            get_pix=0
         endif
         jj=jj+1
      enddo
!
! Set up ivec and jvec so that we only visit pixels that have mask.ge.1
!
      ng=0
      do i=1,nxg
         do j=1,nyg
            if (mask(ivec1(i),jvec1(j)).ge.1) then
               ng=ng+1
               ivec(ng)=ivec1(i)
               jvec(ng)=jvec1(j)
            endif
         enddo
      enddo
   endif
!
! Deallocate memory for temporary arrays.
!
   deallocate(ivec1,jvec1)

end subroutine get_ran_pix
!***********************************************************************************************************************************
