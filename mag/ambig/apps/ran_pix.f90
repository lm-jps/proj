!***********************************************************************************************************************************
! Module containing information related to sequence of random pixels. 
!
! Parameter jump controls the size of the jump between pixels in the sequence of attempted reconfigurations
!
!***********************************************************************************************************************************
module ran_pix
!   integer,parameter :: jump=5,jumpsq=jump*jump,nxmax_jump=nxmax/jump,nymax_jump=nymax/jump
   integer,parameter :: jump=5,jumpsq=jump*jump
!
! Common stuff for the sequence of random pixels
!
   integer :: nxjump,nyjump,ng,ia,ia_prev,ja,ja_prev
!   integer :: ivec(nxmax_jump+1),jvec(nymax_jump+1)
   integer,dimension(:),allocatable :: ivec,jvec
!
! The number of pixels
!
   integer :: nxny
end module ran_pix
!***********************************************************************************************************************************
