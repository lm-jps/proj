!***********************************************************************
! Module containing the field plus padding for the potential field     *
! calculation.                                                         *
!***********************************************************************
module pad
   integer :: ixmin,ixmax,iymin,iymax,npad
   real,dimension(:,:),allocatable :: blpad
end module pad
