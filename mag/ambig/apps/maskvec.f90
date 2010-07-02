!***********************************************************************************************************************************
! Module containing mask array. 
!***********************************************************************************************************************************
module maskvec
   !integer(kind=1),dimension(:,:),allocatable :: mask
   integer(kind=1),dimension(:,:),allocatable :: tmask
   integer,dimension(:,:),pointer :: mask      ! Pointer to input array bitmap_p
   real,dimension(:,:),pointer :: dBt          ! Pointer to input array dBt_p
end module maskvec
!***********************************************************************************************************************************
