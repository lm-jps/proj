subroutine dilate(naxis1, naxis2, map, dilmap, dilation_param)
!
! PURPOSE: Given input "map" array (Nx,Ny) this performs morphological 
!          dilation, assuming a square dilation element of width
!          (2*dilation_param + 1) pixels.
!
! HISTORY: 2016/10/25, BT Welsch: started
!          2017/04/18, BTW: fixed bug in indexing range
!
  integer*4, intent(in) :: naxis1, naxis2
  integer, intent(in) :: map(naxis1, naxis2) ! input map
  integer, intent(inout) :: dilmap(naxis1, naxis2) ! output dilated map
  integer, intent(in) :: dilation_param ! half-width of dilation element

  integer, allocatable :: shifted(:,:) ! 2D

  integer i,j ! counting indices
  integer mcoli,mcolf,mrowi,mrowf ! array bounds on input map
  integer scoli,scolf,srowi,srowf ! array bounds on shifted map

!  integer naxis1,naxis2
  integer alloc_count, dealloc_count ! bookkeeping of allocated arrays

  alloc_count = 0 
  dealloc_count = 0 

!  naxis1 = size(map, 1)
!  naxis2 = size(map, 2)

  allocate( shifted( naxis1, naxis2))
  alloc_count = alloc_count + 1

  dilmap(:,:) = 0

  do i=-dilation_param,dilation_param   ! Shift map array rows

     shifted(:,:) = 0

     if (i.lt.0) then 
        srowi = -i+1
        srowf = naxis1
        mrowi = 1
        mrowf = naxis1 + i
     else 
        srowi = 1
        srowf = naxis1 - i
        mrowi = i+1
        mrowf = naxis1 
     end if
        
     do j=-dilation_param,dilation_param   ! Shift map array columns

        if (j.lt.0) then 
           scoli = -j+1
           scolf = naxis2
           mcoli = 1
           mcolf = naxis2 + j
        else 
           scoli = 1
           scolf = naxis2 - j
           mcoli = j+1
           mcolf = naxis2 
        end if
        
        !print *,''
        !print *,i,j,shape(shifted)
        !print *,'S row & col:',srowf-srowi,scolf-scoli
        !print *,'M row & col:',mrowf-mrowi,mcolf-mcoli
        !print *,'========================================'
        
        shifted(srowi:srowf,scoli:scolf) = map(mrowi:mrowf,mcoli:mcolf)
        dilmap = dilmap + shifted
        
     end do
  end do

  where (dilmap.gt.0) 
     dilmap = 1
  end where

  !  Cleanup remaining allocated variables prior to exiting
  deallocate( shifted)
  dealloc_count = dealloc_count + 1 

  !print *,'Sums of map & dilmap:',sum(map),sum(dilmap) ! diagnostics

  ! Check for matching allocations and deallocations.
  check_alloc: if (alloc_count /= dealloc_count) then
     print *,'Found allocate/ deallocate mismatch.'
  endif check_alloc

end subroutine dilate
