module randomWalk
  use accuracy
  use io
  use random
  implicit none
  
contains

  !> Classical Random Walk through Network 
  !!
  subroutine CRWalk(laplacian,vortex,maxTimestep)

    integer :: ii
    real(dp), intent(in) :: laplacian(:,:)
    integer, intent(in) :: maxTimestep
    real(dp), intent(inout) :: vortex(:)
    
    do ii = 1, maxTimestep
      vortex = matmul(laplacian,vortex)
    end do
    
    end subroutine CRWalk
    
end module randomWalk