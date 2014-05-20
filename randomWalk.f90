!> Provides the routines for each type of random walk.
!!
module randomWalk
  use accuracy
  use io
  use random
  implicit none

contains

  !> Classical Random Walk through network. 
  !!
  !! \details The routine performs a classical random walk through a given
  !! network.
  !!
  !! \param laplacian   Laplacian matrix containing information about the network.
  !! \param vortex      Vector to contain the probability to be in each vortex.
  !! \param timestep    The size of a timestep.
  !! \param maxTimestep The maximum amount of timesteps that will be simulated.
  !! \param time_mode   The time mode of the simulation (either discrete or continous).
  !!
  subroutine CRWalk(laplacian, vortex, timestep, maxTimestep, time_mode)

    integer :: ii, n1
    integer, intent(in) :: maxTimestep
    real(dp), intent(in) :: timestep
    character(*), intent(in) :: time_mode
    real(dp), intent(inout) :: laplacian(:,:)
    real(dp), intent(inout) :: vortex(:)
    
    n1 = size(laplacian, dim=1)
    select case (time_mode)
      case default !! Unknown Mode. End program.
        write(*,*) "Error: Not known modus! Stop program."
        stop
        
      case ("discrete") !! discrete-time Random Walk mode
        laplacian = laplacian * timestep
        do ii = 1, n1
          laplacian(ii,ii) = laplacian(ii,ii) + 1
        end do
        
        do ii = 1, maxTimestep
        vortex = matmul(laplacian,vortex)
        end do
      
      case ("continous") !! continous-time Random Walk mode
      vortex = matmul(laplacian,vortex)
    
    end select
     
     
    do ii = 1, maxTimestep
      vortex = matmul(laplacian,vortex)
    end do

  end subroutine CRWalk

end module randomWalk
