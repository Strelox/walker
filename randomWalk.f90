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
  !! \param start	Start vortex
  !! \param timestep    The size of a timestep.
  !! \param maxTimestep The maximum amount of timesteps that will be simulated.
  !! \param time_mode   The time mode of the simulation (either discrete or continous).
  !!
  subroutine CRWalk(laplacian, vortex, start, timestep, maxTimestep, time_mode, input)

    integer :: ii, n1
    real(dp), allocatable :: rate(:)
    integer, intent(in) :: maxTimestep, start
    real(dp), intent(in) :: timestep, input
    character(*), intent(in) :: time_mode
    real(dp), intent(inout) :: laplacian(:,:)
    real(dp), intent(inout) :: vortex(:)
    
    n1 = size(laplacian, dim=1)
    select case (time_mode)
      case default !! Unknown Mode. End program.
        write(*,*) "Error: Not known modus! Stop program."
        stop
        
      case ("discrete") !! discrete-time Random Walk mode
        !! Transform laplacian in transition matrix
        laplacian = laplacian * timestep
        do ii = 1, n1
          laplacian(ii,ii) = laplacian(ii,ii) + 1
        end do
        laplacian(start, start) = laplacian(start, start) + input*timestep
        do ii = 1, maxTimestep
          vortex = matmul(laplacian,vortex)
        end do
        
      case ("continous") !! continous-time Random Walk mode
      allocate(rate(n1))
      do ii = 1, maxTimestep
        rate = matmul(laplacian,vortex)
        rate(start) = rate(start) + input
        vortex = vortex + rate*timestep
      end do
    end select
     
     
    

  end subroutine CRWalk

end module randomWalk
