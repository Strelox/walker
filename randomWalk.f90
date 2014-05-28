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
  !! \param io_vertices Input/Output vertices 
  !! \param io_rates    Input/Output rates
  !! \param timestep    The size of a timestep.
  !! \param maxTimestep The maximum amount of timesteps that will be simulated.
  !! \param time_mode   The time mode of the simulation (either discrete or continous).
  !!
  subroutine CRWalk(laplacian, vortex, io_vertices, io_rates, timestep, maxTimestep, time_mode)

    !! Declarations
    integer, intent(in) :: maxTimestep, io_vertices(:)
    real(dp), intent(in) :: timestep, io_rates(:)
    character(*), intent(in) :: time_mode
    real(dp), intent(inout) :: laplacian(:,:)
    real(dp), intent(inout) :: vortex(:,:)

    integer :: ii, n1
    real(dp), allocatable :: rate(:)

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
        !! Add Input/Output rates
        do ii = 1, size(io_vertices)
          laplacian(io_vertices(ii), io_vertices(ii)) = laplacian(io_vertices(ii), io_vertices(ii)) + io_rates(ii)*timestep
        end do
        !! Lets walk
        do ii = 1, maxTimestep
          vortex(:,ii+1) = matmul(laplacian,vortex(:,ii))
        end do

      case ("continous") !! continous-time Random Walk mode
      allocate(rate(n1))
      !! Add Input/Output rates
      do ii=1, size(io_vertices)
        laplacian(io_vertices(ii), io_vertices(ii)) = laplacian(io_vertices(ii), io_vertices(ii)) + io_rates(ii)
      end do
      !! Lets walk
      do ii = 1, maxTimestep
        rate = matmul(laplacian,vortex(:,ii))
        vortex(:,ii+1) = vortex(:,ii) + rate*timestep
      end do
    end select

  end subroutine CRWalk

end module randomWalk