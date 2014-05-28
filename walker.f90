!> Let a random (quantum) walker go through a network.
!!
!! \details The provided laplacian matrix of the network will be used to let
!! different types of walkers (Classical Random Walk, Quantum Walk,
!! Quantum Stochastic Walk) perform a random walk through the network. 
!!
program walker
  use accuracy
  use io
  use config
  use random
  use networks
  use randomWalk
  implicit none

  !! Declarations
  integer :: ii, jj, kk, maxTimestep, n1, start, io_n
  integer, allocatable :: distance(:), io_vertices(:)
  real(dp) :: timestep
  real(dp), allocatable :: vortex(:,:), prob_distance(:,:), io_rates(:), laplacian(:,:)
  character(10) :: walk_mode, time_mode
  character(*), parameter :: flaplacian = "laplacian.inp"
  character(*), parameter :: fconfig = "walk.cfg"
  character(*), parameter :: fvortex = "vortex.dat"
  character(*), parameter :: fdistance = "disctance.dat"
  character(*), parameter :: ftotal_vortex = "total_vortex.dat"
  character(*), parameter :: ftotal_distance = "total_distance.dat"
  character(*), parameter :: fprob = "prob_distance.dat"

  call init_random_seed()
  call read_config(fconfig, maxTimestep, walk_mode, time_mode, start, io_vertices, io_rates, timestep)

  select case (walk_mode)
  case default !! Unknown Mode. End program.
    write(*,*) "Error: Not known modus! Stop program."
    stop
  case ("CRW") !! Classical Random Walk mode
    call read_matrix_real(flaplacian, laplacian)
    n1 = size(laplacian, dim=1)
    io_n = size(io_vertices)
    if (maxval(io_vertices) > n1) then 
       write(*,*) "Error: Some start vortex are not in network! Stop program."
       stop
    end if

    !! Initialize start distribution
    allocate(vortex(n1,maxTimestep+1))
    vortex = 0
    vortex(start,1) = 1

    !! do classical random walk
    call CRWalk(laplacian, vortex, io_vertices, io_rates, timestep, maxTimestep, time_mode)
    call write_matrix_real(fvortex, vortex)

    !! Calculating distance of vertices to start
    call calc_distance(laplacian, distance, start)
    call write_x_int(fdistance, distance)

    !! calculate probability dependent on distance from start

    allocate(prob_distance((maxval(distance)+1),maxTimestep+1))
    prob_distance = 0
    do kk = 1, maxTimestep+1
      do ii = 0, maxval(distance)
        do jj = 1, size(distance)
          if (distance(jj) == ii) then
            prob_distance(ii+1,kk) = prob_distance(ii+1,kk) + vortex(jj,kk)
          end if
        end do
      end do
    end do
    call write_matrix_real(fprob, prob_distance)

    !! Integrate probability over time
    do ii = 1, n1
      vortex(ii,1) = sum(vortex(ii,:))
    end do
    do ii = 1, size(prob_distance, dim=1)
      prob_distance(ii,1) = sum(prob_distance(ii,:))
    end do
    
    call write_x(ftotal_vortex, vortex(:,1))
    call write_x(ftotal_distance, prob_distance(:,1))

  case ("QW") !! Quantum Walk
    write(*,*) "Sorry. Not yet implemented. Stop program."
    stop
  end select

end program walker