!> Walker
!! Let a random Walker walk through a network.
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
  integer :: ii, maxTimestep, n1, nWalker, start
  real(dp) :: timestep
  real(dp), allocatable :: vortex(:)
  real(dp), allocatable :: laplacian(:,:)
  character(*), parameter :: file = "laplacian.inp"
  character(10) :: mode
  
  call init_random_seed()
  call read_config("walk.cfg", nWalker, maxTimestep, mode, start, timestep)
  
  select case (mode)
  case default !! Unknown Mode. End program.
    write(*,*) "Error: Not known modus! Stop program."
    stop
  case ("CRW") !! Classical Random Walk mode
     
     call read_matrix_real(file, laplacian)
     n1 = size(laplacian, dim=1)
     allocate(vortex(n1))
     vortex = 0
     vortex(start) = 1
     laplacian = -laplacian * timestep
     do ii = 1, n1
      laplacian(ii,ii) = laplacian(ii,ii) + 1
     end do
     call CRWalk(laplacian,vortex,maxTimestep) 
     call write_x("vortex.dat", vortex)

  case ("QW") !! Quantum Walk
    write(*,*) "Sorry. Not yet implemented. Stop program."
    stop
  end select

end program

