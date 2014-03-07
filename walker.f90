!> Walker
!! Let a random Walker walk through a network.
!!
program walker
  use accuracy
  use io
  use random
  use networks
  use randomWalk
  implicit none

  !! Declarations
  integer :: ii, jj, kk, ll, Time, nwalks, n1, n2, n3
  integer, allocatable :: pos(:,:), counter(:)
  real(dp) :: rand, temp, total
  real(dp), allocatable :: network(:,:,:), current(:,:,:)
  character(*), parameter :: data = "network.inp"

  call init_random_seed()
  nwalks = 20
  Time = 5
  allocate(pos(Time,nwalks))
  pos = 1
  
  call read_array(data, n1, n2, n3, network)
  allocate(counter(n1))
  counter = 0
  allocate(current(n1,n2,n3))
  
  call normal(network)

  call gen_rWalk(network,pos,Time,nwalks,counter)
  
  !! Calculation of the Potential
  counter = counter/nwalks
  do ii = 1, n1
     counter(ii) = counter(ii)/sum(network(ii,:,:))
  end do

  !! Calculation of the current
  do ii = 1, n1
     do jj = 1, n2
        do kk = 1, n3
           current(ii,jj,kk) = (counter(ii) - counter(jj))*network(ii,jj,kk)
        end do
     end do
  end do
  
  call write_matrix("matrix.dat", pos)
  call write_array("current.dat", current)
  
end program walker

