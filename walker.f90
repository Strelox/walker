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
    use networks
    use randomWalk
    implicit none

    !! Declarations
    integer :: ii, jj, kk, ll, maxTimestep, n1
    integer, allocatable :: distance(:)
    real(dp) :: timestep, io_rates(2)
    real(dp), allocatable :: vortex(:), total_vortex(:), prob_distance(:), total_distance(:) 
    real(dp), allocatable :: laplacian(:,:), entropy(:)
    character(10) :: walk_mode, time_mode
    character(*), parameter :: flaplacian = "laplacian.inp"
    character(*), parameter :: fconfig = "walk.cfg"
    character(*), parameter :: fvortex = "vortex.dat"
    character(*), parameter :: fdistance = "distance.dat"
    character(*), parameter :: ftotal_vortex = "total_vortex.dat"
    character(*), parameter :: ftotal_distance = "total_distance.dat"
    character(*), parameter :: fprob = "prob_distance.dat"
    character(*), parameter :: fentropy = "entropy.dat"
    character(*), parameter :: fparticle = "particle.dat"
    character(*), parameter :: frel = "rel_distance.dat"
 
    !! Read Configuration
    call read_config(fconfig, maxTimestep, walk_mode, time_mode, io_rates, timestep)
    
    !! Change io_rates per timestep to io_rate
    io_rates = io_rates/timestep
    
    select case (walk_mode)
    case default !! Unknown Mode. End program.
        write(*,*) "Error: Not known modus! Stop program."
        stop
    case ("CRW") !! Classical Random Walk mode
        !! Read Laplacian matrix
        call read_matrix_real(flaplacian, laplacian)
        n1 = size(laplacian, dim=1)

        !! Initialize start distribution
        allocate(vortex(n1))
        vortex = 0
        vortex(1) = 1
        call write_matrix_real(fvortex, reshape(vortex, [1, n1]))
        
        !! Initialize total vortex
        allocate(total_vortex(n1))
        total_vortex = vortex
        
        
        !! Initialize entropy
        allocate(entropy(maxTimestep+1))
        entropy = 0
        entropy(1) = calc_entropy(vortex)
        
        !! Calculating distance of vertices to start
        call calc_distance(laplacian, distance)
        call write_x_int(fdistance, distance)
        allocate(prob_distance((maxval(distance)+1)))
        prob_distance = 0
        
        !! calculate probability dependent on distance from start
        do kk = 0, maxval(distance)
            do jj = 1, size(distance)
                if (distance(jj) == kk) then
                    prob_distance(kk+1) = prob_distance(kk+1) + vortex(jj)
                end if
            end do
        end do
        call write_matrix_real(fprob, reshape(prob_distance, [1, size(prob_distance)]))
        do ii = 1, size(prob_distance)
            prob_distance(ii) = prob_distance(ii)/count(distance == distance(ii))
        end do
        call write_matrix_real(frel, reshape(prob_distance, [1, size(prob_distance)]))
        prob_distance = 0
        
        !! Initialize total distance
        allocate(total_distance(n1))
        total_distance = vortex
        
        !! Write amount of particle in file
        open(22, file=fparticle, status="replace", form="formatted", action="write")
        write(22,"(ES15.6)") sum(vortex)
        
        !! Create current.dat
        open(25, file="current.dat", status="replace", form="formatted", action="write")
        close(25)
        !! Create recurrent.dat
        open(26, file="recurrent.dat", status="replace", form="formatted", action="write")
        close(26)
        
        !! do classical random walk
        do ii = 1, maxTimestep
            call calc_current(laplacian, vortex, io_rates, distance)
            call CRWalk(laplacian, vortex, io_rates, timestep, time_mode) 
            call write_matrix_real(fvortex, reshape(vortex, [1, size(vortex)]), state_in = "old", pos_in="append")
            total_vortex = total_vortex + vortex    !! Integrate probability over time
            
            entropy(ii+1) = calc_entropy(vortex)    !! Calculate Entropy
            
            !! calculate probability dependent on distance from start
            do kk = 0, maxval(distance)
                do jj = 1, size(distance)
                    if (distance(jj) == kk) then
                        prob_distance(kk+1) = prob_distance(kk+1) + vortex(jj)
                    end if
                end do
            end do
            call write_matrix_real(fprob, reshape(prob_distance, [1, size(prob_distance)]), state_in = "old", pos_in="append")
            total_distance = total_distance + prob_distance
            do ll = 1, size(prob_distance)
                prob_distance(ll) = prob_distance(ll)/count(distance == distance(ll))
            end do
            call write_matrix_real(frel, reshape(prob_distance, [1, size(prob_distance)]), state_in = "old", pos_in="append")
            prob_distance = 0
            
            write(22,"(ES15.6)") sum(vortex) !! Write amount of particle in file
        end do
        close(22)
        
        !! Write results into files
        call write_x(fentropy, entropy)
        call write_x(ftotal_vortex, total_vortex)
        call write_x(ftotal_distance, total_distance)
        
    case ("QW") !! Quantum Walk
        write(*,*) "Sorry. Not yet implemented. Stop program."
        stop
    end select

end program walker