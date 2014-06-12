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
    integer :: ii, jj, kk, ll, maxTimestep, n1, counter
    integer, allocatable :: distance(:)
    real(dp) :: timestep, io_rates(2),entropy
    real(dp), allocatable :: vortex(:), total_vortex(:), prob_distance(:), total_distance(:) 
    real(dp), allocatable :: laplacian(:,:)
    character(10) :: walk_mode, time_mode, simulation_mode
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
    call read_config(fconfig, maxTimestep, walk_mode, time_mode, io_rates, timestep, simulation_mode)
    
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
        
        !! Allocation
        allocate(vortex(n1))
        allocate(total_vortex(n1))
        allocate(total_distance(n1))
        
        !! Calculating distance of vertices to start
        call calc_distance(laplacian, distance)
        call write_x_int(fdistance, distance)
        allocate(prob_distance((maxval(distance)+1)))
        prob_distance = 0
        
        select case (simulation_mode)
        case default !! Unknown simulation mode. Stop program.
            !! Initialize start distribution
            write(*,*) "Error: Not known simulation modus! Stop program."
            stop
        case ("standard") !! Start new simulation
            !! Initialize start distribution
            vortex = 0
            vortex(1) = 1
            call write_matrix_real(fvortex, reshape(vortex, [1, n1]))
            
            !! Initialize total vortex
            total_vortex = vortex
            
            !! Initialize entropy
            entropy = calc_entropy(vortex)
            call write_real(fentropy, entropy)
            
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
            
            !! Write amount of particle in file
            call write_real(fparticle, sum(vortex))
            
            !! Create current.dat
            open(40, file="current.dat", status="replace", form="formatted", action="write")
            close(40)
            !! Create recurrent.dat
            open(41, file="recurrent.dat", status="replace", form="formatted", action="write")
            close(41)
        
        case ("append") !! Simulate old data further
            !! Read old vortex (particle distribution)
            open(42, file=fvortex, status="old", form="formatted", action="read", position="append")
            backspace 42
            read(42,*) (vortex(ii), ii=1, n1)
            close(42)
            
            !! Read old total_vortex
            call read_vec_real(ftotal_vortex, total_vortex, vec_size=n1)
            
            !! Read old total_distance
            call read_vec_real(ftotal_distance, total_distance, vec_size=n1)
             
        end select
        
        !! do classical random walk
        write(*,"(A)", advance="no") "Start CRW: [" !! start progress bar
        counter = 1
        do ii = 1, maxTimestep
            !! progress bar
            if (ii == (counter*maxTimestep/10)) then
                write(*,"(A1)", advance="no") "#"
                counter = counter + 1
            end if
            
            call calc_current(laplacian, vortex, io_rates, distance)
            call CRWalk(laplacian, vortex, io_rates, timestep, time_mode) 
            call write_matrix_real(fvortex, reshape(vortex, [1, size(vortex)]), state_in = "old", pos_in="append")
            total_vortex = total_vortex + vortex    !! Integrate probability over time
            
            entropy = calc_entropy(vortex)    !! Calculate Entropy
            call write_real(fentropy, entropy, state_in="old", pos_in="append")
            
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
            
            call write_real(fparticle, sum(vortex), state_in="old", pos_in="append") !! Write amount of particle in file
        end do
        write(*,"(A)") "] Done."
        
        !! Write results into files
        call write_x(ftotal_vortex, total_vortex)
        call write_x(ftotal_distance, total_distance)
        
    case ("QW") !! Quantum Walk
        write(*,*) "Sorry. Not yet implemented. Stop program."
        stop
    end select

end program walker