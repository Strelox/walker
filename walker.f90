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
    use quantumWalk
    use f95_lapack
    implicit none

    !! Parameters
    integer, parameter :: test_range = 100
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
    character(*), parameter :: fcurrent = "current.dat"
    character(*), parameter :: frecurrent = "recurrent.dat"
    character(*), parameter :: fend_vortex = "end_vortex.dat"
    character(*), parameter :: fend_entropy = "end_entropy.dat"
    character(*), parameter :: fend_particle = "end_particle.dat"
    character(*), parameter :: fend_current = "end_current.dat"
    character(*), parameter :: fend_recurrent = "end_recurrent.dat"
    character(*), parameter :: fhamiltonian = "hamiltonian.inp"
    character(*), parameter :: fdensity = "density_matrix.dat"
    character(*), parameter :: ftrace = "trace.dat"
    character(*), parameter :: feigenval = "eigenval.dat"
    character(*), parameter :: feigenvec = "eigenvec.dat"
    
    !! Declarations
    integer :: ii, jj, kk, ll, maxTimestep, n1, counter
    integer, allocatable :: distance(:)
    real(dp) :: timestep, io_rates(2), entropy, correction
    real(dp), allocatable :: vortex(:), total_vortex(:), prob_distance(:), total_distance(:)
    real(dp), allocatable :: current(:), recurrent(:), laplacian(:,:)
    real(dp) :: test_particle(test_range), test_entropy(test_range)
    character(10) :: walk_mode, time_mode, simulation_mode, write_mode
    complex(dp), allocatable :: hamiltonian(:,:), density_matrix(:,:)
    complex(dp) :: trace
    
    !! Read Configuration
    call read_config(fconfig, maxTimestep, walk_mode, time_mode, io_rates, timestep, simulation_mode, write_mode, correction)
    
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
        
        allocate(current(n1))
        allocate(recurrent(n1))
        
        test_entropy = 0
        test_particle = 0
        
        !! Calculating distance of vertices to start
        call calc_distance(laplacian, distance)     
        call write_x_int(fdistance, distance)
        allocate(prob_distance((maxval(distance)+1)))
        allocate(total_distance((maxval(distance)+1)))
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
            if (write_mode == "all") then
                call write_matrix_real(fvortex, reshape(vortex, [1, n1]))
            end if
            
            !! Initialize total vortex
            total_vortex = vortex
            
            !! Initialize entropy
            entropy = calc_entropy(vortex)
            if (write_mode == "all") then
                call write_real(fentropy, entropy)
            end if
            
            !! calculate probability dependent on distance from start
            do kk = 0, maxval(distance)
                do jj = 1, size(distance)
                    if (distance(jj) == kk) then
                        prob_distance(kk+1) = prob_distance(kk+1) + vortex(jj)
                    end if
                end do
            end do
            if (write_mode == "all") then
                call write_matrix_real(fprob, reshape(prob_distance, [1, size(prob_distance)]))
            end if
            do ii = 1, size(prob_distance)
                prob_distance(ii) = prob_distance(ii)/count(distance == distance(ii))
            end do
            if (write_mode == "all") then
                call write_matrix_real(frel, reshape(prob_distance, [1, size(prob_distance)]))
            end if
            prob_distance = 0
            
            !! Write amount of particle in file
            if (write_mode == "all") then
                call write_real(fparticle, sum(vortex))
            
                !! Create current.dat
                open(40, file="current.dat", status="replace", form="formatted", action="write")
                close(40)
                !! Create recurrent.dat
                open(41, file="recurrent.dat", status="replace", form="formatted", action="write")
                close(41)
            end if
            
        case ("append") !! Simulate old data further
            !! Read old vortex (particle distribution)
            open(42, file=fvortex, status="old", form="formatted", action="read", position="append")
            backspace 42
            read(42,*) (vortex(ii), ii=1, n1)
            close(42)
            
            !! Read old total_vortex
            call read_vec_real(ftotal_vortex, total_vortex, vec_size=n1)
            
            !! Read old total_distance
            call read_vec_real(ftotal_distance, total_distance, vec_size=size(total_distance))
             
        end select
        
        !! do classical random walk
        write(*,"(A)") "____________-________-_________-"
        write(*,"(A)", advance="no") "Start CRW: [" !! start progress bar
        counter = 1
        do ii = 1, maxTimestep
            !! progress bar
            if (ii == (counter*maxTimestep/20)) then
                write(*,"(A1)", advance="no") "#"
                counter = counter + 1
            end if
            
            !! Calculate current
            call calc_current(current, laplacian, vortex, io_rates, distance)
            if (write_mode == "all") then
                call write_matrix_real(fcurrent, reshape(current, [1, size(recurrent)]), state_in = "old", pos_in="append")
            end if
            
            !! Calculate recurrent
            call calc_recurrent(recurrent, laplacian, vortex, io_rates, distance)
            if (write_mode == "all") then
                call write_matrix_real(frecurrent, reshape(recurrent, [1, size(recurrent)]), state_in = "old", pos_in="append")
            end if
            
            call CRWalk(laplacian, vortex, io_rates, timestep, time_mode)
            if (write_mode == "all") then
                call write_matrix_real(fvortex, reshape(vortex, [1, size(vortex)]), state_in = "old", pos_in="append")
            end if
            total_vortex = total_vortex + vortex    !! Integrate probability over time
            
            entropy = calc_entropy(vortex)    !! Calculate Entropy
            if (write_mode == "all") then
                call write_real(fentropy, entropy, state_in="old", pos_in="append")
            end if
            
            if (maxTimestep - ii < test_range) then
                test_entropy(maxTimestep - ii+1) = entropy
                test_particle(maxTimestep - ii+1) = sum(vortex)
            end if
            
            !! calculate probability dependent on distance from start
            do kk = 0, maxval(distance)
                do jj = 1, size(distance)
                    if (distance(jj) == kk) then
                        prob_distance(kk+1) = prob_distance(kk+1) + vortex(jj)
                    end if
                end do
            end do
            if (write_mode == "all") then
                call write_matrix_real(fprob, reshape(prob_distance, [1, size(prob_distance)]), state_in = "old", pos_in="append")
            end if
            total_distance = total_distance + prob_distance
            do ll = 1, size(prob_distance)
                prob_distance(ll) = prob_distance(ll)/count(distance == distance(ll))
            end do
            if (write_mode == "all") then
                call write_matrix_real(frel, reshape(prob_distance, [1, size(prob_distance)]), state_in = "old", pos_in="append")
            end if
            prob_distance = 0
            if (write_mode == "all") then
                call write_real(fparticle, sum(vortex), state_in="old", pos_in="append") !! Write amount of particle in file
            end if
        end do
        write(*,"(A)") "] Done."
        
        !! Write results into files
        call write_x(ftotal_vortex, total_vortex)
        call write_x(ftotal_distance, total_distance)
        call write_x(fend_vortex, vortex)
        call write_real(fend_entropy, entropy)
        call write_real(fend_particle, sum(vortex))
        call write_vec_real(fend_current, current)
        call write_vec_real(fend_recurrent, recurrent)
        
        !! Look for reaching steady state
        if (maxTimestep >= test_range) then
            if (any(test_entropy /= test_entropy(1))) then
                write(*,*) "Entropy did not reach steady state!"
            end if
            if (any(test_particle /= test_particle(1))) then
                write(*,*) "Particle did not reach steady state!"
            end if
        else
            write(*,*) "Could not verify steady state due short simulation time."
        end if
        
    case ("QW") !! Quantum Walkdensity_matrix
    
        !! Creates new eigenvalue file
        open(43, file=feigenval, status="replace", form="formatted", action="write")
        close(43)

        
        !! Read hamiltonian
        call read_matrix_complex(fhamiltonian, hamiltonian)
        n1 = size(hamiltonian, dim=1)
        
        allocate(density_matrix(n1, n1))
        density_matrix = (0.0_dp, 0.0_dp)
        density_matrix(1,1) = (1.0_dp, 0.0_dp)
        
        call write_complex(ftrace, trace_complex(density_matrix)) 
        !! Initialize Progress bar    
        write(*,"(A)") "____________-________-_________-"
        write(*,"(A)", advance="no") "Start QW: [" !! start progress bar
        counter = 1
        
        !! Initial Entropy calculation
        entropy = vonNeumann_entropy(density_matrix)
        call write_real(fentropy, entropy, state_in="replace")
        call write_vec_complex(fdensity, reshape(density_matrix, [(n1**2)]), state_in = "replace", horizontal=.true.)
        
        do ii = 1, maxTimestep
            !! progress bar
            if (ii == (counter*maxTimestep/20)) then
                write(*,"(A1)", advance="no") "#"
                counter = counter + 1
            end if
            
            !! 1-timestep of Quantum Walk
            call QWalk(hamiltonian, density_matrix, timestep, correction)
            call write_vec_complex(fdensity, reshape(density_matrix, [(n1**2)]), state_in="old", &
                                   &pos_in = "append", horizontal=.true.)
            trace = trace_complex(density_matrix)
            if ( (real(trace, dp) > 1.0_dp + timestep) .or. (real(trace, dp) < 1.0_dp - timestep) )then
                write(*,*) "Warning: Trace of density_matrix is not 1!"
            end if
            call write_complex(ftrace, trace_complex(density_matrix), state_in = "old", pos_in = "append") 
            entropy = vonNeumann_entropy(density_matrix, timestep)
            call write_real(fentropy, entropy, state_in="old", pos_in="append")
        end do
        write(*,"(A)") "] Done."
        
    end select

end program walker