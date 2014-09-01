!     Copyright 2014 Frank Stuckenberg
!
!     This file is part of walker.
! 
!     walker is free software: you can redistribute it and/or modify
!     it under the terms of the GNU Affero General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     walker is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU Affero General Public License for more details.
! 
!     You should have received a copy of the GNU Affero General Public License
!     along with walker.  If not, see <http://www.gnu.org/licenses/>.
! 
!     Diese Datei ist Teil von walker.
! 
!     walker ist Freie Software: Sie können es unter den Bedingungen
!     der GNU Affero General Public License, wie von der Free Software Foundation,
!     Version 3 der Lizenz oder (nach Ihrer Wahl) jeder späteren
!     veröffentlichten Version, weiterverbreiten und/oder modifizieren.
! 
!     walker wird in der Hoffnung, dass es nützlich sein wird, aber
!     OHNE JEDE GEWÄHELEISTUNG, bereitgestellt; sogar ohne die implizite
!     Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
!     Siehe die GNU Affero General Public License für weitere Details.
! 
!     Sie sollten eine Kopie der GNU Affero General Public License zusammen mit diesem
!     Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.

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
    character(*), parameter :: fvertex = "vertex.dat"
    character(*), parameter :: fdistance = "distance.dat"
    character(*), parameter :: ftotal_vertex = "total_vertex.dat"
    character(*), parameter :: ftotal_distance = "total_distance.dat"
    character(*), parameter :: fprob = "prob_distance.dat"
    character(*), parameter :: fentropy = "entropy.dat"
    character(*), parameter :: fparticle = "particle.dat"
    character(*), parameter :: frel = "rel_distance.dat"
    character(*), parameter :: fcurrent = "current.dat"
    character(*), parameter :: frecurrent = "recurrent.dat"
    character(*), parameter :: fend_vertex = "end_vertex.dat"
    character(*), parameter :: fend_entropy = "end_entropy.dat"
    character(*), parameter :: fend_particle = "end_particle.dat"
    character(*), parameter :: fend_current = "end_current.dat"
    character(*), parameter :: fend_recurrent = "end_recurrent.dat"
    character(*), parameter :: fhamiltonian = "hamiltonian.inp"
    character(*), parameter :: fdensity = "density_matrix.dat"
    character(*), parameter :: feigenval = "eigenval.dat"
    character(*), parameter :: feigenvec = "eigenvec.dat"
    character(*), parameter :: fenvironment = "environment.inp"
    character(*), parameter :: flog = "error.log"
    character(*), parameter :: ftrace = "trace.dat"
    character(*), parameter :: fend_density = "last_density_matrix.inp"
    
    !! Declarations
    integer :: ii, jj, kk, ll, maxTimestep, nn, counter
    integer, allocatable :: distance(:)
    real(dp) :: timestep, io_rates(2), entropy, p_env, dec_time
    real(dp), allocatable :: vertex(:), total_vertex(:), prob_distance(:), total_distance(:)
    real(dp), allocatable :: current(:), recurrent(:), laplacian(:,:)
    real(dp) :: test_particle(test_range), test_entropy(test_range)
    character(10) :: walk_mode, simulation_mode, write_mode
    complex(dp), allocatable :: hamiltonian(:,:), density_matrix(:,:), environment(:,:,:)
    complex(dp) :: trace
    logical :: neg_ev, complex_diag, neg_diag, not_hermitian, complex_trace, not_norm, correct
    
    !! Open error.log
    open(60, file=flog, status="replace", form="formatted", action="write")
    
    !! Read Configuration
    call read_config(fconfig, maxTimestep, walk_mode, io_rates, timestep, simulation_mode, write_mode, dec_time)
    
    !! Change io_rates per timestep to io_rate
    io_rates = io_rates/timestep
    
    select case (walk_mode)
    case default !! Unknown Mode. End program.
        write(*,*) "Error: Not known modus! Stop program."
        stop
    case ("CRW") !! Classical Random Walk mode
        !! Read Laplacian matrix
        call read_matrix_real(flaplacian, laplacian)
        nn = size(laplacian, dim=1)
        
        !! Allocation
        allocate(vertex(nn))
        allocate(total_vertex(nn))
        
        allocate(current(nn))
        allocate(recurrent(nn))
        
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
            vertex = 0
            vertex(1) = 1
            if (write_mode == "all") then
                call write_vec_real(fvertex, vertex, horizontal=.true.)
            end if
            
            !! Initialize total vertex
            total_vertex = vertex
            
            !! Initialize entropy
            entropy = calc_entropy(vertex)
            if (write_mode == "all") then
                call write_real(fentropy, entropy)
            end if
            
            !! calculate probability dependent on distance from start
            do kk = 0, maxval(distance)
                do jj = 1, size(distance)
                    if (distance(jj) == kk) then
                        prob_distance(kk+1) = prob_distance(kk+1) + vertex(jj)
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
                call write_real(fparticle, sum(vertex))
            
                !! Create current.dat
                open(40, file="current.dat", status="replace", form="formatted", action="write")
                close(40)
                !! Create recurrent.dat
                open(41, file="recurrent.dat", status="replace", form="formatted", action="write")
                close(41)
            end if
            
        case ("append") !! Simulate old data further
            !! Read old vertex (particle distribution)
            open(42, file=fvertex, status="old", form="formatted", action="read", position="append")
            backspace 42
            read(42,*) (vertex(ii), ii=1, nn)
            close(42)
            
            !! Read old total_vertex
            call read_vec_real(ftotal_vertex, total_vertex, vec_size=nn)
            
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
            call calc_current(current, laplacian, vertex, io_rates, distance)
            if (write_mode == "all") then
                call write_matrix_real(fcurrent, reshape(current, [1, size(recurrent)]), state_in = "old", pos_in="append")
            end if
            
            !! Calculate recurrent
            call calc_recurrent(recurrent, laplacian, vertex, io_rates, distance)
            if (write_mode == "all") then
                call write_matrix_real(frecurrent, reshape(recurrent, [1, size(recurrent)]), state_in = "old", pos_in="append")
            end if
            
            call CRWalk(laplacian, vertex, io_rates, timestep)
            if (write_mode == "all") then
                call write_matrix_real(fvertex, reshape(vertex, [1, size(vertex)]), state_in = "old", pos_in="append")
            end if
            total_vertex = total_vertex + vertex    !! Integrate probability over time
            
            entropy = calc_entropy(vertex)    !! Calculate Entropy
            if (write_mode == "all") then
                call write_real(fentropy, entropy, state_in="old", pos_in="append")
            end if
            
            if (maxTimestep - ii < test_range) then
                test_entropy(maxTimestep - ii+1) = entropy
                test_particle(maxTimestep - ii+1) = sum(vertex)
            end if
            
            !! calculate probability dependent on distance from start
            do kk = 0, maxval(distance)
                do jj = 1, size(distance)
                    if (distance(jj) == kk) then
                        prob_distance(kk+1) = prob_distance(kk+1) + vertex(jj)
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
                call write_real(fparticle, sum(vertex), state_in="old", pos_in="append") !! Write amount of particle in file
            end if
        end do
        write(*,"(A)") "] Done."
        
        !! Write results into files
        call write_x(ftotal_vertex, total_vertex)
        call write_x(ftotal_distance, total_distance)
        call write_x(fend_vertex, vertex)
        call write_real(fend_entropy, entropy)
        call write_real(fend_particle, sum(vertex))
        call write_vec_real(fend_current, current)
        call write_vec_real(fend_recurrent, recurrent)
        if (write_mode == "end") then
            call write_vec_real(fprob, prob_distance)
        end if
        
        
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
        
    case ("QW", "QSW") !! Quantum Walk
        test_entropy = 0.0_dp
        test_particle = 0.0_dp
        
        neg_ev = .false.
        complex_diag = .false.
        neg_diag = .false.
        not_hermitian = .false.
        complex_trace = .false.
        not_norm = .false.
        correct = .false.

        !! Read hamiltonian
        call read_matrix_complex(fhamiltonian, hamiltonian)
        nn = size(hamiltonian, dim=1)
        
        allocate(current(nn))
        current = 0.0_dp
        
        
        !! Test hermitian of hamiltonian
            not_hermitian = .not. (is_hermitian(hamiltonian))
            if (not_hermitian .eqv. .true.) then
                write(60, "(A)") "Warning: not_hermitian hamiltonian"
            end if

        !! Read initial density matrix
        call read_matrix_complex("density_matrix.inp", density_matrix)
        
         !! Initial Particle Numer
        trace = trace_complex(density_matrix)
      
        if (aimag(trace) /= 0.0_dp) then
            write(60, "(A,I0)") "Warning: Trace complex in ", 0
            complex_trace = .true.
        end if
        
        !! Creates Vertex
        allocate(vertex(nn))
        do ii=1,nn
            vertex(ii) = real(density_matrix(ii,ii), dp)
        end do
        
        !! Test hermitian of density matrix
            not_hermitian = .not. (is_hermitian(density_matrix))
            if (not_hermitian .eqv. .true.) then
                write(60, "(A,I0)") "Warning: not_hermitian in ", 0
            end if
        
        !! Read Environment term if QSW mode
        if (walk_mode == "QSW") then
            call read_array_complex(fenvironment, environment)
            if  ( (size(environment, dim=1) /= nn ) .or. (size(environment, dim=2) /= nn) ) then
                write(*,*) "Error: Wrong size in environment. Have to be n x n x ?"
                stop
            end if
        end if
        
        hamiltonian = man_exp(hamiltonian, cmplx(0.0_dp, -timestep, dp))
!         hamiltonian = special_2(cmplx(-timestep, 0.0_dp, dp))
        p_env = 1 - exp(-1/dec_time)
        
        selectcase (simulation_mode)
        case ("standard")

             !! Creates new eigenvalue file
            open(43, file=feigenval, status="replace", form="formatted", action="write")
            close(43)
            
            !! Initial Entropy calculation
            entropy = vonNeumann_entropy(density_matrix)
            call write_vec_real(fcurrent, current, state_in="replace", horizontal=.true.)
            call write_real(fparticle, real(trace))
            call write_complex(ftrace, trace)
            call write_vec_real(fvertex, vertex, horizontal=.true.)
            call write_real(fentropy, entropy, state_in="replace")
            call write_vec_complex(fdensity, reshape(density_matrix, [(nn**2)]), state_in = "replace", horizontal=.true.)
       case ("append")
            call read_matrix_complex(fend_density, density_matrix)
        end select
        
        !! Initialize Progress bar    
        write(*,"(A)") "____________-________-_________-"
        write(*,"(A)", advance="no") "Start QW: [" !! start progress bar
        counter = 1
        
        !! Do Quantum (Stochastic) Walk
        do ii = 1, maxTimestep
            !! progress bar
            if (ii == (counter*maxTimestep/20)) then
                write(*,"(A1)", advance="no") "#"
                counter = counter + 1
            end if
            
            !! 1-timestep of Quantum Walk
            if (walk_mode == "QW") then
                call QWalk(hamiltonian, density_matrix, timestep, current, p_env)
            else 
                call QWalk(hamiltonian, density_matrix, timestep, current, p_env, environment)
            end if
            
            call write_vec_real(fcurrent, current, state_in ="old", pos_in="append", horizontal=.true.)
            
            !! Test negative diag and complex diag elements
            do jj = 1, nn
                if (real(density_matrix(jj,jj)) < 0.0_dp) then
                    correct = .true.
                    density_matrix(jj,jj) = (0.0_dp, 0.0_dp)
                    write(60, "(A,I0)") "Warning: neg_diag in ", ii
                    neg_diag = .true.
                end if
                if (aimag(density_matrix(jj,jj)) /= 0.0_dp) then
                    density_matrix(jj,jj) = cmplx(real(density_matrix(jj,jj)), 0.0_dp, dp)
                    write(60, "(A,I0)") "Warning: Diag complex in ", ii
                    complex_diag = .true.
                end if
            end do
            
            !! Adding Particles
            density_matrix(1,1) = density_matrix(1,1) + cmplx(io_rates(1)*timestep, 0.0_dp, dp)
            
            if (real(density_matrix(nn,nn) + io_rates(2)*timestep) >= 0.0_dp) then
                density_matrix(nn,nn) = density_matrix(nn,nn) + cmplx(io_rates(2)*timestep, 0.0_dp, dp)
            end if
            
            !! Writes Particle Number
            trace = trace_complex(density_matrix)
            if (aimag(trace) /= 0.0_dp) then
                write(60, "(A,I0)") "Warning: Trace complex in ", ii
                complex_trace = .true.
            end if
            call write_real(fparticle, real(trace), state_in="old", pos_in="append")
            call write_complex(ftrace, trace, state_in="old", pos_in="append")
            
            call write_vec_complex(fdensity, reshape(density_matrix, [(nn**2)]), state_in="old", &
                                   &pos_in = "append", horizontal=.true.)
           
            
            !! Test hermitian of density matrix
            not_hermitian = .not. (is_hermitian(density_matrix))
            if (not_hermitian .eqv. .true.) then
                write(60, "(A,I0)") "Warning: not_hermitian in ", ii
            end if
            
             !! Calculate and write vertex
            do kk=1,nn
                vertex(kk) = real(density_matrix(kk,kk), dp)
            end do
            call write_vec_real(fvertex, vertex, state_in="old", pos_in="append", horizontal=.true.)
            
            entropy = vonNeumann_entropy(density_matrix, neg_ev)

            if (neg_ev .eqv. .true.) then
                write(60, "(A,I0)") "Warning: neg_ev in ", ii
            end if
            call write_real(fentropy, entropy, state_in="old", pos_in="append")
            
            if (maxTimestep - ii < test_range) then
                test_entropy(maxTimestep - ii+1) = entropy
                test_particle(maxTimestep - ii+1) = sum(vertex)
            end if
            
        end do !! End Quantum Walk
        write(*,"(A)") "] Done."
        
        !! Save the last state of the density matrix
        call write_matrix_complex(fend_density, density_matrix)
        
        if (complex_diag .eqv. .true.) then 
            write(*,*) "Warning: Some probabilities are complex."
        end if
        if (neg_diag .eqv. .true.) then
            write(*,*) "Warning: Some probabilities are negative."
        end if
        if (neg_ev .eqv. .true.) then
            write(*,*) "Warning: Some eigenvalues are negative."
        end if
        if (not_hermitian .eqv. .true.) then
            write(*,*) "Warning: Some density matrices are non hermitian."
        end if
        if (complex_trace .eqv. .true.) then
            write(*,*) "Warning: Some traces are complex."
        end if
        if (not_norm .eqv. .true.) then
            write(*,*) "Warning: Some traces are not one."
        end if
        
        !! Look for reaching steady state
        if (maxTimestep >= test_range) then
            if (any(abs(test_entropy - test_entropy(1)) >= err)) then
                write(*,*) "Entropy did not reach steady state!"
            end if
            if (any(abs(test_particle - test_particle(1)) >= err)) then
                write(*,*) "Particle did not reach steady state!"
            end if
        else
            write(*,*) "Could not verify steady state due short simulation time."
        end if
        
        
    end select

    !! Close error.log
    close(60)
end program walker
