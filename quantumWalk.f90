!> Provides routines for the quantum walk
!!
module quantumWalk
    use accuracy
    use io
    use random
    use f95_lapack
    implicit none
    
contains
    
    !> Calculates the next timestep of the density matrix
    !!
    !! /param hamiltonian       Hamiltonian matrix of the system
    !! /param density_matrix    density matrix of the system
    !! /param timestep          simulation timestep
    !!
    subroutine QWalk(hamiltonian, density_matrix, timestep, io_rates, correction, environment)
        
        !! Declarations
        real(dp), intent(in) :: timestep, correction, io_rates(2)
        complex(dp), intent(in) :: hamiltonian(:,:)
        complex(dp), intent(inout) :: density_matrix(:,:)
        integer :: nn, n3, ii
        complex(dp), allocatable :: repair(:,:), temp(:,:), temp2(:,:)
        complex(dp), optional, intent(in) :: environment(:,:,:)
        
        nn = size(density_matrix, dim=1)
        n3 = size(environment, dim=3)
        
        allocate(temp(nn,nn))
        allocate(temp2(nn,nn))
        temp = vonNeumann_eq(density_matrix, hamiltonian)
        temp2 = vonNeumann_eq(temp, hamiltonian)
        
        if (present(environment) .eqv. .true.) then
            density_matrix = density_matrix + vonNeumann_eq(density_matrix, hamiltonian)*timestep &
                                          & + environment_term(density_matrix, environment)*timestep
        else
            density_matrix = density_matrix + vonNeumann_eq(density_matrix, hamiltonian)*timestep &
                            & + 0.5*temp2*(timestep**2) + (1/6)*vonNeumann_eq(temp2,hamiltonian)&
                            & *timestep**3
        end if
        
        allocate(repair(nn,nn))
        repair = (0.0_dp, 0.0_dp)
        do ii = 1, nn
            repair(ii,ii) = 1/cmplx(nn)
        end do

        density_matrix = (1-correction)*density_matrix + correction*repair
        
        !! Adding Particles
        density_matrix(1,1) = density_matrix(1,1) + io_rates(1)*timestep
        density_matrix(nn,nn) = density_matrix(nn,nn) + io_rates(2)*timestep
        if (real(density_matrix(nn,nn)) < 0.0_dp) then
            density_matrix(nn,nn) = (0.0_dp, 0.0_dp)
        end if
    end subroutine QWalk
    
    function vonNeumann_eq(density_matrix, hamiltonian) result(vonNeumann)
        
        complex(dp), intent(in) :: hamiltonian(:,:)
        complex(dp), intent(inout) :: density_matrix(:,:)
        integer :: nn
        complex(dp), allocatable :: vonNeumann(:,:)
        nn = size(density_matrix, dim=1)
        allocate(vonNeumann(nn,nn))
        
        vonNeumann = ((0.0_dp, 1.0_dp)*matmul(density_matrix, hamiltonian) &
                       & - (0.0_dp, 1.0_dp)*matmul(hamiltonian, density_matrix))
                       
   
    end function vonNeumann_eq
   
    !> Routine for calculating the von Neumann entropy of a density matrix
    !!
    !! /param density_matrix    density matrix of the system
    !! /param err_range         maximum precision (e.g. timestep scale)
    !!
    real(dp) function vonNeumann_entropy(density_matrix, err_range, neg_ev)
        complex(dp), intent(in) :: density_matrix(:,:)
        complex(dp), allocatable :: eigenvec(:,:), diagonal(:,:)
        real(dp), allocatable :: eigenval(:)
        integer :: ii, nn
        real(dp), optional :: err_range
        character(*), parameter :: feigenval = "eigenval.dat"
        logical, optional, intent(inout) :: neg_ev
        
        nn = size(density_matrix, dim=1)

        allocate(eigenvec(nn, nn))
        allocate(eigenval(nn))
        allocate(diagonal(nn, nn))
        
        !! Calculation of the eigenvectors and eigenvalues
        eigenvec = density_matrix
        call la_heev(eigenvec, eigenval, JOBZ= "V")
        
        !! Truncate numerical errors
        if (present(err_range) .eqv. .true.) then
            eigenval = dint(eigenval/err_range)*err_range
        end if
        
        call write_vec_real(feigenval, eigenval, state_in="old", pos_in="append", horizontal=.true.)
        
        !! Creates the diagonal matrix of the density_matrix
        diagonal = (0.0_dp, 0.0_dp)
        do ii  = 1, nn
            if (eigenval(ii) > 0.0_dp) then
                diagonal(ii, ii) = complex(log(eigenval(ii)), dp)
            else
                if (present(neg_ev) .eqv. .true.) then
                    neg_ev = .true.
                end if
            end if
        end do
        
        !! Calculates the matrix multiplication inside the trace in the von Neumann entropy
        eigenvec = matmul(density_matrix, matmul(eigenvec, matmul(diagonal, transpose(conjg(eigenvec)))))
        
        !! Calculates the negative trace
        vonNeumann_entropy = 0
        do ii = 1, nn  
                vonNeumann_entropy = vonNeumann_entropy - real(eigenvec(ii, ii), dp)
        end do
        
    end function vonNeumann_entropy
    
    !! Calculates the trace of a complex matrix
    !!
    !! /param matrix    matrix whose trace will be calculated
    !! 
    !! The Output will also be a complex data type
    !!
    complex(dp) function trace_complex(matrix)
        complex(dp), intent(in) :: matrix(:,:)
        integer :: ii
        
        trace_complex = 0
        do ii = 1, size(matrix, dim=1)
            trace_complex = trace_complex + matrix(ii, ii)
        end do
        
    end function trace_complex
    
    function environment_term(density_matrix, environment) result(term)
        complex(dp), allocatable :: term(:,:) 
        complex(dp) :: environment(:,:,:), density_matrix(:,:)
        integer :: ii, nn, n_env
        
        nn = size(density_matrix, dim=1)
        n_env = size(environment, dim=3)
        
        if ( (size(density_matrix, dim=2) /= nn) .or. (size(environment, dim=1) /= nn) .or. (size(environment, dim=2) /= nn) ) then
            write(*,*) "Error: density matrix and environment superoperator do not match in lindblad_term. Stop."
            stop
        end if
        
        allocate(term(nn,nn))
        term = 0
        
        do ii = 1, n_env
            term =  term + 2*matmul(environment(:,:,ii), matmul(density_matrix, transpose(conjg(environment(:,:,ii))))) &
                    & - matmul(environment(:,:,ii), matmul(transpose(conjg(environment(:,:,ii))), density_matrix)) &
                    & - matmul(transpose(conjg(environment(:,:,ii))), matmul(environment(:,:,ii), density_matrix)) 
        end do
        
     end function environment_term  
     
     logical function is_hermitian(matrix)
        complex(dp), intent(in) :: matrix(:,:)
        
        if ( any( (matrix - transpose(conjg(matrix))) /= 0.0) ) then
            is_hermitian = .true.
        else
            is_hermitian = .false.
        end if
        
    end function is_hermitian
   
end module quantumWalk