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
    subroutine QWalk(hamiltonian, density_matrix, timestep, correction)
    
        real(dp), intent(in) :: timestep, correction
        complex(dp), intent(in) :: hamiltonian(:,:)
        complex(dp), intent(inout) :: density_matrix(:,:)
        integer :: nn, ii
        complex(dp), allocatable :: repair(:,:)

        nn = size(density_matrix, dim=1)

        density_matrix = density_matrix + ((0.0_dp, 1.0_dp)*matmul(density_matrix, hamiltonian)&
        & - (0.0_dp, 1.0_dp)*matmul(hamiltonian, density_matrix))*timestep 

        allocate(repair(nn,nn))
        repair = (0.0_dp, 0.0_dp)
        do ii = 1, nn
            repair(ii,ii) = 1/cmplx(nn)
        end do

        density_matrix = (1-correction)*density_matrix + correction*repair
        
    end subroutine QWalk
    
    !> Routine for calculating the von Neumann entropy of a density matrix
    !!
    !! /param density_matrix    density matrix of the system
    !! /param err_range         maximum precision (e.g. timestep scale)
    !!
    real(dp) function vonNeumann_entropy(density_matrix, err_range)
        complex(dp), intent(in) :: density_matrix(:,:)
        complex(dp), allocatable :: eigenvec(:,:), diagonal(:,:)
        real(dp), allocatable :: eigenval(:)
        integer :: ii, nn
        real(dp), optional :: err_range
        character(*), parameter :: feigenval = "eigenval.dat"
        
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
end module quantumWalk