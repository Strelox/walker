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

!> Provides routines for the quantum walk
!!
module quantumWalk
    use accuracy
    use io
    use f95_lapack
    implicit none
    
    
contains
   
    !> Calculates the next timestep of the density matrix
    !!
    !! /param hamiltonian       Hamiltonian matrix of the system
    !! /param density_matrix    density matrix of the system
    !! /param timestep          simulation timestep
    !!
    subroutine QWalk(dev, density_matrix, timestep, current, p_env, environment)
        
        !! Declarations
        real(dp), intent(in) :: p_env, timestep
        complex(dp), intent(in) :: dev(:,:)
        complex(dp), intent(inout) :: density_matrix(:,:)
        real(dp), intent(inout) :: current(:)
        integer :: nn, ii, jj
        complex(dp), allocatable :: identity(:,:), temp(:,:)
        complex(dp), optional, intent(in) :: environment(:,:,:)

        nn = size(density_matrix, dim=1)
        allocate(temp(nn,nn)) 
        temp = (0.0_dp, 0.0_dp)
        allocate(identity(nn,nn))
        identity = (0.0_dp, 0.0_dp)
        do ii = 1, nn
            identity(ii,ii) = (1.0_dp, 0.0_dp)
        end do
        
        !! Quantum Development
        temp = matmul(dev, matmul(density_matrix, transpose(conjg(dev))))
        

        ! Environment Development

        temp = (1-p_env)*matmul(identity, matmul(temp, identity)) &
                        &+ p_env*environment_term(temp, environment)
        
        identity = (temp - density_matrix)/timestep

        current = 0.0_dp
        do ii = 1, nn-1
            do jj = 1, ii
                current(ii) = current(ii) - real(identity(jj,jj)) - aimag(identity(jj,jj))
            end do
        end do
        
!         current(1) = -0.5_dp*real(identity(1,1))
!         current(2) = -0.5_dp*real(identity(1,1))
!         current(3) = -0.5_dp*real(identity(1,1)) - aimag(identity(2,2))
!         current(4) = -0.5_dp*real(identity(1,1)) - aimag(identity(3,3))
!         current(5) = -0.5_dp*real(identity(1,1)) - aimag(identity(2,2)) - real(identity(4,4)) - aimag(identity(4,4))
!         current(6) = -0.5_dp*real(identity(1,1)) - aimag(identity(3,3)) - real(identity(5,5)) - aimag(identity(5,5))
    
!         current(1) = -1.0_dp/3*real(identity(1,1))
!         current(2) = -1.0_dp/3*real(identity(1,1))
!         current(3) = -1.0_dp/3*real(identity(1,1))
!         current(4) = -1.0_dp/3*real(identity(1,1)) - aimag(identity(2,2))
!         current(5) = -1.0_dp/3*real(identity(1,1)) - aimag(identity(3,3))
!         current(6) = -1.0_dp/3*real(identity(1,1)) - aimag(identity(4,4))
!         current(7) = -1.0_dp/3*real(identity(1,1)) - aimag(identity(2,2)) - real(identity(5,5)) - aimag(identity(5,5))
!         current(8) = -1.0_dp/3*real(identity(1,1)) - aimag(identity(3,3)) - real(identity(6,6)) - aimag(identity(6,6))
!         current(9) = -1.0_dp/3*real(identity(1,1)) - aimag(identity(4,4)) - real(identity(7,7)) - aimag(identity(7,7))
        
        density_matrix = temp
    end subroutine QWalk
    
    function vonNeumann_eq(density_matrix, hamiltonian) result(vonNeumann)
        
        complex(dp), intent(in) :: hamiltonian(:,:)
        complex(dp), intent(inout) :: density_matrix(:,:)
        integer :: nn
        complex(dp), allocatable :: vonNeumann(:,:)
        nn = size(density_matrix, dim=1)
        allocate(vonNeumann(nn,nn))
        
        vonNeumann = (0.0_dp, 1.0_dp)*matmul(density_matrix, hamiltonian) &
                       & - (0.0_dp, 1.0_dp)*matmul(hamiltonian, density_matrix)
                       
   
    end function vonNeumann_eq
   
    !> Routine for calculating the von Neumann entropy of a density matrix
    !!
    !! /param density_matrix    density matrix of the system
    !! /param err_range         maximum precision (e.g. timestep scale)
    !!
    real(dp) function vonNeumann_entropy(density_matrix, neg_ev)
        integer, parameter :: ap = selected_real_kind(15, 300)
        
        complex(dp), intent(in) :: density_matrix(:,:)
        complex(dp), allocatable :: eigenvec(:,:), diagonal(:,:)
        real(dp), allocatable :: eigenval(:)
        integer :: ii, nn, info
        character(*), parameter :: feigenval = "eigenval.dat"
        logical, optional, intent(inout) :: neg_ev
        complex(ap), allocatable ::  dev(:,:), left(:,:), right(:,:)
        real(ap), allocatable ::  singular(:)
        nn = size(density_matrix, dim=1)
        
        allocate(left(nn,nn))
        allocate(right(nn,nn))
        allocate(singular(nn))
        allocate(dev(nn,nn))
        
        allocate(eigenvec(nn, nn))
        allocate(eigenval(nn))
        allocate(diagonal(nn, nn))
        
        if (present(neg_ev)) then
            neg_ev = .false.
        end if
        
        !! Calculation of the eigenvectors and eigenvalues
        eigenvec = density_matrix
        eigenval = 0.0_dp
        dev = cmplx( (density_matrix/trace_complex(density_matrix)), kind=ap)
        
        call zgesvd_f95(A = dev, S= singular, U=left, VT=right, INFO=info)

        
        if (info /= 0) then
            vonNeumann_entropy = 0
            return
        end if
        
        call write_vec_real(feigenval, eigenval, state_in="old", pos_in="append", horizontal=.true.)

        !! Creates the diagonal matrix of the density_matrix
        diagonal = (0.0_dp, 0.0_dp)
        do ii  = 1, nn
            if ((present(neg_ev) .eqv. .true.) .and. (eigenval(ii) < 0.0_dp)) then
                    neg_ev = .true.
                end if
            if (singular(ii) > 0.0_ap) then
                diagonal(ii, ii) = cmplx(log(singular(ii)), 0.0_dp, dp)
            else 
                diagonal(ii, ii) = (0.0_dp, 0.0_dp)
            end if
            if (singular(ii) > 1.0_ap) then
                diagonal(ii,ii) = (0.0_dp, 0.0_dp)
            end if
        end do

        !! Calculates the matrix multiplication inside the trace in the von Neumann entropy
        diagonal = matmul(cmplx(density_matrix, kind=ap), matmul(left, matmul(diagonal, right)))
        
        !! Calculates the negative trace
        vonNeumann_entropy = 0
        do ii = 1, nn  
                vonNeumann_entropy = vonNeumann_entropy - real(diagonal(ii, ii), dp)
        end do

    end function vonNeumann_entropy
    
    function hamilton_exp(var, matrix) result(dev)
        !! Declarations
        complex(dp), intent(in) :: matrix(:,:), var
        integer :: nn, ii, jj, kk, steps, ii_err
        complex(dp), allocatable ::  dev(:,:), basis(:,:,:), potenz(:,:)
        real(dp), allocatable :: coeff(:,:)
        logical :: overflow
        real(dp) :: over
        overflow = .false.
        ii_err = 0
        steps = 500
        
        nn = size(matrix, dim=1)
        allocate(dev(nn,nn))
        allocate(basis(nn,nn,nn))
        allocate(coeff(nn,2))
        allocate(potenz(nn,nn))
        
        if (nn < 2) then
            write(*,*) "Error: Wrong size of input matrix in hamilton_exp. Stop"
            stop
        end if
        
        !! Create basis matrices
        basis = (0.0_dp, 0.0_dp)  
        do ii = 1, nn
            if (modulo(ii,2) == 0) then !! Even Numbers
                do jj = (ii/2), nn-(ii/2)
                    do kk = 0, (ii/2)-1
                        basis(jj+1+kk, jj-kk, ii) = (1.0_dp, 0.0_dp)
                        basis(jj-kk, jj+1+kk, ii) = (1.0_dp, 0.0_dp)
                    end do
                end do
            else                        !! Odd Numbers
                do jj = (ii/2 + 1), nn-(ii/2)
                    do kk = 0, (ii/2 + 1)-1
                        basis(jj+kk, jj-kk, ii) = (1.0_dp, 0.0_dp)
                        basis(jj-kk, jj+kk, ii) = (1.0_dp, 0.0_dp)
                    end do
                end do
            end if
        end do
        
        call write_array_real("basis.dat", real(basis,dp))
        coeff = 0.0_dp
        coeff(1,1) = 1.0_dp
        coeff(2,1) = 1.0_dp
        if (nn > 2) then 
            coeff(3,1) = 1.0_dp
        end if
        
        dev = (0.0_dp, 0.0_dp)
        do ii = 1, nn !! 0. potenz 
            dev(ii,ii) = (1.0_dp, 0.0_dp)
        end do
        
        call write_vec_real("coeff.dat", coeff(:,1), horizontal=.true.)
        over = 10000_dp

        do ii = 0, steps
             !! calculate odd potenz of matrix 1, 3, 5, ..
            potenz = (0.0_dp, 0.0_dp)
            do jj = 1, nn/2
                potenz = potenz + coeff(2*jj,1) *  basis(:,:,2*jj)
            end do
            if ((abs(var**(2*ii+1)) == 0.0_dp) ) then
                exit
            end if
            dev = dev + ((var)**(2*ii+1))/factorial(2*ii+1) * potenz

            !! calculate even potenz of matrix 2, 4, 6, ..
            potenz = (0.0_dp, 0.0_dp)
            do jj = 0, nn/2-modulo(nn+1,2)
                potenz = potenz + coeff(2*jj+1,1) *  basis(:,:,2*jj+1)
            end do
            dev = dev + ((var)**(2*(ii+1)))/factorial(2*(ii+1)) * potenz

            !! calculate next coefficient
            if (nn > 2) then
                coeff(1,2) = coeff(1,1) + coeff(3,1)
                if (nn > 3) then
                    coeff(2,2) = 2*coeff(2,1) + coeff(4,1)
                else
                    coeff(2,2) = 2*coeff(2,1)
                end if
                do jj = 3, nn-2
                    coeff(jj,2) = coeff(jj-2,1) + 2*coeff(jj,1) + coeff(jj+2,1)
                end do
                if (nn > 3) then
                    coeff(nn-1,2) = coeff(nn-3,1) + 2*coeff(nn-1,1)
                else
                    coeff(nn-1,2) = 2*coeff(nn-1,1)
                end if
                coeff(nn,2) = coeff(nn-2,1) + coeff(nn,1)
            end if
            coeff(:,1) = coeff(:,2)
            if (any(coeff(:,1) < 0)) then
                overflow = .true.
                coeff = 0
                if (ii_err == 0) then
                    ii_err = ii
                end if
            end if
            call write_vec_real("coeff.dat", coeff(:,1), horizontal=.true., state_in="old", pos_in="append")
        end do
        
        if (overflow) then
            write(*,*) "Overflow in hamilton_exp in ii =",ii_err
        end if
        
        call write_matrix_complex("dev.dat", dev)
        call write_matrix_complex("udev.dat", matmul(dev,transpose(conjg(dev))))
    end function hamilton_exp


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
        term = (0.0_dp, 0.0_dp)
        
        do ii = 1, n_env
            term = term + matmul(environment(:,:,ii), matmul(density_matrix, transpose(conjg(environment(:,:,ii)))))
        end do
        
     end function environment_term  
     
     logical function is_hermitian(matrix)
        complex(dp), intent(in) :: matrix(:,:)
        
        if ( any( abs((matrix - transpose(conjg(matrix)))) >= err) ) then
            is_hermitian = .false.
        else
            is_hermitian = .true.
        end if
        
    end function is_hermitian
   
    real(dp) function factorial(n)
        integer :: ii
        integer, intent(in) :: n
        
        factorial = 1.0_dp
        if ( n /= 0) then
            do ii = 1, n
                factorial = factorial*ii
            end do
        end if
        
        
    end function factorial
    
    integer(dp) function factorial_int(n)
        integer :: ii
        integer, intent(in) :: n
        
        factorial_int = 1_dp
        if ( n /= 0) then
            do ii = 1, n
                factorial_int = factorial_int*ii
            end do
        end if
        
        
    end function factorial_int
    
    function five(timestep) result(ex)
        real(dp), intent(in) :: timestep
        integer :: ii
        complex(dp) :: left(5,5), singular(5,5), right(5,5), ex(5,5)
        
        singular = (0.0_dp, 0.0_dp)
        left = (0.0_dp, 0.0_dp)
        left(1,1) = (-1.0_dp, 0.0_dp)
        left(2,1) = (1.0_dp, 0.0_dp)
        left(4,1) = (-1.0_dp, 0.0_dp)
        left(5,1) = (1.0_dp, 0.0_dp)
        left(1,2) = (1.0_dp, 0.0_dp)
        left(3,2) = (-1.0_dp, 0.0_dp)
        left(5,2) = (1.0_dp, 0.0_dp)
        left(1,3) = (-1.0_dp, 0.0_dp)
        left(2,3) = (-1.0_dp, 0.0_dp)
        left(4,3) = (1.0_dp, 0.0_dp)
        left(5,3) = (1.0_dp, 0.0_dp)
        left(1,4) = (1.0_dp, 0.0_dp)
        left(2,4) = cmplx(-sqrt(3.0_dp), 0.0_dp, dp)
        left(3,4) = (2.0_dp, 0.0_dp)
        left(4,4) = cmplx(-sqrt(3.0_dp), 0.0_dp, dp)
        left(5,4) = (1.0_dp, 0.0_dp)
        left(1,5) = (1.0_dp, 0.0_dp)
        left(2,5) = cmplx(sqrt(3.0_dp), 0.0_dp, dp)
        left(3,5) = (2.0_dp, 0.0_dp)
        left(4,5) = cmplx(sqrt(3.0_dp), 0.0_dp, dp)
        left(5,5) = (1.0_dp, 0.0_dp)
        
        right = (0.0_dp, 0.0_dp)
        right(1,1) = (-0.25_dp, 0.0_dp)
        right(2,1) = cmplx(1.0_dp/3, 0.0_dp,dp)
        right(3,1) = (-0.25_dp, 0.0_dp)
        right(4,1) = cmplx(1.0_dp/12, 0.0_dp, dp)
        right(5,1) = cmplx(1.0_dp/12, 0.0_dp, dp)
        right(1,2) = (0.25_dp, 0.0_dp)
        right(3,2) = (-0.25_dp, 0.0_dp)
        right(4,2) = cmplx(-1.0_dp/(4.0_dp*sqrt(3.0_dp)), 0.0_dp, dp)
        right(5,2) = cmplx(1.0_dp/(4.0_dp*sqrt(3.0_dp)), 0.0_dp, dp)
        right(1,3) = (0.0_dp, 0.0_dp)
        right(2,3) = cmplx(-1.0_dp/3, 0.0_dp,dp)
        right(4,3) = cmplx(1.0_dp/6, 0.0_dp,dp)
        right(5,3) = cmplx(1.0_dp/6, 0.0_dp,dp)
        right(1,4) = (-0.25_dp, 0.0_dp)
        right(2,4) = (0.0_dp, 0.0_dp)
        right(3,4) = (0.25_dp, 0.0_dp)
        right(4,4) = cmplx(-1.0_dp/(4.0_dp*sqrt(3.0_dp)), 0.0_dp, dp)
        right(5,4) = cmplx(1.0_dp/(4.0_dp*sqrt(3.0_dp)), 0.0_dp, dp)
        right(1,5) = (0.25_dp, 0.0_dp)
        right(2,5) = cmplx(1.0_dp/3, 0.0_dp,dp)
        right(3,5) = (0.25_dp, 0.0_dp)
        right(4,5) = cmplx(1.0_dp/12, 0.0_dp,dp)
        right(5,5) = cmplx(1.0_dp/12, 0.0_dp,dp)
        
        singular(1,1) = (-1.0_dp, 0.0_dp)
        singular(2,2) = (0.0_dp, 0.0_dp)
        singular(3,3) = (1.0_dp, 0.0_dp)
        singular(4,4) = cmplx(-sqrt(3.0_dp), 0.0_dp, dp)
        singular(5,5) = cmplx(sqrt(3.0_dp), 0.0_dp, dp)
        call write_matrix_complex("lef.dat", matmul(left, right) )
        do ii = 1, 5
            singular(ii,ii) = exp(cmplx(0.0_dp, -timestep, dp)*singular(ii,ii))
        end do
        ex = matmul(left, matmul(singular, right))
        
        call write_matrix_complex("dev.dat", ex)
        call write_matrix_complex("udev.dat", matmul(ex, transpose(conjg(ex))) )
        
    end function five
    
    !! Calculates the development matrix of an 2 parallel 4 nodes line network
    function special_2(var) result(ex)
        complex(dp), intent(in) :: var
        complex(dp) :: ex(6,6)
        complex(dp), parameter :: i = (0.0_dp, 1.0_dp)

        ex = 1.0_dp/3 * transpose(reshape( (/  2.0_dp*cos(var)+cos((2.0_dp*var)) , i*(sin(var)+sin(((2.0_dp*var)))) , &
                    & i*(sin(var)+sin((2.0_dp*var))) , cos((2.0_dp*var))-cos(var) , &
                    & cos((2.0_dp*var))-cos(var) , -i*(2.0_dp*sin(var)-sin((2.0_dp*var))),&
                  &  i *(sin(var)+sin((2.0_dp*var))) , 2.0_dp *cos(var)+cos((2.0_dp*var))  , &
                    &  cos((2.0_dp*var))-cos(var) , i* (sin(var)+sin((2.0_dp*var))) , &
                    & -i*(2.0_dp* sin(var)-sin((2.0_dp*var))) , cos((2.0_dp*var))-cos(var), &
                  &  i* (sin(var)+sin((2.0_dp*var))) , cos((2.0_dp*var))-cos(var)  , &
                    &  2.0_dp*cos(var)+cos((2.0_dp*var)) , -i*(2.0_dp*sin(var)-sin((2.0_dp*var))) , &
                    & i*(sin(var)+sin((2.0_dp*var))) , cos((2.0_dp*var))-cos(var) , &
                  &   cos((2.0_dp*var))-cos(var) , i*(sin(var)+sin((2.0_dp*var))) , &
                    & -i*(2.0_dp*sin(var)-sin((2.0_dp*var))) , 2.0_dp*cos(var)+cos((2.0_dp*var)) , &
                    & cos((2.0_dp*var))-cos(var) , i*(sin(var)+sin((2.0_dp*var))) , &
                  &   cos((2.0_dp*var))-cos(var) , -i*(2.0_dp*sin(var)-sin((2.0_dp*var))) , &
                    & i*(sin(var)+sin((2.0_dp*var))) , cos((2.0_dp*var))-cos(var) , &
                    & 2.0_dp*cos(var)+cos((2.0_dp*var)) , i*(sin(var)+sin((2.0_dp*var))) , &
                  &   -i*(2.0_dp*sin(var)-sin((2.0_dp*var))) , cos((2.0_dp*var))-cos(var) ,&
                    & cos((2.0_dp*var))-cos(var) , i*(sin(var)+sin((2.0_dp*var))) , &
                    & i*(sin(var)+sin((2.0_dp*var))) , 2.0_dp*cos(var)+cos((2.0_dp*var)) /), [6,6] ))
        call write_matrix_complex("ex.dat", ex)            
        call write_matrix_complex("uex.dat", matmul(ex, transpose(conjg(ex))))
     end function special_2
     
     function special_3(var) result(ex)
        complex(dp), intent(in) :: var
        complex(dp) :: ex(8,8)
        real(dp) :: coeff(3), ws(8,8)
        integer :: nn, ii
        
        nn = 8
        
        coeff = 1.0_dp
        ws = 0.0_dp
        do ii = 1, nn
            ws(ii, ii) = 1.0_dp
        end do
        
        do ii = 0, 100
            ex = ex + (var**ii)/factorial(ii) *ws
            
            if (modulo(ii+1,2) == 0) then
                
            end if
        end do
    end function special_3
    
    function man_exp(matrix, var) result(ex)
        complex(dp), intent(in) :: var
        complex(dp), intent(in) :: matrix(:,:)
        complex(dp), allocatable :: ex(:,:), ws(:,:)
        integer :: nn, ii, jj
        
        nn = size(matrix, dim=1)
        allocate(ex(nn,nn))
        allocate(ws(nn,nn))
        
        ws = matrix
        ex = (0.0_dp, 0.0_dp)

        do ii = 0, 100
            ws = (0.0_dp, 0.0_dp)
            do jj = 1, nn
                ws(jj,jj) = (1.0_dp, 0.0_dp)
            end do
            do jj = 1, ii
                ws = matmul(ws, matrix)
            end do
            ex = ex + (var**ii)/factorial(ii) * ws
        end do
        
        call write_matrix_complex("man.dat", ex)
        call write_matrix_complex("uman.dat", matmul(ex, transpose(conjg(ex))))
    end function man_exp
end module quantumWalk
