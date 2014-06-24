!> Provides the routines for each type of random walk.
!!
module randomWalk
    use accuracy
    use io
    use random
    implicit none

contains

    !> Classical Random Walk through network. 
    !!
    !! \details The routine performs a classical random walk through a given
    !! network.
    !!
    !! \param laplacian     Laplacian matrix containing information about the network.
    !! \param vortex        Vector to contain the probability to be in each vortex. 
    !! \param io_rates      Input/Output rates
    !! \param timestep      The size of a timestep.
    !! \param maxTimestep   The maximum amount of timesteps that will be simulated.
    !! \param time_mode     The time mode of the simulation (either discrete or continous).
    !!
    subroutine CRWalk(laplacian, vortex, io_rates, timestep, time_mode)

        !! Declarations
        real(dp), intent(in) :: timestep, io_rates(:)
        character(*), intent(in) :: time_mode
        real(dp), intent(inout) :: laplacian(:,:)
        real(dp), intent(inout) :: vortex(:)
        integer :: n1
        real(dp), allocatable :: rate(:)
        n1 = size(laplacian, dim=1)
        
        select case (time_mode)
            case default !! Unknown Mode. End program.
                write(*,*) "Error: Not known modus! Stop program."
                stop

            case ("discrete") !! discrete-time Random Walk mode
                write(*,*) "Sorry discrete-time mode is broken! Stop program!"
                stop
!                 !! Transform laplacian in transition matrix
!                 laplacian = laplacian * timestep
!                 do ii = 1, n1
!                     laplacian(ii,ii) = laplacian(ii,ii) + 1
!                 end do
!                 
!                 !! Add Input/Output rates
!                 do ii = 1, size(io_vertices)
!                     laplacian(io_vertices(ii), io_vertices(ii)) = laplacian(io_vertices(ii), io_vertices(ii)) + io_rates(ii)*timestep
!                 end do
!                 
!                 !! Lets walk        
!                 vortex = matmul(laplacian,vortex)                    

            case ("continous") !! continous-time Random Walk mode
         
            !! Calculate rate
            allocate(rate(n1))
            rate = matmul(laplacian,vortex)
            
            !! Lets walk
            vortex = vortex + rate*timestep
            
            !! Add input/output rates
            vortex(1) = vortex(1) + io_rates(1)*timestep
            vortex(size(vortex)) = vortex(size(vortex)) + io_rates(2)*timestep
            
            if (minval(vortex) < 0) then
!                 write(*,*) "Warning: Probability gets negative! I fix the issue but you may look for systematic errors."
                where(vortex < 0)
                    vortex = 0
                end where
            end if
        end select

    end subroutine CRWalk
    
    !> Calculates entropy of a probability vector
    !!
    !! \param vec Probability vector
    !!
    real(dp) function calc_entropy(vec)
        real(dp), intent(in) :: vec(:)
        integer :: ii
        
        calc_entropy = 0
        do ii = 1, size(vec)
            if (vec(ii) /= 0) then
                calc_entropy = calc_entropy - vec(ii)/sum(vec) * log(vec(ii)/sum(vec))
            end if
        end do
        
    end function calc_entropy
    
    !> Calculates current
    subroutine calc_current(current, laplacian, vortex, io_rates, distance)
        real(dp), intent(in) :: laplacian(:,:), vortex(:), io_rates(2) 
        integer, intent(in) :: distance(:)
        integer :: ii, jj
        real(dp), allocatable, intent(out) :: current(:)

        !! Calculates current from start to end
        allocate(current(size(vortex)+1))
        current = 0
        current(1) = io_rates(1)
        do ii = 1, size(vortex)-1
            do jj = 1, size(vortex)
                if (laplacian(ii,jj) /= 0) then
                    if (distance(jj) == distance(ii) + 1) then
                        current(ii+1) = current(ii+1) + laplacian(jj,ii)*vortex(ii) - laplacian(ii,jj) * vortex(jj) 
                    end if
                end if
            end do
        end do
        current(ubound(current)) = -io_rates(2)
        
    end subroutine calc_current
    
    subroutine calc_recurrent(recurrent, laplacian, vortex, io_rates, distance)
        real(dp), intent(in) :: laplacian(:,:), vortex(:), io_rates(2) 
        integer, intent(in) :: distance(:)
        integer :: ii, jj
        real(dp), allocatable, intent(out) :: recurrent(:)
        
        !! Calculates recurrent from end to start
        allocate(recurrent(size(vortex)+1))
        recurrent = 0
        recurrent(1) = -io_rates(1)
        do ii = 2, size(vortex)
            do jj = 1, size(vortex)
                if (laplacian(ii,jj) /= 0) then
                    if (distance(jj) == distance(ii) - 1) then
                        recurrent(ii) = recurrent(ii) + laplacian(jj,ii)*vortex(ii) - laplacian(ii,jj) * vortex(jj) 
                    end if
                end if
            end do
        end do
        recurrent(ubound(recurrent)) = io_rates(2)
      
    end subroutine calc_recurrent
                        
        
end module randomWalk