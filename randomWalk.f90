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

!> Provides the routines for the classical random walk.
!!
module randomWalk
    use accuracy
    use io
    implicit none

contains

    !> Classical Random Walk through network. 
    !!
    !! \details The routine performs a classical random walk through a given
    !! network.
    !!
    !! \param laplacian     Laplacian matrix containing information about the network.
    !! \param vertex        Vector to contain the probability to be in each vertex. 
    !! \param io_rates      Input/Output rates
    !! \param timestep      The size of a timestep.
    !! \param maxTimestep   The maximum amount of timesteps that will be simulated.
    !!
    subroutine CRWalk(laplacian, vertex, io_rates, timestep)

        !! Declarations
        real(dp), intent(in) :: timestep, io_rates(:)
        real(dp), intent(inout) :: laplacian(:,:)
        real(dp), intent(inout) :: vertex(:)
        integer :: n1
        real(dp), allocatable :: rate(:)
        n1 = size(laplacian, dim=1)
        

        !! Calculate rate
        allocate(rate(n1))
        rate = matmul(laplacian,vertex)
        
        !! Lets walk
        vertex = vertex - rate*timestep
        
        !! Add input/output rates
        vertex(1) = vertex(1) + io_rates(1)*timestep
        vertex(size(vertex)) = vertex(size(vertex)) + io_rates(2)*timestep
           
        if (minval(vertex) < 0) then
            where(vertex < 0)
                vertex = 0
            end where
        end if

    end subroutine CRWalk
    
    !> Calculates entropy of a probability vector
    !!
    !! \param vec Probability vector
    !!
    real(dp) function calc_entropy(vec)
        real(dp), intent(in) :: vec(:)
        integer :: ii
        
        calc_entropy = 0.0_dp
        do ii = 1, size(vec)
            if (vec(ii) /= 0) then
                calc_entropy = calc_entropy - vec(ii)/sum(vec) * log(vec(ii)/sum(vec))
            end if
        end do
        
    end function calc_entropy
    
    !> Calculates current
    subroutine calc_current(current, laplacian, vertex, io_rates, distance)
        real(dp), intent(in) :: laplacian(:,:), vertex(:), io_rates(2) 
        integer, intent(in) :: distance(:)
        integer :: ii, jj
        real(dp), allocatable, intent(out) :: current(:)

        !! Calculates current from start to end
        allocate(current(size(vertex)+1))
        current = 0.0_dp
        current(1) = io_rates(1)
        do ii = 1, size(vertex)-1
            do jj = 1, size(vertex)
                if (laplacian(ii,jj) /= 0.0_dp) then
                    if (distance(jj) == distance(ii) + 1) then
                        current(ii+1) = current(ii+1) + laplacian(jj,ii)*vertex(ii) - laplacian(ii,jj) * vertex(jj) 
                    end if
                end if
            end do
        end do
        current(ubound(current)) = -io_rates(2)
        
    end subroutine calc_current
    
    !> Calculates current from end to start (recurrent)
    subroutine calc_recurrent(recurrent, laplacian, vertex, io_rates, distance)
        real(dp), intent(in) :: laplacian(:,:), vertex(:), io_rates(2) 
        integer, intent(in) :: distance(:)
        integer :: ii, jj
        real(dp), allocatable, intent(out) :: recurrent(:)
        
        !! Calculates recurrent from end to start
        allocate(recurrent(size(vertex)+1))
        recurrent = 0
        recurrent(1) = -io_rates(1)
        do ii = 2, size(vertex)
            do jj = 1, size(vertex)
                if (laplacian(ii,jj) /= 0) then
                    if (distance(jj) == distance(ii) - 1) then
                        recurrent(ii) = recurrent(ii) + laplacian(jj,ii)*vertex(ii) - laplacian(ii,jj) * vertex(jj) 
                    end if
                end if
            end do
        end do
        recurrent(ubound(recurrent)) = io_rates(2)
      
    end subroutine calc_recurrent
                        
        
end module randomWalk