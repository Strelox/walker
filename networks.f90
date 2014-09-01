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

!> Provides network related routines
!!
module networks
    use accuracy
    use io
    implicit none 

contains
	
    !> Normalize the entries in a network array.
    !!
    !! \param network Array of the network that will be normalized.
    !!
    subroutine normal(network)
        real(dp), intent(inout) :: network(:,:,:)
        integer :: ii, jj, kk, n1, n2, n3
        n1 = size(network, dim=1)
        n2 = size(network, dim=2)
        n3 = size(network, dim=3)

        do ii=1, n1
            do jj=1, n2
                    do kk=1, n3
                         if(network(ii,jj,kk) /= 0) then
                                network(ii,jj,kk) = 1/network(ii,jj,kk)
                         end if
                    end do
             end do
        end do

    end subroutine normal
    
    
    !> Calculating distances of vetices to start
    !!
    !! \param laplacian Laplacian matrix containing the information about the network.
    !! \param distance    Vector containing the distance to the start point of each vortex.
    !!
    subroutine calc_distance(laplacian, distance)
        real(dp), intent(in) :: laplacian(:,:)
        integer, allocatable, intent(out) :: distance(:)
        integer :: ii, jj, kk, n1
        
        n1 = size(laplacian, dim=1)
        allocate(distance(n1))
        distance = 0
  
        do kk = 2, n1   ! Calculates connections to start node
            if (laplacian(1,kk) < 0.0) then
                distance(kk) = 1
            end if
        end do
        
        do ii = 1, n1   ! Distance to start
            do jj = 2, n1   ! First one ist start node
                do kk = 2, n1
                    if ( (distance(jj) == ii) .and. (laplacian(jj,kk) < 0.0) .and. (distance(kk) == 0) ) then
                            distance(kk) = distance(jj) + 1
                    end if
                end do
            end do
        end do
    end subroutine calc_distance
    
end module networks