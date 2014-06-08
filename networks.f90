!> Provides network related routines
!!
module networks
    use accuracy
    use io
    use random
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
        integer :: ii, jj, n1
        
        n1 = size(laplacian, dim=1)
        allocate(distance(n1))
        distance = 0
         
        do ii = 1, n1
             do jj = 1, n1
                    if ((laplacian(ii,jj) > 0.0) .and. (jj /= 1)) then
                         if ((distance(jj) == 0).or.(distance(jj) > distance(ii) + 1)) then
                                distance(jj) = distance(ii) + 1
                         end if
                    end if
             end do
        end do
        
    end subroutine calc_distance
    
end module networks