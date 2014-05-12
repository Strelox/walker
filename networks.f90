!> Provides network related routines
!!
module networks
  use accuracy
  use io
  use random
  implicit none 

contains

  !! Normalization (???)
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

end module networks
