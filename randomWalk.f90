module randomWalk
  use accuracy
  use io
  use random
  implicit none
  
contains
  
  !> General Random Walk through Network
  !!
  subroutine gen_rWalk(network,pos,Time,nwalks,counter)

    real(dp), intent(in) :: network(:,:,:)
    integer, intent(in) :: Time, nwalks
    integer, intent(inout) :: pos(:,:), counter(:)
    integer :: ii, jj, kk, ll, n1, n2, n3
    real(dp) :: rand, temp, total
    n1 = size(network, dim=1)
    n2 = size(network, dim=2)
    n3 = size(network, dim=3)
    
    !! Random Walk
    do ii = 2, Time
       do jj = 1, nwalks
          !! Is the walker already out of the wire?
          if (pos(ii-1, jj) == n1) then
             pos(ii, jj) = n1
             cycle
          end if

          !! Call a random number 0<random<=1
          do
             call random_number(rand)
             if (rand /= 0.0) then
                exit
             end if
          end do

          !! Calculation of next node
          total = sum(network(pos(ii-1,jj),:,:))
          temp = 0
          loop: do kk = 1, n3
             do ll = 1, n2
                temp = temp + network(pos(ii-1, jj), ll, kk)/total
                if (temp >= rand) then
                   pos(ii, jj) = ll
                   counter(ll) = counter(ll) + 1
                   exit loop
                end if
             end do
          end do loop

       end do
    end do

  end subroutine gen_rWalk

  
end module randomWalk

