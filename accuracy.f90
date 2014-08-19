!> Sets accuracy of calculation
!!
!! \param dp  sets double precision
!! \param err  sets lower accuracy bound
!!
module accuracy
  implicit none

!   integer, parameter :: dp = selected_real_kind(15, 300)
  integer, parameter :: dp = selected_real_kind(33,4931)
    real(dp), parameter :: err = 1e-20

end module accuracy
