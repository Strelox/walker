!> Sets accuracy of calculation
!!
!! \param dp  sets double precision
!! \param err  sets lower accuracy bound
!!
module accuracy
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 300)
  real, parameter :: err = 1e-12

end module accuracy
