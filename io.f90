!> Provides input/output related routines
!!
module io
  use accuracy
  implicit none

contains
  !> Reads an 3 dim array
  !!
  subroutine read_array(data, n1, n2, n3, array)

    !! Declarations
    integer :: ii, jj, kk, status
    character(*), intent(in) :: data
    integer, intent(out) :: n1, n2, n3
    real(dp), allocatable, intent(out) :: array(:,:,:)

    !! Reads size of array
    open(11, file=data, status="old", form="formatted", action="read")
    read(11,*, iostat=status) n1, n2, n3
    if (status /= 0) then
       write(*,"(A)") "incorrect input at size"
       stop
    end if

    !! Reads array
    allocate(array(n1,n2,n3))
    do kk = 1, n3
       do ii = 1, n1
          read(11,*, iostat=status) (array(ii,jj,kk), jj=1,n2)
          if (status /= 0) then
             write(*,"(A)") "incorrect input in array"
             stop
          end if
       end do
       read(11,*)
    end do

    close(11)
  end subroutine read_array

  !> Reads a square matrix from a file
  !!
  subroutine read_matrix(data, array)

    !! Declarations
    integer :: ii, jj, kk, status, n1, n2
    character(*), intent(in) :: data
    real(dp), allocatable, intent(out) :: array(:,:)

    !! Reads size of matrix
    open(11, file=data, status="old", form="formatted", action="read")
    read(11,*, iostat=status) n1, n2
    if (status /= 0) then
       write(*,"(A)") "incorrect input at size"
       stop
    end if

    !! Reads matrix
    allocate(array(n1,n2))
    do ii = 1, n1
       read(11,*, iostat=status) (array(ii,jj), jj=1,n2)
       if (status /= 0) then
          write(*,"(A)") "incorrect input in array"
          stop
       end if
    end do

    close(11)
  end subroutine read_matrix

  !> Writes a matrix into a file
  !!
  subroutine write_matrix(data, array)

    integer :: ii, jj, nn, mm
    character(*), intent(in) :: data
    integer, intent(in) :: array(:,:)

    nn = size(array(:,1))
    mm = size(array(1,:))

    open(12, file=data, status="replace", form="formatted", action="write")
    do ii = 1, nn
       write(12, *) ii, (array(ii,jj), jj=1, mm)
    end do

    close(12)
  end subroutine write_matrix


  !> Writes an array into a file
  !!
  subroutine write_array(data, array)

    integer :: ii, jj, kk, n1, n2, n3
    character(*), intent(in) :: data
    real(dp), intent(in) :: array(:,:,:)

    n1 = size(array, DIM = 1)
    n2 = size(array, DIM = 2)
    n3 = size(array, DIM = 3)

    open(13, file=data, status="replace", form="formatted", action="write")
    do kk = 1, n3
       do ii = 1, n1
          write(13, *) ii, (array(ii,jj,kk), jj=1, n2)
       end do
       write(13, *)
    end do

    close(13)
  end subroutine write_array
end module io



