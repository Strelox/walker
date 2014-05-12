!> Provides input/output related routines
!!
module io
  use accuracy
  implicit none

contains

  !> Reads Configuration of the random Walk
  subroutine read_config(data, nWalker, maxTimestep, mode, start, timestep)
    character(*), intent(in) :: data
    integer, intent(out) :: nWalker, maxTimestep, start
    real(dp), intent(out) :: timestep
    character(*), intent(out) :: mode
    integer :: status

    open(14, file=data, status="old", form="formatted", action="read")
    read(14,*, iostat=status) nWalker
    if (status /= 0) then
       write(*,*) "Error: Wrong number of walker format in config file! End program."
       stop
    end if

    read(14,*, iostat=status) maxTimestep
    if (status /= 0) then
       write(*,*) "Error: Wrong maxTimestep steps format in config file! End program."
       stop
    end if

    read(14,*, iostat=status) mode
    if (status /= 0) then
       write(*,*) "Error: Wrong mode format in config file! End program."
       stop
    end if

    read(14,*, iostat=status) start
    if (status /= 0) then
       write(*,*) "Error: Wrong start format in config file! End program."
       stop
    end if
    
    read(14,*, iostat=status) timestep
    if (status /= 0) then
       write(*,*) "Error: Wrong timestep format in config file! End program."
       stop
    end if
    
    close(14)
  end subroutine read_config

  !> Reads a (integer) matrix from a file
  !!
  subroutine read_matrix_int(data, array)

    !! Declarations
    integer :: ii, jj, status, n1, n2
    character(*), intent(in) :: data
    integer, allocatable, intent(out) :: array(:,:)

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
  end subroutine read_matrix_int
  
  !> Reads a (real) matrix from a data
  !!
  subroutine read_matrix_real(data, array)

    !! Declarations
    integer :: ii, jj, status, n1, n2
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
  end subroutine read_matrix_real
  
  !> Reads a (complex) matrix from a file
  !!
  subroutine read_matrix_complex(data, array)

    !! Declarations
    integer :: ii, jj, status, n1, n2
    character(*), intent(in) :: data
    complex(dp), allocatable, intent(out) :: array(:,:)

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
  end subroutine read_matrix_complex
  
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
  
  !> Writes a (real) matrix into a file
  !!
  subroutine write_matrix_int(data, array)

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
  end subroutine write_matrix_int

  !> Writes a (real) matrix into a file
  !!
  subroutine write_matrix_real(data, array)

    integer :: ii, jj, nn, mm
    character(*), intent(in) :: data
    real(dp), intent(in) :: array(:,:)

    nn = size(array(:,1))
    mm = size(array(1,:))

    open(12, file=data, status="replace", form="formatted", action="write")
    do ii = 1, nn
       write(12, *) ii, (array(ii,jj), jj=1, mm)
    end do

    close(12)
  end subroutine write_matrix_real


  !> Writes a (real) matrix into a file
  !!
  subroutine write_matrix_complex(data, array)

    integer :: ii, jj, nn, mm
    character(*), intent(in) :: data
    complex(dp), intent(in) :: array(:,:)

    nn = size(array(:,1))
    mm = size(array(1,:))

    open(12, file=data, status="replace", form="formatted", action="write")
    do ii = 1, nn
       write(12, *) ii, (array(ii,jj), jj=1, mm)
    end do

    close(12)
  end subroutine write_matrix_complex
  
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



  !> Writing data by the line
  !!
  !! \param vec  Containing the data to write out
  !! \param data  Name of the outputfile
  !!
  subroutine write_x(data, vec)
    real(dp), intent(in) :: vec(:)
    character(*), intent(in) :: data
    integer :: ii

    open(12, file=data, status="replace", form="formatted", action="write")
    do ii = 1, size(vec)
       write(12, "(ES20.12)") vec(ii)
    end do
    close(12)
  end subroutine write_x

  subroutine write_x_int(data, vec)
    integer, intent(in) :: vec(:)
    character(*), intent(in) :: data
    integer :: ii

    open(12, file=data, status="replace", form="formatted", action="write")
    do ii = 1, size(vec)
       write(12, "(I6)") vec(ii)
    end do
    close(12)
  end subroutine write_x_int

  !> Writes data in XY-format
  !!
  !! \param xx Containing x-data
  !! \param yy Containing y-data
  !! \param data Name of the outputfile
  !!
  subroutine write_xy(data, xx, yy)
    real(dp), intent(in) :: xx(:), yy(:)
    character(*), intent(in) :: data
    integer :: ii

    open(13, file=data, status="replace", form="formatted", action="write")
    do ii = 1, size(xx)
       write(13, "(2(ES20.12))") xx(ii), yy(ii)
    end do
    close(13)

  end subroutine write_xy


  !> Writing data in NXY-format
  !!
  !! \param xx  Containing x-data
  !! \param yy  Containing y-data column-wise
  !! \param data  Name of the outputfile
  !!
  subroutine write_nxy(data, xx, yy)
    real(dp), intent(in) :: xx(:), yy(:,:)
    character(*), intent(in) :: data
    integer :: ii, nx, ny
    character(25) :: form

    nx = size(xx)
    ny = size(yy(1,:))
    write(form, "(A1, I0, A)") "(", (ny+1), "(ES20.12, 2X))"

    open(14, file=data, status="replace", form="formatted", action="write")
    do ii = 1, nx
       write(14, form) xx(ii), yy(ii, 1:ny)
    end do
    close(14)

  end subroutine write_nxy

end module io



