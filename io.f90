!> Provides input/output related routines
!!
module io
  use accuracy
  implicit none

contains

!-------------------Input-------------------------------------------------

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

    !! Reads matrix50
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
  
!---------------------------Output------------------------------

  !> Writes a (integer) matrix into a file
  !!
  subroutine write_matrix_int(data, array, lineNr, size_out, state_in, pos_in)

    integer :: ii, jj, nn, mm
    character(*), intent(in) :: data
    character(30) :: form, state, pos
    real(dp), intent(in) :: array(:,:) 
    logical, optional, intent(in) :: lineNr, size_out
    character(*), optional, intent(in) :: state_in, pos_in
    
    nn = size(array(:,1))
    mm = size(array(1,:))
    
    if (present(state_in) .eqv. .true.) then
      if ((state_in /= "new") .or. (state_in /= "old") .or. (state_in /= "replace") .or. (state_in /= "scratch")) then
        write(*,*) "Wrong file status in write_matrix_int. End program."
        stop
      end if
      state = state_in
    else
      state = "replace"
    end if
    
    if (present(pos_in) .eqv. .true.) then
      if (( pos_in == "append") .or. (pos_in == "rewind") .or. (pos_in == "asis")) then
        pos = pos_in
      else
        write(*,*) "Error: Wrong position input in write_matrix_int. Stop program."
      end if
    else
      pos = "rewind"
    end if
    
    open(12, file=data, status=state, form="formatted", action="write", position=pos)
    
    if ((present(size_out) .eqv. .true.) .and. (size_out .eqv. .true. )) then
      write(12, "(2(I0))") nn, mm
    end if
    
    if ((present(lineNr) .eqv. .true. ) .and. ( lineNr .eqv. .true.)) then
      write(form, "(A3, I0, A)") "(I,", (nn), "(I15, 2X))"
      do ii = 1, nn
      write(12, form) ii, (array(ii,jj), jj=1, mm)
      end do
    else
      write(form, "(A1, I0, A)") "(", (nn), "(I15, 2X))"
      do ii = 1, nn
      write(12, *) (array(ii,jj), jj=1, mm)
      end do
    end if
    
    close(12)
  end subroutine write_matrix_int

  !> Writes a (real) matrix into a file
  !!
  subroutine write_matrix_real(data, array, lineNr, size_out, state_in, pos_in)
    
    integer :: ii, jj, nn, mm
    character(30) :: form, state, pos
    character(*), intent(in) :: data
    real(dp), intent(in) :: array(:,:) 
    logical, optional, intent(in) :: lineNr, size_out
    character(*), optional, intent(in) :: state_in, pos_in
    
    nn = size(array, dim=1)
    mm = size(array, dim=2)
    
    if (present(state_in) .eqv. .true.) then
      if ((state_in /= "new") .and. (state_in /= "old") .and. (state_in /= "replace") .and. (state_in /= "scratch")) then
        write(*,*) "Wrong file status in write_matrix_real. End program."
        stop
      end if
      state = state_in
    else
      state = "replace"
    end if
    
    if (present(pos_in) .eqv. .true.) then
      if (( pos_in == "append") .or. (pos_in == "rewind") .or. (pos_in == "asis")) then
        pos = pos_in
      else
        write(*,*) "Error: Wrong position input in write_matrix_real. Stop program."
      end if
    else
      pos = "rewind"
    end if
    
    open(12, file=data, status=state, form="formatted", action="write", position=pos)
    
    if ((present(size_out) .eqv. .true.) .and. (size_out .eqv. .true. ))then
      write(12, *) nn, mm
    end if

    if ((present(lineNr) .eqv. .true. ) .and. ( lineNr .eqv. .true.)) then
      write(form, "(A3, I0, A)") "(I,", (nn), "(ES15.6, 2X))"
      do ii = 1, nn
      write(12, form) ii, (array(ii,jj), jj=1, mm)
      end do
    else 
      write(form, "(A1, I0, A)") "(", (mm), "(ES15.6, 2X))"
      do ii = 1, nn
      write(12, form) (array(ii,jj), jj=1, mm)
      end do
    end if
    
    close(12)
  end subroutine write_matrix_real

  !> Writes an (real) array into a file
  !!
  subroutine write_array_real(data, array)

    integer :: ii, jj, kk, n1, n2, n3
    character(*), intent(in) :: data
    real(dp), intent(in) :: array(:,:,:)

    n1 = size(array, DIM = 1)
    n2 = size(array, DIM = 2)
    n3 = size(array, DIM = 3)

    open(13, file=data, status="replace", form="formatted", action="write")
    do kk = 1, n3
       do ii = 1, n1
          write(13, *) (array(ii,jj,kk), jj=1, n2)
       end do
       write(13, *)
    end do

    close(13)
  end subroutine write_array_real

  !> Writing data by the line
  !!
  !! \param vec  Containing the data to write out
  !! \param data  Name of the outputfile
  !!
  subroutine write_x(data, vec, size_out, state_in, pos_in)
    real(dp), intent(in) :: vec(:)
    character(*), intent(in) :: data
    character(*), optional, intent(in) :: state_in, pos_in
    logical, optional :: size_out
    integer :: ii
    character(10) :: state, pos

    if (present(state_in) .eqv. .true.) then
      if ((state_in /= "new") .or. (state_in /= "old") .or. (state_in /= "replace") .or. (state_in /= "scratch")) then
        write(*,*) "Wrong file status in write_x. End program."
        stop
      end if
      state = state_in
    else
      state = "replace"
    end if
    
    if (present(pos_in) .eqv. .true.) then
      if (( pos_in == "append") .or. (pos_in == "rewind") .or. (pos_in == "asis")) then
        pos = pos_in
      else
        write(*,*) "Error: Wrong position input in write_x. Stop program."
      end if
    else
      pos = "rewind"
    end if
    
    open(12, file=data, status=state, form="formatted", action="write", position=pos)
    if ((present(size_out)) .and. (size_out .eqv. .true.)) then
      write(12, "(I0)") size(vec)
    end if
    
    do ii = 1, size(vec)
       write(12, "(ES30.8)") vec(ii)
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