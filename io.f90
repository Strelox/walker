!> Provides input/output related routines
!!
module io
    use accuracy
    implicit none

contains

!-------------------Input-------------------------------------------------
    
    !> Reads a (integer) vector from a file
    !!
    !! \param data          file to read from
    !! \param vec           vector to write the data in
    !! \param vec_size      size of the vector (vec) if it is not written in the file
    !! \param horizontal    read data horizontal instead of vertical
    !!
    subroutine read_vec_int(data, vec, vec_size, horizontal)
        
        integer :: ii, status, n1
        character(*), intent(in) :: data
        integer, allocatable, intent(out) :: vec(:)
        integer, optional, intent(in) :: vec_size
        logical, optional, intent(in) :: horizontal
        
        open(11, file=data, status="old", form="formatted", action="read")
        if (present(vec_size) .eqv. .true.) then
            n1 = vec_size
        else
            read(11,*, iostat=status) n1
            if (status /= 0) then
                write(*,"(A)") "Incorrect input at size in read_vec_int"
                stop
            end if
        end if
        
        allocate(vec(n1))
        if ((present(horizontal) .eqv. .true.) .and. (horizontal .eqv. .true.)) then
            read(11,*, iostat=status) (vec(ii), ii=1, n1)
            if (status /= 0) then
                write(*,"(A)") "Incorrect input at file in read_vec_int"
                stop
            end if
        else
            do ii = 1, n1
                read(11,*, iostat=status) vec(ii)
                if (status /= 0) then
                    write(*,"(A)") "Incorrect input at file in read_vec_int"
                    stop
                end if
            end do
        end if
        
        close(11)
    end subroutine read_vec_int
    
    !> Reads a (real) vector from a file
    !!
    !! \param data          file to read from
    !! \param vec           vector to write the data in
    !! \param vec_size      size of the vector (vec) if it is not written in the file
    !! \param horizontal    read data horizontal instead of vertical
    !!
    subroutine read_vec_real(data, vec, vec_size, horizontal)
        
        integer :: ii, status, n1
        character(*), intent(in) :: data
        real(dp), allocatable, intent(out) :: vec(:)
        integer, optional, intent(in) :: vec_size
        logical, optional, intent(in) :: horizontal
        
        open(12, file=data, status="old", form="formatted", action="read")
        if (present(vec_size) .eqv. .true.) then
            n1 = vec_size
        else
            read(12,*, iostat=status) n1
            if (status /= 0) then
                write(*,"(A)") "Incorrect input at size in read_vec_real"
                stop
            end if
        end if
        
        allocate(vec(n1))
        if ((present(horizontal) .eqv. .true.) .and. (horizontal .eqv. .true.)) then
            read(12,*, iostat=status) (vec(ii), ii=1, n1)
            if (status /= 0) then
                write(*,"(A)") "Incorrect input at file in read_vec_real"
                stop
            end if
        else
            do ii = 1, n1
                read(12,*, iostat=status) vec(ii)
                if (status /= 0) then
                    write(*,"(A)") "Incorrect input at file in read_vec_real"
                    stop
                end if
            end do
        end if
    
        close(12)
    end subroutine read_vec_real
    
    !> Reads a (integer) matrix from a file
    !!
    !! \param data          file to read from
    !! \param matrix        matrix to write the data in
    !!
    subroutine read_matrix_int(data, matrix)

        !! Declarations
        integer :: ii, jj, status, n1, n2
        character(*), intent(in) :: data
        integer, allocatable, intent(out) :: matrix(:,:)

        !! Reads size of matrix
        open(13, file=data, status="old", form="formatted", action="read")
        read(13,*, iostat=status) n1, n2
        if (status /= 0) then
             write(*,"(A)") "incorrect input at size"
             stop
        end if

        !! Reads matrix
        allocate(matrix(n1,n2))
        do ii = 1, n1
             read(13,*, iostat=status) (matrix(ii,jj), jj=1,n2)
             if (status /= 0) then
                    write(*,"(A)") "incorrect input in matrix"
                    stop
             end if
        end do

        close(13)
    end subroutine read_matrix_int
    
    !> Reads a (real) matrix from a data
    !!
    !! \param data          file to read from
    !! \param matrix        matrix to write the data in
    !!
    subroutine read_matrix_real(data, matrix)

        !! Declarations
        integer :: ii, jj, status, n1, n2
        character(*), intent(in) :: data
        real(dp), allocatable, intent(out) :: matrix(:,:)

        !! Reads size of matrix
        open(14, file=data, status="old", form="formatted", action="read")
        read(14,*, iostat=status) n1, n2
        if (status /= 0) then
             write(*,"(A)") "incorrect input at size"
             stop
        end if

        !! Reads matrix
        allocate(matrix(n1,n2))
        do ii = 1, n1
             read(14,*, iostat=status) (matrix(ii,jj), jj=1,n2)
             if (status /= 0) then
                    write(*,"(A)") "incorrect input in matrix"
                    stop
             end if
        end do

        close(14)
    end subroutine read_matrix_real
    
    !> Reads a (complex) matrix from a file
    !!
    !! \param data          file to read from
    !! \param matrix        matrix to write the data in
    !!
    subroutine read_matrix_complex(data, matrix)

        !! Declarations
        integer :: ii, jj, status, n1, n2
        character(*), intent(in) :: data
        complex(dp), allocatable, intent(out) :: matrix(:,:)

        !! Reads size of matrix
        open(15, file=data, status="old", form="formatted", action="read")
        read(15,*, iostat=status) n1, n2
        if (status /= 0) then
             write(*,"(A)") "incorrect input at size"
             stop
        end if

        !! Reads matrix50
        allocate(matrix(n1,n2))
        do ii = 1, n1
             read(15,*, iostat=status) (matrix(ii,jj), jj=1,n2)
             if (status /= 0) then
                    write(*,"(A)") "incorrect input in matrix"
                    stop
             end if
        end do

        close(15)
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
        open(16, file=data, status="old", form="formatted", action="read")
        read(16,*, iostat=status) n1, n2, n3
        if (status /= 0) then
             write(*,"(A)") "incorrect input at size"
             stop
        end if

        !! Reads array
        allocate(array(n1,n2,n3))
        do kk = 1, n3
             do ii = 1, n1
                    read(16,*, iostat=status) (array(ii,jj,kk), jj=1,n2)
                    if (status /= 0) then
                         write(*,"(A)") "incorrect input in array"
                         stop
                    end if
             end do
             read(16,*)
        end do

        close(16)
    end subroutine read_array
    
!---------------------------Output------------------------------

    !> Writes a (integer) scalar into a file
    !!
    !! \param data          file to write in
    !! \param scalar        number which is to write
    !! \param state_in      write in new or old file or replace existing file
    !! \param pos_in        write from begin of the file or append to it
    !!
    subroutine write_int(data, scalar, state_in, pos_in)
        !! Declarations
        character(*), intent(in) :: data
        integer, intent(in)  :: scalar
        character(*), optional, intent(in) :: state_in, pos_in
        character(30) :: pos, state
        
        if (present(state_in) .eqv. .true.) then
            if ((state_in /= "new") .or. (state_in /= "old") .or. (state_in /= "replace") .or. (state_in /= "scratch")) then
                write(*,*) "Wrong file status in write_int. End program."
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
                write(*,*) "Error: Wrong position input in writ_int. Stop program."
            end if
        else
            pos = "rewind"
        end if
        
        open(17, file=data, status=state, form="formatted", action="write", position=pos)
        write(17, "(I0)") scalar
        close(17)     
        
    end subroutine write_int
    
    !> Writes a (real) scalar into a file
    !!
    !! \param data          file to write in
    !! \param scalar        number which is to write
    !! \param state_in      write in new or old file or replace existing file
    !! \param pos_in        write from begin of the file or append to it
    !!
    subroutine write_real(data, scalar, state_in, pos_in)
        !! Declarations
        character(*), intent(in) :: data
        real(dp), intent(in)  :: scalar
        character(*), optional, intent(in) :: state_in, pos_in
        character(30) :: pos, state
        
        if (present(state_in) .eqv. .true.) then
            if ((state_in /= "new") .and. (state_in /= "old") .and. (state_in /= "replace") .and. (state_in /= "scratch")) then
                write(*,*) "Wrong file status in write_real. End program."
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
                write(*,*) "Error: Wrong position input in writ_real. Stop program."
            end if
        else
            pos = "rewind"
        end if
        
        open(18, file=data, status=state, form="formatted", action="write", position=pos)
        write(18, "(ES15.6)") scalar
        close(18)     
        
    end subroutine write_real
    
    !> Writes a (integer) vector into a file
    !!
    !! \param data          file to write in
    !! \param vec           vector which is to write
    !! \param size_out      should the size of the vector be written in the first line of the file
    !! \param state_in      write in new or old file or replace existing file
    !! \param pos_in        write from begin of the file or append to it
    !! \param horizontal    should the vector be written horizontal instead of vertical?
    !!
    subroutine write_vec_int(data, vec, size_out, state_in, pos_in, horizontal)

        integer :: ii, nn
        character(*), intent(in) :: data
        character(30) :: form, state, pos
        integer, intent(in) :: vec(:) 
        logical, optional, intent(in) :: size_out, horizontal
        character(*), optional, intent(in) :: state_in, pos_in
        
        nn = size(vec)
        
        if (present(state_in) .eqv. .true.) then
            if ((state_in /= "new") .or. (state_in /= "old") .or. (state_in /= "replace") .or. (state_in /= "scratch")) then
                write(*,*) "Wrong file status in write_vec_int. End program."
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
                write(*,*) "Error: Wrong position input in write_vec_int. Stop program."
            end if
        else
            pos = "rewind"
        end if
        
        open(19, file=data, status=state, form="formatted", action="write", position=pos)
        
        if ((present(size_out) .eqv. .true.) .and. (size_out .eqv. .true. )) then
            write(19, "(I0)") nn
        end if
        
        if (horizontal .eqv. .true.) then
            write(form, "(A1, I0, A)") "(", (nn), "(I15, 2X))"
            write(19, form) (vec(ii), ii=1, nn)
        else
            do ii = 1, nn
            write(19, "(A15)") vec(ii)
            end do
        end if
        
        close(19)
    end subroutine write_vec_int
    
    !> Writes a (real) vector into a file
    !!
    !! \param data          file to write in
    !! \param vec           vector which is to write
    !! \param size_out      should the size of the vector be written in the first line of the file
    !! \param state_in      write in new or old file or replace existing file
    !! \param pos_in        write from begin of the file or append to it
    !! \param horizontal    should the vector be written horizontal instead of vertical?
    !!
    subroutine write_vec_real(data, vec, size_out, state_in, pos_in, horizontal)

        integer :: ii, nn
        character(*), intent(in) :: data
        character(30) :: form, state, pos
        real(dp), intent(in) :: vec(:) 
        logical, optional, intent(in) :: size_out, horizontal
        character(*), optional, intent(in) :: state_in, pos_in
        
        nn = size(vec)
        
        if (present(state_in) .eqv. .true.) then
            if ((state_in /= "new") .or. (state_in /= "old") .or. (state_in /= "replace") .or. (state_in /= "scratch")) then
                write(*,*) "Wrong file status in write_vec_int. End program."
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
                write(*,*) "Error: Wrong position input in write_vec_int. Stop program."
            end if
        else
            pos = "rewind"
        end if
        
        open(20, file=data, status=state, form="formatted", action="write", position=pos)
        
        if ((present(size_out) .eqv. .true.) .and. (size_out .eqv. .true. )) then
            write(20, "(I0)") nn
        end if
        
        if (horizontal .eqv. .true.) then
            write(form, "(A1, I0, A)") "(", (nn), "(ES15.6, 2X))"
            write(20, form) (vec(ii), ii=1, nn)
        else
            do ii = 1, nn
            write(20, "(ES15.6)") vec(ii)
            end do
        end if
        
        close(20)
    end subroutine write_vec_real
    
    !> Writes a (integer) matrix into a file
    !!
    !! \param data          file to write in
    !! \param matrix        matrix which is to write
    !! \param lineNr        should the line number be written in the first column?
    !! \param size_out      should the size of the matrix be written in the first line of the file
    !! \param state_in      write in new or old file or replace existing file
    !! \param pos_in        write from begin of the file or append to it
    !!
    subroutine write_matrix_int(data, array, lineNr, size_out, state_in, pos_in)

        integer :: ii, jj, nn, mm
        character(*), intent(in) :: data
        character(30) :: form, state, pos
        integer, intent(in) :: array(:,:) 
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
        
        open(21, file=data, status=state, form="formatted", action="write", position=pos)
        
        if ((present(size_out) .eqv. .true.) .and. (size_out .eqv. .true. )) then
            write(21, "(2(I0))") nn, mm
        end if
        
        if ((present(lineNr) .eqv. .true. ) .and. ( lineNr .eqv. .true.)) then
            write(form, "(A3, I0, A)") "(I,", (nn), "(I15, 2X))"
            do ii = 1, nn
            write(21, form) ii, (array(ii,jj), jj=1, mm)
            end do
        else
            write(form, "(A1, I0, A)") "(", (nn), "(I15, 2X))"
            do ii = 1, nn
            write(21, *) (array(ii,jj), jj=1, mm)
            end do
        end if
        
        close(21)
    end subroutine write_matrix_int

    !> Writes a (real) matrix into a file
    !!
    !! \param data          file to write in
    !! \param matrix        matrix which is to write
    !! \param lineNr        should the line number be written in the first column?
    !! \param size_out      should the size of the matrix be written in the first line of the file
    !! \param state_in      write in new or old file or replace existing file
    !! \param pos_in        write from begin of the file or append to it
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
        
        open(22, file=data, status=state, form="formatted", action="write", position=pos)
        
        if ((present(size_out) .eqv. .true.) .and. (size_out .eqv. .true. ))then
            write(22, *) nn, mm
        end if

        if ((present(lineNr) .eqv. .true. ) .and. ( lineNr .eqv. .true.)) then
            write(form, "(A3, I0, A)") "(I,", (nn), "(ES15.6, 2X))"
            do ii = 1, nn
            write(22, form) ii, (array(ii,jj), jj=1, mm)
            end do
        else 
            write(form, "(A1, I0, A)") "(", (mm), "(ES15.6, 2X))"
            do ii = 1, nn
            write(22, form) (array(ii,jj), jj=1, mm)
            end do
        end if
        
        close(22)
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

        open(23, file=data, status="replace", form="formatted", action="write")
        do kk = 1, n3
             do ii = 1, n1
                    write(23, *) (array(ii,jj,kk), jj=1, n2)
             end do
             write(23, *)
        end do

        close(23)
    end subroutine write_array_real

    !> Writing data by the line
    !!
    !! \param vec           Containing the data to write out
    !! \param data          Name of the outputfile
    !! \param size_out      should the size of the matrix be written in the first line of the file
    !! \param state_in      write in new or old file or replace existing file
    !! \param pos_in        write from begin of the file or append to it
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
        
        open(24, file=data, status=state, form="formatted", action="write", position=pos)
        if ((present(size_out)) .and. (size_out .eqv. .true.)) then
            write(24, "(I0)") size(vec)
        end if
        
        do ii = 1, size(vec)
             write(24, "(ES30.8)") vec(ii)
        end do
        close(24)
    end subroutine write_x
    
    !> Writing a (integer) vector by the line
    !!
    !! \param vec    Containing the data to write out
    !! \param data    Name of the outputfile
    !!
    subroutine write_x_int(data, vec)
        integer, intent(in) :: vec(:)
        character(*), intent(in) :: data
        integer :: ii

        open(25, file=data, status="replace", form="formatted", action="write")
        do ii = 1, size(vec)
             write(25, "(I6)") vec(ii)
        end do
        close(25)
    end subroutine write_x_int

    !> Writes data in XY-format
    !!
    !! \param xx    Containing x-data
    !! \param yy    Containing y-data
    !! \param data  Name of the outputfile
    !!
    subroutine write_xy(data, xx, yy)
        real(dp), intent(in) :: xx(:), yy(:)
        character(*), intent(in) :: data
        integer :: ii

        open(26, file=data, status="replace", form="formatted", action="write")
        do ii = 1, size(xx)
             write(26, "(2(ES20.12))") xx(ii), yy(ii)
        end do
        close(26)

    end subroutine write_xy


    !> Writing data in NXY-format
    !!
    !! \param xx     Containing x-data
    !! \param yy     Containing y-data column-wise
    !! \param data   Name of the outputfile
    !!
    subroutine write_nxy(data, xx, yy)
        real(dp), intent(in) :: xx(:), yy(:,:)
        character(*), intent(in) :: data
        integer :: ii, nx, ny
        character(25) :: form

        nx = size(xx)
        ny = size(yy(1,:))
        write(form, "(A1, I0, A)") "(", (ny+1), "(ES20.12, 2X))"

        open(27, file=data, status="replace", form="formatted", action="write")
        do ii = 1, nx
             write(27, form) xx(ii), yy(ii, 1:ny)
        end do
        close(27)

    end subroutine write_nxy

end module io