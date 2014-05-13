module config
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

end module config