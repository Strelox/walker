!> Provides routines to read a Configuration file
!!
module config
  use accuracy
  implicit none

contains

  !> Reads Configuration of the random Walk
  !!
  !! \param data        File containing the configuration.
  !! \param nWalker     Amount of walkers going through the network
  !! \param maxTimestep The maximum amount of timesteps that will be simulated.
  !! \param walk_mode   Mode of the simulated walk (CRW, QW or QSW).
  !! \param time_mode   The time mode of the simulation (either discrete or continous).
  !! \param start       Vortex where the walk will start.
  !! \param timestep    The size of a timestep.
  !!
  subroutine read_config(data, nWalker, maxTimestep, walk_mode, time_mode, start, timestep)
    character(*), intent(in) :: data
    integer, intent(out) :: nWalker, maxTimestep, start
    real(dp), intent(out) :: timestep
    character(*), intent(out) :: walk_mode, time_mode
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

    read(14,*, iostat=status) walk_mode
    if (status /= 0) then
       write(*,*) "Error: Wrong walk_mode format in config file! End program."
       stop
    end if

    read(14,*, iostat=status) time_mode
    if (status /= 0) then
       write(*,*) "Error: Wrong time_mode format in config file! End program."
       stop
    end if

    read(14,*, iostat=status) start
    if (status /= 0) then
       write(*,*) "Error: Wrong start format in config file! End program."
       stop
    end if

    if (start <= 0) then
       write(*,*) "Error: Start vortex is negative or zero. End program."
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
