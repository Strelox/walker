!> Provides routines to read a Configuration file
!!
module config
    use accuracy
    implicit none

contains

    !> Reads Configuration of the random Walk
    !!
    !! \param data              File containing the configuration.
    !! \param maxTimestep       The maximum amount of timesteps that will be simulated.
    !! \param walk_mode         Mode of the simulated walk (CRW, QW or QSW).
    !! \param time_mode         The time mode of the simulation (either discrete or continous).
    !! \param io_rates          Input/Output rates.
    !! \param timestep          The size of a timestep.
    !! \param simulation_mode   New simulation or append to old data.
    !!
    subroutine read_config(data, maxTimestep, walk_mode, time_mode, io_rates, timestep, simulation_mode, write_mode, correction)
        character(*), intent(in) :: data
        integer, intent(out) :: maxTimestep
        real(dp), intent(out) :: timestep, correction
        real(dp), intent(out) :: io_rates(2)
        character(*), intent(out) :: walk_mode, time_mode, simulation_mode, write_mode
        integer :: status

        open(14, file=data, status="old", form="formatted", action="read")
        
        !! Read max Timestep
        read(14,*, iostat=status) maxTimestep
        if (status /= 0) then
             write(*,*) "Error: Wrong maxTimestep steps format in config file! End program."
             stop
        else if (maxTimestep < 0) then
            write(*,*) "Error: max. Timestep has to be non-negative in config file! End program."
            stop
        end if

        !! Read walk mode
        read(14,*, iostat=status) walk_mode
        if (status /= 0) then
             write(*,*) "Error: Wrong walk_mode format in config file! End program."
             stop
        else if ((walk_mode /= "CRW") .and. (walk_mode /= "QW")) then
            write(*,*) "Error: Unknown walk mode in config file! End program."
            stop
        end if

        !! Read time mode
        read(14,*, iostat=status) time_mode
        if (status /= 0) then
             write(*,*) "Error: Wrong time_mode format in config file! End program."
             stop
        else if ((time_mode /= "discrete") .and. (time_mode /= "continous")) then
            write(*,*) "Error: Unknown time mode in config file! End program."
            stop
        end if
        
        !! Read input/output rates
        read(14,*, iostat=status) io_rates(1), io_rates(2)
        if (status /= 0) then
             write(*,*) "Error: Wrong input/output rates format (io_rates) in config file! End program."
             stop
        end if

        !! Read timestep
        read(14,*, iostat=status) timestep
        if (status /= 0) then
             write(*,*) "Error: Wrong timestep format in config file! End program."
             stop
        else if (timestep < 0) then
            write(*,*) "Error: timestep has to be non-negative in config file! End program."
            stop
        end if
        
        !! Read simulation mode
        read(14,*, iostat=status) simulation_mode
        if (status /= 0) then
             write(*,*) "Error: Wrong simulation_mode format in config file! End program."
             stop
        else if ((simulation_mode /= "standard") .and. (simulation_mode /= "append")) then
            write(*,*) "Error: Unknown simulation mode in config file! End program."
            stop
        end if
        
        !! Read write mode
        read(14,*, iostat=status) write_mode
        if (status /= 0) then
             write(*,*) "Error: Wrong write_mode format in config file! End program."
             stop
        else if ((write_mode /= "all") .and. (write_mode /= "end")) then
            write(*,*) "Error: Unknown write mode in config file! End program."
            stop
        end if
        
        !! Read correction
        read(14,*, iostat=status) correction
        if (status /= 0) then
             write(*,*) "Error: Wrong correction format in config file! End program."
             stop
        else if ( (correction < 0.0_dp) .or. (correction > 1.0_dp)) then
            write(*,*) "Error: Dummy correction in config file! End program."
            stop
        end if
        
        close(14)
    end subroutine read_config

end module config
