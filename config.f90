!     Copyright 2014 Frank Stuckenberg
!
!     This file is part of walker.
! 
!     walker is free software: you can redistribute it and/or modify
!     it under the terms of the GNU Affero General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     walker is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU Affero General Public License for more details.
! 
!     You should have received a copy of the GNU Affero General Public License
!     along with walker.  If not, see <http://www.gnu.org/licenses/>.
! 
!     Diese Datei ist Teil von walker.
! 
!     walker ist Freie Software: Sie können es unter den Bedingungen
!     der GNU Affero General Public License, wie von der Free Software Foundation,
!     Version 3 der Lizenz oder (nach Ihrer Wahl) jeder späteren
!     veröffentlichten Version, weiterverbreiten und/oder modifizieren.
! 
!     walker wird in der Hoffnung, dass es nützlich sein wird, aber
!     OHNE JEDE GEWÄHELEISTUNG, bereitgestellt; sogar ohne die implizite
!     Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
!     Siehe die GNU Affero General Public License für weitere Details.
! 
!     Sie sollten eine Kopie der GNU Affero General Public License zusammen mit diesem
!     Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.

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
    !! \param io_rates          Input/Output rates.
    !! \param timestep          The size of a timestep.
    !! \param simulation_mode   New simulation or append to old data.
    !!
    subroutine read_config(data, maxTimestep, walk_mode, io_rates, timestep, simulation_mode, write_mode, dec_time)
        character(*), intent(in) :: data
        integer, intent(out) :: maxTimestep
        real(dp), intent(out) :: timestep, dec_time
        real(dp), intent(out) :: io_rates(2)
        character(*), intent(out) :: walk_mode, simulation_mode, write_mode
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
        else if ((walk_mode /= "CRW") .and. (walk_mode /= "QW") .and. (walk_mode /= "QSW")) then
            write(*,*) "Error: Unknown walk mode in config file! End program."
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
        
        !! Read dec_time
        read(14,*, iostat=status) dec_time
        if (status /= 0) then
             write(*,*) "Error: Wrong dec_time format in config file! End program."
             stop
        end if
        
        close(14)
    end subroutine read_config

end module config
