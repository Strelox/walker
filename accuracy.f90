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

!> Sets accuracy of calculation
!!
!! \param dp  sets double precision
!! \param err  sets lower accuracy bound
!!
module accuracy
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 300)
!   integer, parameter :: dp = selected_real_kind(33,4931)
    real(dp), parameter :: err = 1e-20

end module accuracy
