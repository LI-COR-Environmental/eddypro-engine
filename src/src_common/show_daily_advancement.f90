!***************************************************************************
! show_daily_advancement.f90
! --------------------------
! Copyright (C) 2011-2014, LI-COR Biosciences
!
! This file is part of EddyPro (TM).
!
! EddyPro (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyPro (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!
!***************************************************************************
!
! \brief	   Prints out date (up to the day) of currently processed period
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ShowDailyAdvancement(init_message, date)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*) :: init_message
    type (DateType), intent(in) :: date
    !> local variables
    character(4) :: lyear
    character(2) :: cday


    write(lyear, '(i4)') date%year
    write(cday, '(i2)') date%day
    select case (date%month)
        case(1)
            write(*, *) init_message // '  ' // cday // ' January ' // lyear
        case(2)
            write(*, *) init_message // '  ' // cday // ' February ' // lyear
        case(3)
            write(*, *) init_message // '  ' // cday // ' March ' // lyear
        case(4)
            write(*, *) init_message // '  ' // cday // ' April ' // lyear
        case(5)
            write(*, *) init_message // '  ' // cday // ' May ' // lyear
        case(6)
            write(*, *) init_message // '  ' // cday // ' June ' // lyear
        case(7)
            write(*, *) init_message // '  ' // cday // ' July ' // lyear
        case(8)
            write(*, *) init_message // '  ' // cday // ' August ' // lyear
        case(9)
            write(*, *) init_message // '  ' // cday // ' September ' // lyear
        case(10)
            write(*, *) init_message // '  ' // cday // ' October ' // lyear
        case(11)
            write(*, *) init_message // '  ' // cday // ' November ' // lyear
        case(12)
            write(*, *) init_message // '  ' // cday // ' December ' // lyear
    end select
end subroutine ShowDailyAdvancement
