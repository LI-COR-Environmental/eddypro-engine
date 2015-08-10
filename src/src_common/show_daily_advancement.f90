!***************************************************************************
! show_daily_advancement.f90
! --------------------------
! Copyright (C) 2011-2015, LI-COR Biosciences
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
subroutine DisplayProgress(progress_type, init_message, tstamp, adv)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*) :: init_message
    character(*) :: progress_type
    character(*) :: adv
    type (DateType), intent(in) :: tstamp
    !> local variables
    character(10) :: date
    character(5) :: time
    character(15) :: months(12)
    data months(1:12) / 'January', 'February', 'March', &
        'April', 'May', 'June', 'July', 'August', &
        'September', 'October', 'November', 'December' /


    call DateTypeToDateTime(tstamp, date, time)
    select case(trim(adjustl(progress_type)))
    case ('daily')
        write(*, '(a)', advance = adv) init_message // date(9:10) // ' ' // &
            trim(adjustl(months(tstamp%month))) // ' ' // date(1:4) // ' '

    case ('avrg_interval')
        write(*, '(a)', advance = adv) init_message // time(1:5)
    end select

end subroutine DisplayProgress
