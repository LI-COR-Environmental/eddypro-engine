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
    character(15) :: months(12)
    data months(1:12) / 'January', 'February', 'March', &
        'April', 'May', 'June', 'July', 'August', &
        'September', 'October', 'November', 'December' /


    write(lyear, '(i4)') date%year
    write(cday, '(i2)') date%day
    write(*, '(a)')
    write(*, '(a)', advance = 'no') init_message // '  ' // cday // ' ' // &
        trim(adjustl(months(date%month))) // ', ' // lyear // ' '

end subroutine ShowDailyAdvancement
