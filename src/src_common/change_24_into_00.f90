!***************************************************************************
! change_24_into_00.f90
! ---------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Change 24:00 into 00:00 of day after
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Change24Into00(Datestring, doy_format)
    use m_common_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: doy_format
    character(*), intent(inout) :: Datestring
    !> local variables
    integer :: int_year
    integer :: int_doy
    integer :: max_days
    integer :: i
    logical :: isleap
    character(5) :: tmp_date


    !> Determine whether current year is leap
    read(Datestring(1:4), '(i4)') int_year
    isleap = leapyear(int_year)
    if (isleap) then
        max_days = 366
    else
        max_days = 365
    end if

    if (doy_format) then
        if (isleap) then
            if (Datestring(10:11) == '24') then
                Datestring(10:11) = '00'
                do i = 1, lndays
                    tmp_date = Datestring(5:6) // '-' // Datestring(7:8)
                    if (DayOfLeapYear(i) == tmp_date) then
                        if (i /= lndays) then
                            Datestring(5:6) = DayOfLeapYear(i + 1)(1:2)
                            Datestring(7:8) = DayOfLeapYear(i + 1)(4:5)
                        else
                            Datestring(5:6) = DayOfLeapYear(1)(1:2)
                            Datestring(7:8) = DayOfLeapYear(1)(4:5)
                            int_year = int_year + 1
                            call int2char(int_year, Datestring(1:4), 0)
                        end if
                        exit
                    end if
                end do
            end if
        else
            if (Datestring(10:11) == '24') then
                Datestring(10:11) = '00'
                do i = 1, ndays
                    tmp_date = Datestring(5:6) // '-' // Datestring(7:8)
                    if (DayOfYear(i) == tmp_date) then
                        if (i /= ndays) then
                            Datestring(5:6) = DayOfYear(i + 1)(1:2)
                            Datestring(7:8) = DayOfYear(i + 1)(4:5)
                        else
                            Datestring(5:6) = DayOfYear(1)(1:2)
                            Datestring(7:8) = DayOfYear(1)(4:5)
                            int_year = int_year + 1
                            call int2char(int_year, Datestring(1:4), 0)
                        end if
                        exit
                    end if
                end do
            end if
        end if
    else
        if (Datestring(10:11) == '24') then
            Datestring(10:11) = '00'
            read(Datestring(5:7), '(i3)') int_doy
            int_doy = int_doy + 1
            if (int_doy == max_days + 1) then
                int_doy = 1
                int_year = int_year + 1
                call int2char(int_year, Datestring(1:4), 0)
            end if
            call int2char(int_doy, Datestring(5:7), 0)
        end if
    end if
end subroutine Change24Into00
