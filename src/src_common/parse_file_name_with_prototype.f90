!*******************************************************************************
! parse_file_name_with_template.f90
! ---------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
!*******************************************************************************
!
! \brief       Retrieve timestamp string from file names based on Template
!              input date patterns: yyyy, yy, ddd, dd, mm
!              input time patterns: HH MM
!              normal output pattern:    yyyymmdd-HHMM
!              output pattern with DOY : yyyyddd--HHMM
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine ParseFileNameWithTemplate(Filename, Template, DateString)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: Filename
    character(*), intent(in) :: Template
    character(*), intent(out) :: DateString
    !> local variables
    character(4) :: loc_year
    character(2) :: loc_month
    character(2) :: loc_day
    character(3) :: loc_doy
    character(2) :: loc_hour
    character(2) :: loc_minute
    integer :: start


    !> year
    if (index(Template, 'yyyy') /= 0) then
        start = index(Template, 'yyyy')
        loc_year(1:4) = Filename(start: start + 3)
    else
        if (index(Template, 'yy') /= 0) then
            start = index(Template, 'yy')
            if (Filename(start: start + 1) > '70') then
                loc_year(1:4) = '19' // Filename(start: start + 1)
            else
                loc_year(1:4) = '20' // Filename(start: start + 1)
            end if
        else
            loc_year(1:4) = 'xxxx'
        end if
    end if

    !> month
    if (index(Template, 'mm') /= 0) then
        start = index(Template, 'mm')
        loc_month(1:2) = Filename(start : start + 3)
    else
        loc_month(1:2) = 'xx'
    end if
    !> day or DOY
    if (index(Template, 'ddd') /= 0) then
        start = index(Template, 'ddd')
        loc_doy(1:3) = Filename(start : start + 2)
        loc_day(1:2) = 'xx'
    else
        if (index(Template, 'dd') /= 0) then
            start = index(Template, 'dd')
            loc_day(1:2) = Filename(start : start + 1)
            loc_doy(1:3) = 'xxx'
        else
            loc_day(1:2) = 'xx'
            loc_doy(1:2) = 'xxx'
        end if
    end if

    !> hour
    if (index(Template, 'HH') /= 0) then
        start = index(Template, 'HH')
        loc_hour(1:2) = Filename(start : start + 1)
    else
        loc_hour(1:2) = 'xx'
    end if

    !> minute
    if (index(Template, 'MM') /= 0) then
        start = index(Template, 'MM')
        loc_minute(1:2) = Filename(start : start + 1)
    else
        loc_minute(1:2) = 'xx'
    end if

    !> Build up date string
    if (loc_doy(1:3) /= 'xxx') then
        Datestring = loc_year(1:4) // loc_doy(1:3) // '--' &
        // loc_hour(1:2) // loc_minute(1:2)
    elseif (loc_day(1:2) /= 'xx') then
        Datestring = loc_year(1:4) // loc_month(1:2) // loc_day(1:2) // '-' &
            // loc_hour(1:2) // loc_minute(1:2)
    end if
end subroutine ParseFileNameWithTemplate

